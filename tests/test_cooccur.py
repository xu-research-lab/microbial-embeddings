"""Tests for membed.cooccur_embedding.

The workflow expectations below are hand-computed from tests/data/test_raw.biom:

    counts = [[0, 10, 20],   # f0, absent in s0
              [10, 0, 40],   # f1, absent in s1
              [15, 20, 10],  # f2
              [40, 30, 0],   # f3, absent in s2
              [35, 40, 30]]  # f4

Each sample sums to 100, so the 'totalsum' scaling divides by 100, and the
non-zero values within a sample are all distinct, so within-sample ranks have
no ties. `cooccur_workflow` multiplies every similarity by n_samples (3)
before writing, so the stored value is the raw score times 3.
"""

import numpy as np
import pytest
from numpy import testing as npt

from membed.cooccur_embedding import (
    SUPPORTED_METRICS,
    RANK_DISTANCE_EPS,
    read_biom,
    gen_even_pairs,
    cooccur_workflow,
    build_x_max_file_workflow,
    get_feature_dict,
    cooccur_jaccard_dense,
    cooccur_braycurtis_similarity_dense,
    cooccur_abundance_dense,
    cooccur_faith_dense,
    cooccur_weighted_russell_rao_dense,
)

COUNTS = np.array([[0, 10, 20],
                   [10, 0, 40],
                   [15, 20, 10],
                   [40, 30, 0],
                   [35, 40, 30]], dtype=float)

# Stored values per metric, keyed by 0-based lower-triangular pair (i, j), i > j.
_EXPECTED = {
    'russell_rao': {(1, 0): 1, (2, 0): 2, (3, 0): 1, (4, 0): 2,
                    (2, 1): 2, (3, 1): 1, (4, 1): 2,
                    (3, 2): 2, (4, 2): 3, (4, 3): 2},
    'faith': {(1, 0): 1, (2, 0): 2, (3, 0): 1, (4, 0): 2,
              (2, 1): 2, (3, 1): 1, (4, 1): 2,
              (3, 2): 2, (4, 2): 3, (4, 3): 2},
    'jaccard': {(1, 0): 1, (2, 0): 2, (3, 0): 1, (4, 0): 2,
                (2, 1): 2, (3, 1): 1, (4, 1): 2,
                (3, 2): 2, (4, 2): 3, (4, 3): 2},
    'abundance_percentile': {(1, 0): 1 / 4, (2, 0): 3 / 8, (3, 0): 1 / 8,
                             (4, 0): 7 / 16, (2, 1): 1 / 4, (3, 1): 1 / 16,
                             (4, 1): 11 / 16, (3, 2): 5 / 8, (4, 2): 3 / 4,
                             (4, 3): 9 / 8},
    'abundance_totalsum': {(1, 0): 0.16, (2, 0): 0.18, (3, 0): 0.08,
                           (4, 0): 0.25, (2, 1): 0.165, (3, 1): 0.07,
                           (4, 1): 0.345, (3, 2): 0.2925, (4, 2): 0.36,
                           (4, 3): 0.6025},
    'braycurtis_percentile': {(1, 0): 3 / 4, (2, 0): 3 / 4, (3, 0): 3 / 10,
                              (4, 0): 9 / 13, (2, 1): 3 / 5, (3, 1): 1 / 4,
                              (4, 1): 4 / 5, (3, 2): 1, (4, 2): 1,
                              (4, 3): 18 / 17},
    'braycurtis_totalsum': {(1, 0): 3 / 4, (2, 0): 4 / 5, (3, 0): 3 / 10,
                            (4, 0): 2 / 3, (2, 1): 12 / 19, (3, 1): 1 / 4,
                            (4, 1): 24 / 31, (3, 2): 21 / 23, (4, 2): 9 / 10,
                            (4, 3): 39 / 35},
    'weighted_russell_rao': {
        (1, 0): 1 / (0.5 + RANK_DISTANCE_EPS),
        (2, 0): 2 / (0.25 + RANK_DISTANCE_EPS),
        (3, 0): 1 / (0.5 + RANK_DISTANCE_EPS),
        (4, 0): 1 / (0.75 + RANK_DISTANCE_EPS) + 1 / (0.25 + RANK_DISTANCE_EPS),
        (2, 1): 1 / (0.25 + RANK_DISTANCE_EPS) + 1 / (0.75 + RANK_DISTANCE_EPS),
        (3, 1): 1 / (0.75 + RANK_DISTANCE_EPS),
        (4, 1): 1 / (0.5 + RANK_DISTANCE_EPS) + 1 / (0.25 + RANK_DISTANCE_EPS),
        (3, 2): 1 / (0.5 + RANK_DISTANCE_EPS) + 1 / (0.25 + RANK_DISTANCE_EPS),
        (4, 2): 1 / (0.25 + RANK_DISTANCE_EPS) + 2 / (0.5 + RANK_DISTANCE_EPS),
        (4, 3): 2 / (0.25 + RANK_DISTANCE_EPS)},
}

COOCCUR_DTYPE = np.dtype([('otu1', 'i4'), ('otu2', 'i4'), ('value', 'f8')])


def _read_cooccur(cooccur_file):
    data = np.fromfile(cooccur_file, dtype=COOCCUR_DTYPE)
    return {(int(r['otu1']), int(r['otu2'])): float(r['value']) for r in data}


def _gen_even_pairs_helper(d, n):
    x = np.zeros(shape=(d, d), dtype='int')
    for beg, end, i, t in gen_even_pairs(d, n):
        for row in range(beg, end):
            for col in range(row):
                x[row, col] = i
    return x


def test_read_biom(biom_table):
    sid, fid, table = biom_table('test_raw.biom')
    assert list(sid) == ['s0', 's1', 's2']
    assert list(fid) == ['f0', 'f1', 'f2', 'f3', 'f4']
    npt.assert_array_equal(table.matrix_data.toarray(), COUNTS)


def test_read_biom_invalid_file(tmp_path):
    with pytest.raises(ValueError):
        read_biom(tmp_path / 'does-not-exist.biom')


@pytest.mark.parametrize(
    'd, n, exp',
    [(5,
      1,
      [[0, 0, 0, 0, 0],
       [1, 0, 0, 0, 0],
       [1, 1, 0, 0, 0],
       [1, 1, 1, 0, 0],
       [1, 1, 1, 1, 0]]),
     (5,
      2,
      [[0, 0, 0, 0, 0],
       [1, 0, 0, 0, 0],
       [1, 1, 0, 0, 0],
       [1, 1, 1, 0, 0],
       [2, 2, 2, 2, 0]]),
     (5,
      3,
      [[0, 0, 0, 0, 0],
       [1, 0, 0, 0, 0],
       [1, 1, 0, 0, 0],
       [2, 2, 2, 0, 0],
       [3, 3, 3, 3, 0]])])
def test_gen_even_pairs(d, n, exp):
    obs = _gen_even_pairs_helper(d, n)
    npt.assert_array_equal(obs, exp)


@pytest.mark.parametrize(
    'd, n',
    [(99, 2),
     (5000, 5),
     (9000, 9)])
def test_gen_even_pairs_complex(d, n):
    obs = _gen_even_pairs_helper(d, n)
    total_pairs = (d**2 - d) / 2
    pack_pairs = total_pairs / n
    # check all combinatorial pairs are accessed
    assert np.sum(obs > 0) == total_pairs
    pack_sizes = np.array([np.sum(obs == i + 1) for i in range(n)])
    deviates = pack_sizes / pack_pairs - 1
    # check all pack are more or less equal sizes - no deviates larger than 10%
    assert np.all(deviates < 0.1)


@pytest.mark.parametrize(
    'u, v, exp',
    [([0, 0, 1, 0, 1],
      [1, 1, 1, 0, 1],
      1 / 2),
     ([0, 0, 1, 0, 0],
      [1, 1, 1, 0, 0],
      1 / 3),
     ([1, 1, 1, 1, 1],
      [1, 1, 1, 1, 1],
      1.00),
     ([1, 1, 1, 1, 1],
      [0, 0, 0, 0, 0],
      0.00)])
def test_cooccur_jaccard_dense(u, v, exp):
    u = np.array(u)
    v = np.array(v)
    assert cooccur_jaccard_dense(u, v) == pytest.approx(exp)


@pytest.mark.parametrize(
    'u, v, exp',
    [([0, 0, 1, 0, 1],
      [1, 1, 1, 0, 1],
      (2 + 1 / 5) / 5),
     ([0, 0, 1, 0, 0],
      [1, 1, 1, 0, 0],
      (1 + 2 / 5) / 5),
     ([1, 1, 1, 1, 1],
      [1, 1, 1, 1, 1],
      1.00),
     ([1, 1, 1, 1, 1],
      [0, 0, 0, 0, 0],
      0.00)])
def test_cooccur_faith_dense(u, v, exp):
    u = np.array(u)
    v = np.array(v)
    assert cooccur_faith_dense(u, v) == pytest.approx(exp)


@pytest.mark.parametrize(
    'u, v, exp',
    [([0.10, 0.00, 0.20, 0.00, 0.30],
      [0.00, 0.10, 0.20, 0.10, 0.10],
      3 / 11),
     ([0.10, 0.10, 0.20, 0.10, 0.30],
      [0.10, 0.10, 0.20, 0.10, 0.30],
      0.5),
     ([0.10, 0.10, 0.20, 0.10, 0.30],
      [0.00, 0.00, 0.00, 0.00, 0.00],
      0.00)])
def test_cooccur_braycurtis_similarity_dense(u, v, exp):
    u = np.array(u)
    v = np.array(v)
    assert cooccur_braycurtis_similarity_dense(u, v) == pytest.approx(exp)


@pytest.mark.parametrize(
    'u, v, exp',
    [([0.10, 0.00, 0.20, 0.00, 0.30],
      [0.00, 0.10, 0.20, 0.10, 0.10],
      0.056),
     ([0.10, 0.10, 0.20, 0.10, 0.30],
      [0.10, 0.10, 0.20, 0.10, 0.30],
      0.16),
     ([0.10, 0.10, 0.20, 0.10, 0.30],
      [0.00, 0.00, 0.00, 0.00, 0.00],
      0.00)])
def test_cooccur_abundance_dense(u, v, exp):
    u = np.array(u)
    v = np.array(v)
    assert cooccur_abundance_dense(u, v) == pytest.approx(exp)


@pytest.mark.parametrize(
    'u, v, exp',
    [([0, 0.5, 0, 1.0],
      [0.5, 0.5, 1.0, 0.0],
      1 / (4 * RANK_DISTANCE_EPS)),
     ([0, 0, 0, 0],
      [1, 1, 1, 1],
      0.00),
     ([0.25, 0.75, 0, 1.0],
      [0.75, 0.75, 1.0, 0.5],
      (2 / (0.5 + RANK_DISTANCE_EPS) + 1 / RANK_DISTANCE_EPS) / 4),
     ([0.5, 0.5, 0, 0.5],
      [0.5, 0.5, 0, 0.5],
      3 / (4 * RANK_DISTANCE_EPS))])
def test_cooccur_weighted_russell_rao_dense(u, v, exp):
    u = np.array(u, dtype=float)
    v = np.array(v, dtype=float)
    assert cooccur_weighted_russell_rao_dense(u, v) == pytest.approx(exp)


@pytest.mark.parametrize('metric', SUPPORTED_METRICS)
def test_cooccur_workflow(data_dir, tmp_path, metric):
    cooccur_file = tmp_path / f'{metric}.co'
    cooccur_workflow(data_dir / 'test_raw.biom', cooccur_file,
                     metric=metric, cpus=1)

    obs = _read_cooccur(cooccur_file)
    exp = _EXPECTED[metric]

    # Symmetric expansion: both (i, j) and (j, i) for every non-zero pair,
    # with 1-based feature indices.
    assert len(obs) == 2 * len(exp)
    for (i, j), value in exp.items():
        assert obs[(i + 1, j + 1)] == pytest.approx(value, rel=1e-5)
        assert obs[(j + 1, i + 1)] == pytest.approx(value, rel=1e-5)


def test_cooccur_workflow_parallel(data_dir, tmp_path):
    cooccur_file = tmp_path / 'par.co'
    cooccur_workflow(data_dir / 'test_raw.biom', cooccur_file,
                     metric='jaccard', cpus=2)
    obs = _read_cooccur(cooccur_file)
    for (i, j), value in _EXPECTED['jaccard'].items():
        assert obs[(i + 1, j + 1)] == pytest.approx(value)


def test_cooccur_workflow_unknown_metric(data_dir, tmp_path):
    with pytest.raises(ValueError, match='Unknown metric'):
        cooccur_workflow(data_dir / 'test_raw.biom', tmp_path / 'x.co',
                         metric='bogus')


def test_get_feature_dict(data_dir, tmp_path):
    feature_dict = tmp_path / 'feature-dict.csv'
    get_feature_dict(data_dir / 'test_raw.biom', feature_dict)
    assert feature_dict.read_text().splitlines() == [
        'f0 2', 'f1 2', 'f2 3', 'f3 2', 'f4 3']


def test_build_x_max_file_workflow(tmp_path):
    records = np.array(
        [(1, 2, 1.0), (2, 1, 1.0),
         (1, 3, 2.0), (3, 1, 2.0),
         (2, 3, 4.0), (3, 2, 4.0),
         (1, 4, 8.0), (4, 1, 8.0)], dtype=COOCCUR_DTYPE)
    cooccur_file = tmp_path / 'table.co'
    records.tofile(cooccur_file)

    x_max_file = tmp_path / 'xmax'
    build_x_max_file_workflow(cooccur_file, x_max_file, percentile_num=50)
    # 50th percentile of [1, 1, 2, 2, 4, 4, 8, 8], and '.npy' is appended.
    assert np.load(f'{x_max_file}.npy') == pytest.approx(3.0)
