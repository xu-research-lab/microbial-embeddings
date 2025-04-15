from datetime import datetime

import pytest
import numpy as np
from numpy import testing as npt
from line_profiler import LineProfiler

from membed.glove import (gen_even_pairs, build_cooccur_matrix, cooccur_workflow,
                         cooccur_jaccard_dense,
                          cooccur_dice_dense, cooccur_braycurtis_dense,
                          cooccur_abundance_dense, cooccur_faith_dense,
                          cooccur_phi_dense)

@pytest.mark.skipif(
    "not config.getoption('--table')",
    reason="Only run time and memory profiling when --table is given")
def test_profile(table, dense, metric, cpus, lpo, mpo):
    functions = [cooccur_abundance_sparse, cooccur_abundance_dense, build_cooccur_matrix]
    lp = LineProfiler()
    for fun in functions:
        lp.add_function(fun)
    lp_wrapper = lp(cooccur_workflow)

    import tracemalloc
    tracemalloc.start(10)  # Save up to 10 stack frames
    before = tracemalloc.take_snapshot()
    _ = lp_wrapper(table, cpus=cpus, normalize=True, dense=dense, metric=metric)
    after = tracemalloc.take_snapshot()
    stats_line = after.compare_to(before, 'lineno')
    stats_trace = after.compare_to(before, 'traceback')

    ct = datetime.now().strftime("%Y%m%d%H%M%S")
    if mpo is None:
        mpo = f'{table}.{metric}.{dense}.{cpus}.{ct}.mpo'

    with open(mpo, 'w') as o:
        for line in stats_line[:50]:
            print(line, file=o)
        #print(f"\n*** Trace for 3 largest memory blocks - ({largest.count} blocks, {largest.size/1024} Kb) ***")
        for line in stats_trace[0].traceback.format():
            print(line, file=o)

    if lpo is None:
        lpo = f'{table}.{metric}.{dense}.{cpus}.{ct}.lpo'

    with open(lpo, 'w') as o:
        lp.print_stats(o)


def _gen_even_pairs_helper(d, n):
    x = np.zeros(shape=(d, d), dtype='int')
    for beg, end, i, t in gen_even_pairs(d, n):
        for row in range(beg, end):
            for col in range(row):
                x[row, col] = i
    return x


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
    pack_sizes = np.array([np.sum(obs==i+1) for i in range(n)])
    deviates = pack_sizes / pack_pairs - 1
    # check all pack are more or less equal sizes - no deviates larger than 10%
    assert np.all(deviates < 0.1)


@pytest.mark.parametrize(
    'normalize, percentile, dense, exp',
    [(False,
      False,
      True,
      np.array([[  0., 200.,   0.],
                [100.,   0., 300.],
                [200., 100., 100.],
                [550., 200.,   0.],
                [150., 500., 600.]])),
     (False,
      True,
      True,
      np.array([[0.00, 5/8, 0.00],
                [1/4, 0.00, 2/3],
                [3/4, 1/4, 1/3],
                [1.00, 5/8, 0.00],
                [1/2, 1.00, 1.00]])),
     (True,
      False,
      True,
      np.array([[0.00, 0.20, 0.00],
                [0.10, 0.00, 0.30],
                [0.20, 0.10, 0.10],
                [0.55, 0.20, 0.00],
                [0.15, 0.50, 0.60]]))
     ])
def test_read_biom(biom_table, normalize, percentile, dense, exp):
    obs = biom_table('test_raw.biom', normalize=normalize, percentile=percentile, dense=dense)[2]
    npt.assert_array_equal(obs, exp)


@pytest.mark.parametrize(
    'u, v, exp',
    [([0, 0, 1, 0, 1],
      [1, 1, 1, 0, 1],
      1/2),
     ([0, 0, 1, 0, 0],
      [1, 1, 1, 0, 0],
      1/3),
     ([1, 1, 1, 1, 1],
      [1, 1, 1, 1, 1],
      1.00),
     ([1, 1, 1, 1, 1],
      [0, 0, 0, 0, 0],
      0.00)
     ])
def test_cooccur_jaccard_dense(u, v, exp):
    u = np.array(u)
    v = np.array(v)
    assert cooccur_jaccard_dense(u, v) == pytest.approx(exp)


@pytest.mark.parametrize(
    'u, v, exp',
    [([0, 0, 1, 0, 1],
      [1, 1, 1, 0, 1],
      (2+1/5)/5),
     ([0, 0, 1, 0, 0],
      [1, 1, 1, 0, 0],
      (1+2/5)/5),
     ([1, 1, 1, 1, 1],
      [1, 1, 1, 1, 1],
      1.00),
     ([1, 1, 1, 1, 1],
      [0, 0, 0, 0, 0],
      0.00)
     ])
def test_cooccur_faith_dense(u, v, exp):
    u = np.array(u)
    v = np.array(v)
    assert cooccur_faith_dense(u, v) == pytest.approx(exp)


@pytest.mark.parametrize(
    'u, v, exp',
    [([0, 0, 1, 0, 1],
      [1, 1, 1, 0, 1],
      (2*1-0*2)/np.sqrt((2+0)*(2+2)*(0+1)*(2+1))),
     ([0, 0, 1, 0, 0],
      [1, 1, 1, 0, 0],
      (1*2-0*2)/np.sqrt((1+0)*(1+2)*(0+2)*(2+2))),
     ([1, 1, 0, 1, 0],
      [1, 0, 1, 0, 1],
      (1*0-2*2)/np.sqrt((1+2)*(1+2)*(2+0)*(2+0)))])
def test_cooccur_phi_dense(u, v, exp):
    u = np.array(u)
    v = np.array(v)
    assert cooccur_phi_dense(u, v) == pytest.approx(exp)


@pytest.mark.parametrize(
    'u, v, exp',
    [([0, 0, 1, 0, 1],
      [1, 1, 1, 0, 1],
      2/3),
     ([0, 0, 1, 0, 0],
      [1, 1, 1, 0, 0],
      0.50),
     ([1, 1, 1, 1, 1],
      [1, 1, 1, 1, 1],
      1.00),
     ([1, 1, 1, 1, 1],
      [0, 0, 0, 0, 0],
      0.00)
     ])
def test_cooccur_dice_dense(u, v, exp):
    u = np.array(u)
    v = np.array(v)
    assert cooccur_dice_dense(u, v) == pytest.approx(exp)


@pytest.mark.parametrize(
    'u, v, exp',
    [([0.10, 0.00, 0.20, 0.00, 0.30],
      [0.00, 0.10, 0.20, 0.10, 0.10],
      3/11),
     ([0.10, 0.10, 0.20, 0.10, 0.30],
      [0.10, 0.10, 0.20, 0.10, 0.30],
      0.5),
     ([0.10, 0.10, 0.20, 0.10, 0.30],
      [0.00, 0.00, 0.00, 0.00, 0.00],
      0.00)
     ])
def test_cooccur_braycurtis_dense(u, v, exp):
    u = np.array(u)
    v = np.array(v)
    assert cooccur_braycurtis_dense(u, v) == pytest.approx(exp)


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
      0.00)
     ])
def test_cooccur_abundance_dense(u, v, exp):
    u = np.array(u)
    v = np.array(v)
    assert cooccur_abundance_dense(u, v) == pytest.approx(exp)


@pytest.mark.parametrize(
    'normalize, percentile, dense, metric, cpus, exp',
    [(False,
      False,
      False,
      'binary',
      1,
      np.array([[0.00, 0.00, 0.00, 0.00, 0.00],
                [0.00, 0.00, 0.00, 0.00, 0.00],
                [1/3, 2/3, 0.00, 0.00, 0.00],
                [1/3, 1/3, 2/3, 0.00, 0.00],
                [1/3, 2/3, 1.00, 2/3, 0.00]])),
      (False,
       True,
       True,
       'abundance',
       1,
       np.array([[0.00, 0.00, 0.00, 0.00, 0.00],
                 [0.00, 0.00, 0.00, 0.00, 0.00],
                 [5/96, 25/216, 0.00, 0.00, 0.00],
                 [5/24, 1/48, 23/96, 0.00, 0.00],
                 [25/192, 91/432, 79/432, 41/192, 0.00]])),
      (True,
       False,
       True,
       'abundance',
       1,
       np.array([[0.00, 0.00, 0.00, 0.00, 0.00],
                 [0.00, 0.00, 0.00, 0.00, 0.00],
                 [0.03, 17/300, 0.00, 0.00, 0.00],
                 [2/30, 11/600, 11/150, 0.00, 0.00],
                 [7/150, 61/600, 101/1200, 23/300, 0.00]])),
      (False,
       False,
       True,
       'jaccard',
       1,
       np.array([[0.00, 0.00, 0.00, 0.00, 0.00],
                 [0.00, 0.00, 0.00, 0.00, 0.00],
                 [1/9, 2/9, 0.00, 0.00, 0.00],
                 [1/6, 1/9, 2/9, 0.00, 0.00],
                 [1/9, 2/9, 1/3, 2/9, 0.00]])),
      (True,
       False,
       True,
       'abundance',
       2,
       np.array([[0.00, 0.00, 0.00, 0.00, 0.00],
                 [0.00, 0.00, 0.00, 0.00, 0.00],
                 [0.03, 17/300, 0.00, 0.00, 0.00],
                 [2/30, 11/600, 11/150, 0.00, 0.00],
                 [7/150, 61/600, 101/1200, 23/300, 0.00]])),
      (False,
       False,
       True,
       'buildCooccuranceMatrix',
       2,
       np.array([[0.00, 0.00, 0.00, 1.00, 1.00],
                 [0.00, 0.00, 1.00, 0.00, 2.00],
                 [0.00, 1.00, 0.00, 2.00, 1.00],
                 [1.00, 0.00, 2.00, 0.00, 0.00],
                 [1.00, 2.00, 1.00, 0.00, 0.00]]))
     ])
def test_cooccur_workflow(data_dir, normalize, percentile, dense, metric, cpus, exp):
    cooccur = cooccur_workflow(data_dir / 'test_raw.biom', normalize=normalize,
                               percentile=percentile, dense=dense, metric=metric, cpus=cpus)
    npt.assert_array_almost_equal(cooccur.toarray(), exp)




