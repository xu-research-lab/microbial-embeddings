"""Co-occurrence based GloVe embeddings for microbiome feature tables.

This module turns a BIOM feature table into dense vector representations of
microbial features, reusing the GloVe word-embedding model with samples playing
the role of documents and features the role of words.

The pipeline has four stages, each exposed as a sub-command in `cli.py`:

1. `get_feature_dict` - write the feature vocabulary and its occurrence counts.
2. `cooccur_workflow` - score every feature pair with one of
   `SUPPORTED_METRICS` and write the co-occurrence records.
3. `build_x_max_file_workflow` - pick the x_max weighting cutoff from those
   records.
4. `train_glove_model` - shell out to the reference GloVe C binaries in
   `glove_build/` to fit the embeddings.

Stages 2-4 exchange co-occurrence records as a headerless binary array of dtype
[('otu1', 'i4'), ('otu2', 'i4'), ('value', 'f8')], the format the GloVe
binaries read. Feature indices in that file are 1-based, since GloVe reserves 0.

Note that BIOM stores tables as (observation, sample), so throughout this module
rows are features and columns are samples.
"""

import os
import sys
import subprocess
from math import ceil, sqrt
from logging import getLogger

import numpy as np
import pandas as pd

import biom
from scipy.sparse import dok_array, coo_array, tril

from joblib import Parallel, delayed
from tqdm import tqdm


logger = getLogger(__name__)

# Metric names accepted by `cooccur_workflow`.
SUPPORTED_METRICS = (
    'russell_rao',
    'weighted_russell_rao',
    'abundance_percentile',
    'abundance_totalsum',
    'braycurtis_percentile',
    'braycurtis_totalsum',
    'jaccard',
    'faith',
)


def read_biom(biom_file):
    """Read BIOM format file and return sample ids, feature ids and the table.

    Parameters:
        biom_file (str or file-like object): Path to the BIOM file or an already
            opened file handle.

    Returns:
        tuple: A 3-tuple of
            - sid (numpy.ndarray): Sample identifiers from the BIOM table.
            - fid (numpy.ndarray): Feature/observation identifiers (e.g. taxonomy
              IDs or gene names).
            - table (biom.Table): The loaded table itself. Densifying or
              normalising `table.matrix_data` is left to the caller, since that
              depends on the chosen metric.

    Raises:
        ValueError: If the file cannot be parsed as a BIOM table.
    """

    try:
        table = biom.load_table(biom_file)
    except Exception as e:
        raise ValueError(f"Failed to load BIOM data: {str(e)}") from e

    logger.info('Loaded %d features, %d samples' % table.shape)

    sid = table.ids(axis='sample')
    fid = table.ids(axis='observation')


    return sid, fid, table


def cooccur_jaccard_dense(u, v):
    """Calculate Jaccard similarity between two features.

    Parameters:
        u, v (numpy.ndarray): 1-D boolean arrays of length n_samples, giving the
            presence/absence of two features across samples.

    Returns:
        float: |u & v| / |u | v|, in [0, 1]. Higher means the two features are
            present in more of the same samples.
    """
    union = (u | v)
    intersect = (u & v)
    return intersect.sum() / union.sum() 


def cooccur_faith_dense(u, v):
    """Calculate Faith's phylogenetic similarity between two features.

    Parameters:
        u, v (numpy.ndarray): 1-D boolean arrays of length n_samples, giving the
            presence/absence of two features across samples.

    Returns:
        float: (a + d / n) / n, where a is the number of samples holding both
            features, d the number holding neither, and n = n_samples. The
            joint-absence term is damped by n rather than the textbook factor
            of 2, so shared presences dominate the score.
    """
    a = (u & v)
    d = (1 - u) & (1 - v)
    n = len(u)
    return np.sum(a + d / n) / n


def cooccur_braycurtis_similarity_dense(u, v):
    """Calculate Bray-Curtis *similarity* between two features.

    Parameters:
        u, v (numpy.ndarray): 1-D arrays of length n_samples holding the scaled
            abundance of two features across samples. `cooccur_workflow` supplies
            either within-sample ranks or relative abundances.

    Returns:
        float: sum(min(u, v)) / sum(u + v), in [0, 0.5]. Note this is half the
            usual Bray-Curtis similarity; the corresponding dissimilarity is
            1 - 2 * this value. Higher means more similar, which is the
            orientation the co-occurrence matrix expects.
    """
    m = u - v
    min_uv = u - np.maximum(m, 0)
    return min_uv.sum() / (u + v).sum()


def cooccur_abundance_dense(u, v):
    """Calculate co-occurrence strength weighted by abundance agreement.

    Parameters:
        u, v (numpy.ndarray): 1-D arrays of length n_samples holding the
            abundance of two features across samples, already scaled to [0, 1]
            (as `cooccur_workflow` does via rankdata/norm). Outside that range
            the (1 - |u - v|) factor turns negative.

    Returns:
        float: mean(min(u, v) * (1 - |u - v|)). A sample contributes most when
            both features are abundant *and* similarly abundant.
    """
    m = u - v
    min_uv = u - np.maximum(m, 0)
    return (np.sum(min_uv - np.abs(m) * min_uv))  / len(u) 
    # return (np.sum(min_uv - np.abs(m) * min_uv))


#: Guards 1/d against rank ties (d == 0), which are common in count data.
RANK_DISTANCE_EPS = 1e-4


def cooccur_weighted_russell_rao_dense(u, v):
    """Calculate weighted Russell-Rao similarity (Eq. 2).

    Parameters:
        u, v (numpy.ndarray): 1-D arrays of length n_samples holding the
            max-normalised within-sample rank of two features (0 = absent), so
            ranks lie in (0, 1] and their distance d lies in [0, 1).

    Returns:
        float: The sum of 1 / (d + RANK_DISTANCE_EPS) over co-occurring samples,
            divided by n_samples. Only samples where both features are present
            contribute; the closer their within-sample ranks, the more that
            sample counts, capped at 1 / RANK_DISTANCE_EPS on an exact tie.
    """
    n = len(u)
    co = (u > 0) & (v > 0)              # only co-occurring samples count
    d = np.abs(u - v)
    w = 1.0 / (d[co] + RANK_DISTANCE_EPS)
    return w.sum() / n


def gen_even_pairs(n, n_packs):
    """Divide element pairs in the lower triangular matrix (total n*(n-1)/2 pairs) into evenly distributed task chunks for parallel processing.

    Parameters:
        n (int): Total number of items (rows/columns in the matrix).
        n_packs (int): Number of chunks (tasks) to divide the pairs into.

    Yields:
        tuple: A tuple containing four elements:
            - beg (int): Starting row index of the chunk (0-based, inclusive).
            - end (int): Ending row index of the chunk (0-based, exclusive).
            - pack_number (int): Chunk identifier (1 to `n_packs`).
            - pairs (int): Number of pairs in the chunk.
    """
    # Total elements in lower triangular matrix
    total = n * (n - 1) // 2
    pack = total // n_packs
    beg = end = 0
    for i in range(1, n_packs + 1):
        if i == n_packs:
            # The per-chunk boundary is rounded up, so the running total can
            # drift; pin the final chunk to n or trailing rows are dropped.
            end = n
        else:
            end = (1 + sqrt(1 - 4 * (beg - beg**2 - 2 * pack))) / 2
            end = ceil(end)
            if end > n:
                end = n
        yield beg, end, i, (beg + end - 1) * (end - beg) // 2
        beg = end


def _x_max_path(x_max_file):
    """Normalise the x_max path so save and load always agree.

    `np.save` silently appends '.npy' when the name lacks it, so loading the
    original name would raise FileNotFoundError. Applying this on both sides
    keeps the two in sync whatever the caller passes.

    Parameters:
        x_max_file (str or path-like): Path as supplied by the caller, with or
            without suffix.

    Returns:
        str: The same path, guaranteed to end in '.npy'.
    """
    x_max_file = os.fspath(x_max_file)
    return x_max_file if x_max_file.endswith('.npy') else f'{x_max_file}.npy'


def build_x_max_file_workflow(cooccur_file, x_max_file, percentile_num=10):
    """Calculate and save x_max value for GloVe training from co-occurrence data.

    Parameters:
        cooccur_file (str):
            Input path of the co-occurrence records written by `cooccur_workflow`.
            Holds all non-zero pairs (otu1, otu2, value) with symmetric expansion
            (both (i, j) and (j, i) are stored), as a headerless binary array of
            dtype [('otu1', 'i4'), ('otu2', 'i4'), ('value', 'f8')].
        x_max_file (str):
            Output path for the x_max value; '.npy' is appended when missing.
            Stores the `percentile_num`-th percentile of co-occurrence values. In GloVe models:
            - The x_max hyperparameter caps the weighting function's upper bound
            - When co-occurrence value x exceeds x_max, weight w(x) is fixed at 1.0
            - This prevents high-frequency pairs from dominating training and mitigates outlier effects.
        percentile_num (float, optional):
            Percentile value (0-100) to use for x_max calculation. Default is 10 (10th percentile).
            Typical GloVe usage recommends 99th percentile (set percentile_num=99).

    Returns:
        None: The percentile is written to `x_max_file` as a 0-d .npy array.
    """

    dt = np.dtype([('otu1', 'i4'), ('otu2', 'i4'), ('value', 'f8')])
    data = np.fromfile(cooccur_file, dtype=dt)

    percentile_val = np.percentile(data['value'], percentile_num)
    np.save(_x_max_path(x_max_file), percentile_val)


def get_feature_dict(biom_file, feature_dict):
    """Count the number of non-zero occurrences for each feature across all samples in the dataset.

    Parameters:
        biom_file (str): Path to input BIOM file containing feature abundance data.
        feature_dict (str): Output file path to save the feature dictionary. File format:
            Two columns: [feature_id] [non_zero_count], separated by spaces.

    Returns:
        None: The vocabulary is written to `feature_dict`, in the order GloVe
            expects for its -vocab-file argument.
    """
    table = biom.load_table(biom_file)
    tmp = pd.DataFrame({
        'id': table.ids(axis='observation'),
        'value': table.nonzero_counts(axis='observation')
    })
    tmp.to_csv(feature_dict, sep=" ", index=False, header=None)


def build_cooccur_matrix(X, metric, cpus=1):
    """Build a sparse co-occurrence matrix of pairwise feature similarities, in parallel.

    Parameters:
        X (np.ndarray):
            A 2D array of shape (n_features, n_samples), i.e. one row per feature,
            matching the layout of `biom.Table.matrix_data`. Similarities are
            computed between row pairs, so the result is feature-by-feature.
        metric (callable):
            Function taking two 1D arrays (the rows of `X` for two features) and
            returning a float similarity. Higher must mean more co-occurring;
            values of exactly 0 are dropped to keep the result sparse.
        cpus (int, optional):
            Number of CPU cores to use for parallel processing. Defaults to 1 (sequential).

    Returns:
        cooccur (scipy.sparse.dok_array):
            A sparse matrix of shape (n_features, n_features) holding similarity
            scores. Only the lower-triangular part (excluding the diagonal) is
            filled to avoid redundancy.
    """

    logger.debug(f'Use {cpus} processes.')
    logger.debug(f'X data type: {type(X)}')
    n = X.shape[0]
    disable = logger.getEffectiveLevel() > 10

    def _fill_cooccur_mat(beg, end, i, pairs):
        logger.debug(f'Begin row {beg} and end row {end} have {pairs} pairs')
        cooccur_chunk = dok_array((n, n), dtype='float32')
        with tqdm(total=pairs, position=i, desc=f'CPU core {i}', disable=disable) as pbar:
            for row in range(beg, end):
                for col in range(0, row):
                    pbar.update()
                    v = metric(X[row], X[col])
                    if v != 0:
                        cooccur_chunk[row, col] = v
        return cooccur_chunk

    with Parallel(n_jobs=cpus) as parallel:
        chunks = parallel(
            delayed(_fill_cooccur_mat)(beg, end, i, pairs)
            for beg, end, i, pairs in gen_even_pairs(n, cpus))
        cooccur = np.sum(chunks)

    return cooccur


def cooccur_workflow(biom_file,
                     cooccur_file,
                     metric='russell_rao',
                     cpus=1,):
    """Workflow for co-occurrence analysis with multiple similarity metrics.

    Parameters:
        biom_file (str): Input BIOM file path.
        cooccur_file (str): Output path for co-occurrence data. Written as a binary
            packed array of dtype [('otu1', 'i4'), ('otu2', 'i4'), ('value', 'f8')],
            symmetrically expanded so both (i, j) and (j, i) are present. Defaults to
            '<biom_file>.cooccur' when None.
        metric (str): Similarity metric, one of `SUPPORTED_METRICS`. The
            '_percentile' variants scale each sample by its within-sample ranks,
            the '_totalsum' variants by relative abundance.
        cpus (int): Number of CPU cores for parallel processing. Ignored by
            'russell_rao', which is computed as a single matrix product.

    Returns:
        None: The co-occurrence records are written to `cooccur_file`. Feature
            indices are 1-based, since GloVe reserves index 0.

    Raises:
        ValueError: If `metric` is not one of `SUPPORTED_METRICS`.
    """

    sid, fid, table = read_biom(biom_file)
    
    logger.info('Compute cooccur matrix using %r...' % metric)

    if metric == 'russell_rao':
        data = table.matrix_data.toarray()
        data = data.astype(bool).astype('int')
        logger.debug('Done computing cooccurence.')

    elif metric == 'weighted_russell_rao':
        table.rankdata(axis='sample', inplace=True)
        data = table.matrix_data.multiply(1 / table.max(axis='sample')).toarray()
        metric = cooccur_weighted_russell_rao_dense
        logger.debug('Done computing cooccurence.')

    elif metric == 'abundance_percentile':
        table.rankdata(axis='sample', inplace=True)
        data = table.matrix_data.multiply(1 / table.max(axis='sample')).toarray()
        metric = cooccur_abundance_dense

    elif metric == 'abundance_totalsum':
        table.norm(axis='sample', inplace=True)
        data = table.matrix_data.toarray()
        metric = cooccur_abundance_dense

    elif metric == 'braycurtis_percentile':
        table.rankdata(axis='sample', inplace=True)
        data = table.matrix_data.multiply(1 / table.max(axis='sample')).toarray()
        metric = cooccur_braycurtis_similarity_dense

    elif metric == 'braycurtis_totalsum':
        table.norm(axis='sample', inplace=True)
        data = table.matrix_data.toarray()
        metric = cooccur_braycurtis_similarity_dense

    elif metric == 'jaccard':
        data = table.matrix_data.toarray()
        data = data.astype(bool)
        metric = cooccur_jaccard_dense

    elif metric == 'faith':
        data = table.matrix_data.toarray()
        data = data.astype(bool)
        metric = cooccur_faith_dense

    else:
        raise ValueError('Unknown metric %r. Supported metrics: %s'
                         % (metric, ', '.join(SUPPORTED_METRICS)))

    if metric == 'russell_rao':
        cooccur = tril(np.dot(data, data.T), k=-1) / len(sid)
    else:
        cooccur = build_cooccur_matrix(data, metric=metric, cpus=cpus)

    logger.debug('Cooccur matrix consumed mem:%.2fMB' %
                 (sys.getsizeof(cooccur) / (1024 * 1024)))
    
    if cooccur_file is None:
        cooccur_file = f'{biom_file}.cooccur'

    cooccur = cooccur * len(sid) 
    cooccur = coo_array(cooccur)

    dt = np.dtype([('otu1', 'int32'), ('otu2', 'int32'), ('value', 'float64')])
    row = np.append(cooccur.row, cooccur.col)
    col = np.append(cooccur.col, cooccur.row)
    value = np.append(cooccur.data, cooccur.data)

    temp = np.empty((len(row), ), dt)
    temp['otu1'] = row + 1
    temp['otu2'] = col + 1
    temp['value'] = value

    temp.tofile(cooccur_file)


def train_glove_model(cooccur_file,
                      x_max_file,
                      feature_dict,
                      result,
                      lr=0.05,
                      cpus=1,
                      n_iter=10,
                      embedding_size=100):
    """Train GloVe model using original C implementation.

    Parameters:
        cooccur_file (str):
            Path to co-occurrence data file generated by `cooccur_workflow`.
            Format: binary packed array of dtype
            [('otu1', 'i4'), ('otu2', 'i4'), ('value', 'f8')].
        x_max_file (str):
            Path to file containing the x_max value (float) for GloVe's weighting
            function, generated by `build_x_max_file_workflow`. A missing '.npy'
            suffix is added automatically.
        feature_dict (str):
            Path to feature dictionary mapping feature IDs to metadata.
            Format: space-separated [feature_id] [non_zero_count] (from `get_feature_dict`).
        result (str):
            Output *directory*, created if absent. Files written inside it:
            - Embeddings: embeddings_<embedding_size>.txt
            - Temporary file: embeddings_temp.shuf.bin
            - Gradients: embeddings_gradsq
        lr (float, optional):
            Learning rate (eta) for AdaGrad optimization. Defaults to 0.05.
        cpus (int, optional):
            Number of CPU cores for parallel processing. Defaults to 1.
        n_iter (int, optional):
            Number of training iterations over the dataset. Defaults to 10.
        embedding_size (int, optional):
            Dimension of output embeddings. Defaults to 100.

    Returns:
        None: The trained vectors are written under `result` by the GloVe binary.

    Raises:
        subprocess.CalledProcessError: If the shuffle or glove binary fails.
    """
    x_max = float(np.load(_x_max_path(x_max_file)))
    glove_dir = os.path.dirname(__file__)

    os.makedirs(result, exist_ok=True)
    result = os.path.join(result, 'embeddings')
    shuf_file = f'{result}_temp.shuf.bin'

    logger.info('Cooccurrence shuffle -> %s', shuf_file)
    with open(cooccur_file, 'rb') as fin, open(shuf_file, 'wb') as fout:
        subprocess.run([os.path.join(glove_dir, 'glove_build', 'shuffle'),
                        '-verbose', '2', '-memory', '40'],
                       stdin=fin, stdout=fout, check=True)

    logger.info('Train GloVe: %d dims, %d iterations, x_max=%g',
                embedding_size, n_iter, x_max)
    subprocess.run([os.path.join(glove_dir, 'glove_build', 'glove'),
                    '-input-file', shuf_file,
                    '-vocab-file', feature_dict,
                    '-save-file', f'{result}_{embedding_size}',
                    '-gradsq-file', f'{result}_gradsq',
                    '-verbose', '2',
                    '-vector-size', str(embedding_size),
                    '-threads', str(cpus),
                    '-iter', str(n_iter),
                    '-eta', str(lr),
                    '-x-max', str(x_max)],
                   check=True)