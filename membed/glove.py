
from math import ceil, sqrt
import sys
import time
import yaml
import subprocess
from logging import getLogger

import os
import biom
import numpy as np
import pandas as pd

from scipy.stats import hypergeom
from scipy.sparse import dok_array, coo_array, tril
from scipy.spatial import distance

from joblib import Parallel, delayed
from scipy.io import mmwrite
from tqdm import tqdm

import torch
import math
from collections import Counter
from torch import nn
from torch.utils import data as tud
from torch.multiprocessing import Process

import wandb
import os

logger = getLogger(__name__)


def gpu(i=0):
    return torch.device(f'cuda:{i}')


def try_all_gpus(numb):
    return [gpu(i) for i in [numb]]


def read_biom(biom_file, *, normalize=False):
    """Read BIOM format file and return sample and feature data.

    Parameters
    ----------
    biom_file : str or file-like object
        Path to the BIOM file or an already opened file handle.
    normalize : bool, optional (default: False)
        If True, normalize the data matrix by dividing each sample by its total sum.
        Ensures each sample's values sum to 1. Mutually exclusive with percentile.

    Returns
    -------
    tuple: (sample_ids, feature_ids, data_matrix)
        sample_ids : list
            List of sample identifiers from the BIOM table
        feature_ids : list
            List of feature/observation identifiers (e.g., taxonomy IDs or gene names)
        data_matrix : scipy.sparse.csr_matrix or numpy.ndarray
            Processed data matrix. Format depends on `dense` parameter:
            - Sparse CSR matrix when dense=False (default)
            - Dense array when dense=True
        table： 
    """
 
    try:
        table = biom.load_table(biom_file)
    except Exception as e:
        raise ValueError(f"Failed to load BIOM data: {str(e)}") from e


    logger.info('Loaded %d features, %d samples' % table.shape)
    sid = table.ids(axis='sample')
    fid = table.ids(axis='observation')

    if normalize:
        table.norm(axis='sample', inplace=True)

    data = table.matrix_data

    return sid, fid, data.toarray(), table


def read_fasta_ids(fasta_file):
    """Read sequence IDs from FASTA file.

    Parameters
    ----------
    fasta_file : str
        Path to the input FASTA file. File content should follow standard
        FASTA format specifications where sequences are prefixed with '>'.

    Yields
    ------
    str
        Sequence identifier parsed from FASTA headers. Extracted using the
        following rules:
        1. Takes the header line after '>' character
        2. Splits at the first whitespace (removes trailing descriptions)
        3. Returns the remaining identifier string
    """
    with open(fasta_file, 'r') as f:
        for l in f:
            l = l.strip()
            if l.startswith('>'):
                yield l.split(' ', 1)[0][1:]


def cooccur_jaccard_dense(u, v):
    """Calculate Jaccard similarity between two vectors.
    """
    union = (u | v)
    intersect = (u & v)
    return intersect.sum() / union.sum()


def cooccur_dice_dense(u, v):
    """Calculate Dice similarity between two vectors.
    """
    union = (u | v)
    intersect = (u & v)
    return 2 * intersect.sum() / np.sum(union + intersect)


def cooccur_faith_dense(u, v):
    """Calculate Faith's phylogenetic similarity.
    """
    a = (u & v)
    d = (1 - u) & (1 - v)
    n = len(u)
    return np.sum(a + d / n) / n


def cooccur_phi_dense(u, v):
    """Calculate Phi coefficient between two vectors.
    """
    a = np.sum(u & v)
    d = np.sum((1 - u) & (1 - v))
    b = np.sum(u) - a
    c = np.sum(v) - a
    m = np.sqrt((a + b) * (a + c) * (b + d) * (c + d))
    return (a * d - b * c) / m if m else 0


def cooccur_braycurtis_dense(u, v):
    """Calculate Bray-Curtis dissimilarity.
    """
    m = u - v
    min_uv = u - np.maximum(m, 0)
    return min_uv.sum() / (u + v).sum()


def cooccur_abundance_dense(u, v):
    """Calculate co-occurrence strength based on abundance.
    """
    m = u - v
    min_uv = u - np.maximum(m, 0)
    return (np.sum(min_uv - np.abs(m) * min_uv))  / len(u)
    # return (np.sum(min_uv - np.abs(m) * min_uv)) + 1

# def cooccur_abundance_dense_len(u, v):
#     """Calculate co-occurrence strength based on abundance.
#     """
#     m = u - v
#     min_uv = u - np.maximum(m, 0)
#     return (np.sum(min_uv - np.abs(m) * min_uv))  / len(u)
#     # return (np.sum(min_uv - np.abs(m) * min_uv)) + 1


def cooccur_effect_size(n, N, M, cooccur):
    """Calculate effect size for co-occurrence.
    """
    exp = hypergeom.expect(args=(M, n, N), loc=0)
    effect_sizes = (cooccur - exp) / M
    return effect_sizes


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
        end = (1 + sqrt(1 - 4 * (beg - beg**2 - 2 * pack))) / 2
        end = ceil(end)
        if end > n:
            end = n
        yield beg, end, i, (beg + end - 1) * (end - beg) // 2
        beg = end


def build_cooccur_matrix(X, metric, cpus=1):
    """Build a sparse co-occurrence matrix by computing pairwise similarities between all sample pairs. Supports parallel acceleration.

    Parameters:
        X (np.ndarray or scipy.sparse.spmatrix):
            A 2D array of shape (n_samples, n_features), e.g., species abundance data.
        metric (str or callable):
            Similarity metric to compute between rows. If "cooccur_effect_size", a predefined effect size
            calculation is used; otherwise, a custom function accepting two 1D arrays and returning a float.
        cpus (int, optional):
            Number of CPU cores to use for parallel processing. Defaults to 1 (sequential).

    Returns:
        cooccur (scipy.sparse.coo_array):
            A sparse matrix in COOrdinate format, containing similarity scores. Only the lower-triangular
            part (excluding the diagonal) is filled to avoid redundancy.
    """

    logger.debug(f'Use {cpus} processes.')
    logger.debug(f'X data type: {type(X)}')
    n = X.shape[0]
    m = X.shape[1]
    disable = logger.getEffectiveLevel() > 10

    if metric is cooccur_effect_size:
        data = tril(np.dot(X, X.T))
        data = data.toarray()
        num = np.sum(X, axis=1)

        def _fill_cooccur_mat(beg, end, i, pairs):
            logger.debug(f'Begin row {beg} and end row {end} have {pairs} pairs')
            cooccur_chunk = dok_array((n, n), dtype='float32')
            with tqdm(total=pairs, position=i, desc=f'CPU core {i}', disable=disable) as pbar:
                for row in range(beg, end):
                    for col in range(0, row):
                        pbar.update()
                        v = metric(num[row], num[col], m, data[row, col])
                        cooccur_chunk[row, col] = v
            return cooccur_chunk

        with Parallel(n_jobs=cpus) as parallel:
            chunks = parallel(
                delayed(_fill_cooccur_mat)(beg, end, i, pairs)
                for beg, end, i, pairs in gen_even_pairs(n, cpus))
            cooccur = np.sum(chunks)
    else:
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


def glove_input(data, cooccur_file, x_max_file, percentile_num=10):
    """Prepare GloVe input from co-occurrence data.

    Parameters:
        data (coo_array): Co-occurrence matrix.
        cooccur_file (str):
            Output path for co-occurrence pairs. Stores all non-zero co-occurrence pairs (word1, word2, value)
            with symmetric expansion (both (i,j) and (j,i) are stored). File format: binary packed array of
            dtype [('word1', 'i4'), ('word2', 'i4'), ('value', 'f8')].
        x_max_file (str): Path to save x_max value.
            Output path for x_max value.  Stores the 99th percentile of co-occurrence values.  In GloVe models:
            - The x_max hyperparameter caps the weighting function's upper bound
            - When co-occurrence value x exceeds x_max, weight w(x) is fixed at 1.0
            - This prevents high-frequency pairs from dominating training and mitigates outlier effects.
    """
    dt = np.dtype([('word1', 'int32'), ('word2', 'int32'), ('value', 'float64')])
    row = np.append(data.row, data.col)
    col = np.append(data.col, data.row)
    value = np.append(data.data, data.data)

    temp = np.empty((len(row), ), dt)
    for i in range(len(row)):
        temp[i] = (row[i] + 1, col[i] + 1, value[i])

    temp.tofile(cooccur_file)

    percentile_90 = np.percentile(data.data, percentile_num)
    np.save(x_max_file, percentile_90)


def get_feature_dict(biom_file, feature_dict):
    """Count the number of non-zero occurrences for each feature across all samples in the dataset.

    Parameters:
        biom_file (str): Path to input BIOM file containing feature abundance data.
        feature_dict (str): Output file path to save the feature dictionary. File format:
            Two columns: [feature_id] [non_zero_count], separated by spaces.
    """
    table = biom.load_table(biom_file)
    tmp = pd.DataFrame({
        'id': table.ids(axis='observation'),
        'value': table.nonzero_counts(axis='observation')
    })
    tmp.to_csv(feature_dict, sep=" ", index=False, header=None)


def buildCooccuranceMatrix(biom_file):
    """Build a co-occurrence matrix with position-based weighting for feature pairs.

    Parameters:
        biom_file (str): Path to BIOM file containing feature abundance data.

    Returns:
        np.ndarray: Symmetric co-occurrence matrix of shape (n_features, n_features), 
                    where values are normalized by total sample count.
    """
    table = biom.load_table(biom_file)
    fid = table.ids(axis='observation')
    fid_dict = {feature: ind for ind, feature in enumerate(list(fid))}
    feature_size = len(fid_dict)
    cooccurance_matrix = np.zeros((feature_size, feature_size), dtype=float)
    table = table.matrix_data.toarray()
    num_samples = table.shape[1]

    for i in tqdm(range(num_samples)):
        n = table[:, i]
        m = fid[np.argsort(-n)][0:len(n[n > 0])]
        feature_id = [fid_dict.get(feature) for feature in m]
        maxlength = len(feature_id)
        for j, center_word_id in enumerate(feature_id[:-1]):
            window_indices = np.array(list(range(j + 1, maxlength)))
            for k in window_indices:
                context_word_id = feature_id[k]
                cooccurance_matrix[center_word_id][context_word_id] += 1.0 / (k-j)

    print(">>>>> Save co-occurance matrix completed.")
    return cooccurance_matrix / num_samples


def cooccur_workflow(biom_file,
                     cooccur_file,
                     x_max_file,
                     normalize=False,
                     percentile=False,
                     dense=True,
                     metric='russell_rao',
                     cpus=1,
                     percentile_num=1):
    """Workflow for co-occurrence analysis with multiple similarity metrics.

    Parameters:
        biom_file (str): Input BIOM file path.
        cooccur_file (str): Output path for co-occurrence data (GloVe format).
        x_max_file (str): Output path to save x_max value (used in GloVe training).
        normalize (bool): If True, normalize feature counts by sample totals.
        percentile (bool): If True, apply percentile normalization to data.
        dense (bool): Use dense matrix representation (faster for small datasets).
        metric (str): Similarity metric. Options: 'russell_rao', 'glove', 'jaccard', etc.
        cpus (int): Number of CPU cores for parallel processing.

    Generates table.co and xmax_file.npy for GloVe model training.
    """
    logger.info(
        f'Set normalize={normalize}, percentile={percentile}, dense={dense} for biom table.'
    )
    if metric == 'russellrao':
        if normalize or percentile:
            raise ValueError(
                'You do NOT need to normalize or percentile transform the biom table!'
            )
        if dense:
            logger.warning(
                'Are you sure to convert sparse array to dense?! Sparse array is much faster for russellrao metric.'
            )

    sid, fid, data,table = read_biom(biom_file,
                               normalize=normalize,)
    
    logger.info('Compute cooccur matrix using %r...' % metric)

    if metric == 'russell_rao': 
        data = table.matrix_data.toarray()
        data = data.astype(bool).astype('int')
        cooccur = tril(np.dot(data, data.T), k=-1) / len(sid)
        logger.debug('Done computing cooccurence.')

    elif metric == 'abundance_percentile':

        table.rankdata(axis='sample', inplace=True)
        data = table.matrix_data.multiply(1 / table.max(axis='sample')).toarray()
        metric = cooccur_abundance_dense 

    elif metric == 'braycurtis':
        if not dense:
            raise ValueError(f'{metric} is not supported for sparse array.')
        metric = cooccur_braycurtis_dense

    elif metric == 'jaccard':
        data = data.astype(bool)
        if not dense:
            raise ValueError(f'{metric} is not supported for sparse array.')
        metric = cooccur_jaccard_dense

    elif metric == 'dice':
        data = data.astype(bool)
        if not dense:
            raise ValueError(f'{metric} is not supported for sparse array.')
        metric = cooccur_dice_dense

    elif metric == 'faith':
        data = data.astype(bool)
        if not dense:
            raise ValueError(f'{metric} is not supported for sparse array.')
        metric = cooccur_faith_dense

    elif metric == 'phi':
        data = data.astype(bool)
        if not dense:
            raise ValueError(f'{metric} is not supported for sparse array.')
        metric = cooccur_phi_dense

    elif metric == 'effect_size':
        data = data.astype(bool).astype('int')
        if not dense:
            raise ValueError(f'{metric} is not supported for sparse array.')
        metric = cooccur_effect_size

    elif callable(metric):
        pass

    else:
        metric = getattr(distance, metric)

    if metric == 'glove':
        cooccur = buildCooccuranceMatrix(biom_file)
    else:
        cooccur = build_cooccur_matrix(data, metric=metric, cpus=cpus)

    logger.debug('Cooccur matrix consumed mem:%.2fMB' %
                 (sys.getsizeof(cooccur) / (1024 * 1024)))
    
    if cooccur_file is None:
        cooccur_file = f'{biom_file}.cooccur'
    
    cooccur = coo_array(cooccur)
    glove_input(cooccur, cooccur_file, x_max_file, percentile_num)

def train_glove_model(cooccur_file,
                      x_max_file,
                      feature_dict,
                      result,
                      lr=0.05,
                      cpus=1,
                      iter=10,
                      embedding_size=100):
    """Train GloVe model using original C implementation.

    Parameters:
        cooccur_file (str): 
            Path to co-occurrence data file generated by `cooccur_workflow`.
            Expected format: space-separated triples [feature_i, feature_j, cooccur_value].
        x_max_file (str): 
            Path to file containing the x_max value (float) for GloVe's weighting function.
            Typically generated by `cooccur_workflow`.
        feature_dict (str): 
            Path to feature dictionary mapping feature IDs to metadata.
            Format: space-separated [feature_id] [non_zero_count] (from `get_feature_dict`).
        result (str): 
            Output path prefix for saving:
            - Embeddings: embeddings.txt
            - Temporary file: embeddings_temp.shuf.bin
            - Gradients: embeddings_gradsq
        lr (float, optional): 
            Learning rate (eta) for AdaGrad optimization. Defaults to 0.05.
        cpus (int, optional): 
            Number of CPU cores for parallel processing. Defaults to 1.
        iter (int, optional): 
            Number of training iterations over the dataset. Defaults to 10.
        embedding_size (int, optional): 
            Dimension of output embeddings. Defaults to 100.
    """
    x_max = float(np.load(x_max_file))
    glove_dir = os.path.dirname(__file__)
    result = f"{result}/embeddings"
    subprocess.call("echo Cooccurent shuffle", shell=True)
    subprocess.call(f"{glove_dir}/glove_build/shuffle -verbose 2 -memory 20 \
                    < {cooccur_file} > {result}_temp.shuf.bin",
                    shell=True)
    subprocess.call(f"{glove_dir}/glove_build/glove -input-file {result}_temp.shuf.bin \
                    -vocab-file {feature_dict} -save-file {result} \
                    -gradsq-file {result}_gradsq -verbose 2 \
                    -vector-size {embedding_size} -threads {cpus} \
                    -iter {iter}  -eta {lr} -x-max {x_max}",
                    shell=True)
