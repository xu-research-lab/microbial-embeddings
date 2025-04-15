"""
glove (:mod:`membed.glove`)
===================================
.. currentmodule::
Classes
^^^^^^^
.. autosummary::
   :toctree: generated

   GloVeModel

Functions
^^^^^^^^^
.. autosummary::
   :toctree: generated

   read_biom
   read_fasta_ids
   normalize_table
   merge_cooccur
   compute_cooccur
   build_cooccur_matrix
   cooccur_workflow
   train_glove_model

"""
from math import ceil, sqrt
import sys
from logging import getLogger

import yaml
import biom
import numpy as np
from numba import jit
from scipy.sparse import dok_array, coo_array, tril
from scipy.spatial import distance

from joblib import Parallel, delayed
from scipy.io import mmwrite, mmread
from tqdm import tqdm

import torch
import torch.nn as nn
import torch.nn.init as init
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset

import wandb

logger = getLogger(__name__)


class GloVeDataSet(Dataset):
    def __init__(self, cooccur):
        '''cooccur is `scipy.sparse.coo_array`.'''
        if not cooccur.data.any():
            raise ValueError('There is no co-occurring pairs in the cooccurrence table.')
        self._cooccur = cooccur

    def __getitem__(self, index):
        return self._cooccur.row[index], self._cooccur.col[index], self._cooccur.data[index]

    def __len__(self):
        return len(self._cooccur.data)


class GloVeModel(nn.Module):
    """Implement GloVe model with Pytorch.

    Parameters
    ----------
    feature_size : int
        the number of features
    embedding_size : int
        dimension of feature vector.
    x_max : int
        define weighting function when computing the cost for two feature pairs;
        see the GloVe paper for more details.
        Returns the cost associated with the given weight assignments and
        updates the weights by online AdaGrad in place.
    alpha : float
        define weighting function when computing the cost for two feature pairs;
        see the GloVe paper for more details.
        Returns the cost associated with the given weight assignments and
        updates the weights by online AdaGrad in place.
    feature1 : 1-D list-like array of int.
                   The index for microbes.
    feature2 : 1-D list-like array of int.
               The index for microbes.
    cooccur : 1-D list-like array of flaat.
               The co-occurrence for microbes.
    epochs : int
        iteration times while moving toward a minimum of a loss function.
    batch_size : int
        how many samples per batch to load.

    Attributes
    ----------
    embedding_size : int
        dimension of feature vector.
    feature_size : int
        the number of features.
    x_max : int
        define weighting function when computing the cost for two feature pairs;
        see the GloVe paper for more details.
        Returns the cost associated with the given weight assignments and
        updates the weights by online AdaGrad in place.
    alpha : float
        define weighting function when computing the cost for two feature pairs;
        see the GloVe paper for more details.
        Returns the cost associated with the given weight assignments and
        updates the weights by online AdaGrad in place.
    """

    def __init__(self, feature_size, embedding_size, x_max=100, alpha=3/4):
        super().__init__()

        self.embedding_size = embedding_size
        self.feature_size = feature_size

        self.x_max = x_max
        self.alpha = alpha

        self._glove_dataset = None

        self._v = nn.Embedding(feature_size, embedding_size)
        self._w = nn.Embedding(feature_size, embedding_size)
        self._v_bias = nn.Embedding(feature_size, 1)
        self._w_bias = nn.Embedding(feature_size, 1)

        for param in self.parameters():
            logger.debug('PyTorch model parameter: %r %r' % (type(param), param.size()))
            init.uniform_(param, a=-1, b=1)

    def train(self, cooccur, epochs, device, batch_size=512, learning_rate=0.05):
        """Training GloVe model

        Parameters
        ----------
        cooccur : 1-D list-like array of flaat.
                   The co-occurrence for microbes.
        epochs : int
            iteration times while moving toward a minimum of a loss function.
        device : str
            cpu or gpu.
        batch_size : int
            how many samples per batch to load.
        learning_rate : float
            step size at each iteration while moving toward a minimum of a loss function.
        """
        self._glove_dataset = GloVeDataSet(cooccur)

        optimizer = optim.Adam(self.parameters(), lr=learning_rate)
        scheduler = optim.lr_scheduler.ExponentialLR(optimizer, gamma=0.9)
        glove_dataloader = DataLoader(self._glove_dataset, batch_size, shuffle=True)
        for epoch in range(epochs):
            total_loss = 0
            with tqdm(glove_dataloader, unit='batch') as batch_progress:
                batch_progress.set_description(f'Epoch {epoch}')
                for n, batch in enumerate(batch_progress, 1):
                    optimizer.zero_grad()
                    i_s, j_s, cooccur_s = batch
                    i_s, j_s, cooccur_s = i_s.to(device), j_s.to(device), cooccur_s.to(device)
                    loss = self._loss(i_s, j_s, cooccur_s)
                    loss_val = loss.item()
                    total_loss += loss_val
                    loss.backward()
                    optimizer.step()
                    # batch_progress.set_postfix(loss=f'{loss_val:.5f}', ave_loss=f'{total_loss/n:.5f}')
                scheduler.step()
                wandb.log({"ave_loss": total_loss/cooccur.shape[0]})

    def _loss(self, feature_idx1, feature_idx2, cooccur):
        """"
        Run a single iteration of GloVe training using the given co-occurrence data and the previously computed weight vectors
        biases and accompanying gradient histories.

        Parameters
        ----------
        feature_idx1 : list
            the colindex of the self._feature_embeddings
        feature_idx2 : list
            the colindex the self._feature_symmetry_embeddings
        counts: list
            co-occurrence number.

        Returns
        -------
        float
            cost of a single iteration of GloVe training.
        """
        x_max, alpha = self.x_max, self.alpha

        v = self._v(feature_idx1)
        w = self._w(feature_idx2)
        v_bias = self._v_bias(feature_idx1)
        w_bias = self._w_bias(feature_idx2)
        # TODO: weight can be computed only once and moved to dataset loader (if memory is enough)
        weight = torch.pow(cooccur / x_max, alpha)
        weight[weight > 1] = 1

        distance_expr = (torch.sum(v * w, dim=1) + v_bias + w_bias + torch.log(cooccur)) ** 2

        losses = torch.sum(weight * distance_expr)
        return losses

    def fetch_embedding(self):
        return self._v.weight.data.cpu().numpy() + self._w.weight.data.cpu().numpy()


class PhiDataSet(Dataset):
    def __init__(self, cooccur):
        '''cooccur is `scipy.sparse.coo_array`.'''
        if not cooccur.data.any():
            raise ValueError('There is no co-occurring pairs in the cooccurrence table.')
        data = cooccur.toarray()
        n = int(data.shape[1] * (data.shape[1] - 1) / 2)
        self._cooccur = np.zeros(shape=(n, 1))
        self._row = np.zeros(shape=(n, 1))
        self._col = np.zeros(shape=(n, 1))
        m = 0
        for i in range(data.shape[1]):
            for j in range(i - 1):
                self._cooccur[m] = data[j, i]
                self._row[m] = j
                self._col[m] = i
                m += 1

    def __getitem__(self, index):
        return self._cooccur[index], self._cooccur[index], self._cooccur[index]

    def __len__(self):
        return len(self._cooccur.data)


class PhiModel(nn.Module):
    """Implement GloVe model with Pytorch.

    Parameters
    ----------
    feature_size : int
        the number of features
    embedding_size : int
        dimension of feature vector.
    x_max : int
        define weighting function when computing the cost for two feature pairs;
        see the GloVe paper for more details.
        Returns the cost associated with the given weight assignments and
        updates the weights by online AdaGrad in place.
    alpha : float
        define weighting function when computing the cost for two feature pairs;
        see the GloVe paper for more details.
        Returns the cost associated with the given weight assignments and
        updates the weights by online AdaGrad in place.
    feature1 : 1-D list-like array of int.
                   The index for microbes.
    feature2 : 1-D list-like array of int.
               The index for microbes.
    cooccur : 1-D list-like array of flaat.
               The co-occurrence for microbes.
    epochs : int
        iteration times while moving toward a minimum of a loss function.
    batch_size : int
        how many samples per batch to load.

    Attributes
    ----------
    embedding_size : int
        dimension of feature vector.
    feature_size : int
        the number of features.
    x_max : int
        define weighting function when computing the cost for two feature pairs;
        see the GloVe paper for more details.
        Returns the cost associated with the given weight assignments and
        updates the weights by online AdaGrad in place.
    alpha : float
        define weighting function when computing the cost for two feature pairs;
        see the GloVe paper for more details.
        Returns the cost associated with the given weight assignments and
        updates the weights by online AdaGrad in place.
    """

    def __init__(self, feature_size, embedding_size, lambda_1, lambda_2):
        super().__init__()

        self.lambda_1 = lambda_1
        self.lambda_2 = lambda_2

        self.embedding_size = embedding_size
        self.feature_size = feature_size

        self._phi_dataset = None

        self._w = nn.Embedding(feature_size, embedding_size)

        for param in self.parameters():
            logger.debug('PyTorch model parameter: %r %r' % (type(param), param.size()))
            init.uniform_(param, a=-1, b=1)

    def train(self, cooccur, epochs, device, batch_size=512, learning_rate=0.05):
        """Training GloVe model

        Parameters
        ----------
        cooccur : 1-D list-like array of flaat.
                   The co-occurrence for microbes.
        epochs : int
            iteration times while moving toward a minimum of a loss function.
        device : str
            cpu or gpu.
        batch_size : int
            how many samples per batch to load.
        learning_rate : float
            step size at each iteration while moving toward a minimum of a loss function.
        """
        self._phi_dataset = PhiDataSet(cooccur)

        optimizer = optim.Adam(self.parameters(), lr=learning_rate)
        scheduler = optim.lr_scheduler.ExponentialLR(optimizer, gamma=0.9)
        glove_dataloader = DataLoader(self._phi_dataset, batch_size, shuffle=True)
        for epoch in range(epochs):
            total_loss = 0
            with tqdm(glove_dataloader, unit='batch') as batch_progress:
                batch_progress.set_description(f'Epoch {epoch}')
                for n, batch in enumerate(batch_progress, 1):
                    optimizer.zero_grad()
                    i_s, j_s, cooccur_s = batch
                    i_s, j_s, cooccur_s = i_s.to(device), j_s.to(device), cooccur_s.to(device)
                    loss = self._loss(i_s, j_s, cooccur_s)
                    loss_val = loss.item()
                    total_loss += loss_val
                    loss.backward()
                    optimizer.step()
                    # batch_progress.set_postfix(loss=f'{loss_val:.5f}', ave_loss=f'{total_loss/n:.5f}')
                scheduler.step()
                wandb.log({"ave_loss": total_loss/cooccur.shape[0]})

    def _loss(self, feature_idx1, feature_idx2, cooccur):
        """"
        Run a single iteration of GloVe training using the given co-occurrence data and the previously computed weight vectors
        biases and accompanying gradient histories.

        Parameters
        ----------
        feature_idx1 : list
            the colindex of the self._feature_embeddings
        feature_idx2 : list
            the colindex the self._feature_symmetry_embeddings
        counts: list
            co-occurrence number.

        Returns
        -------
        float
            cost of a single iteration of GloVe training.
        """
        w1 = self._w(feature_idx1)
        w2 = self._w(feature_idx2)
        loss = torch.norm(cooccur - w1 * w2, p='fro') + self.lambda_2 / 2 * (torch.norm(w1, p=1) + torch.norm(w2, p=1))
        return loss

    def fetch_embedding(self):
        return self._w.weight.data.cpu().numpy()


def read_biom(biom_file, *, normalize=False, percentile=False, dense=False):
    """Read in a biom table file.

    Parameters
    ----------
    biom_file : str or file object
        file path to the biom table.
    normalize : bool, optional
        True: percentage of one feature in all features in one sample
        False: Do not do any processing
    percentile : bool, optional
        True: Correct for batch effects using percentile-normalization
        False: Do not do any processing
    dense : bool, optional
        True: Convert to dense array
        False: Return sparse array
    Returns
    -------
    sid : list of str
        the sample ids.
    fid : list of str
        the feature ids.
    data : numpy array (2d) of float
        the table. microbe-by-sample table.
    """
    if normalize and percentile:
        raise ValueError('You do NOT need to normalize AND percentile transform the biom table.')

    if hasattr(biom_file, 'read'):
        table = biom.parse_table(biom_file)
    else:
        table = biom.load_table(biom_file)
    logger.info('Loaded %d features, %d samples' % table.shape)
    sid = table.ids(axis='sample')
    fid = table.ids(axis='observation')

    if normalize:
        table.norm(axis='sample', inplace=True)
        data = table.matrix_data
        if not dense:
            data = data.tocsr()
    elif percentile:
        table.rankdata(axis='sample', inplace=True)
        # divide by the max rank in each sample;
        # use `multiply()` to keep array in csr sparse form.
        data = table.matrix_data.multiply(1 / table.max(axis='sample'))
        if not dense:
            # change back to csr sparse array because muliply() changes it to coo array.
            data = data.tocsr()
    else:
        data = table.matrix_data
    return sid, fid, data.toarray() if dense else data


def read_fasta_ids(fasta_file):
    """Read in the fasta file to fetch its sequence ids.

    Parameters
    ----------
    fasta_file : str or file object
        file path to the biom table.

    Returns
    -------
    out :

    """
    with open(fasta_file, 'r') as f:
        for l in f:
            l = l.strip()
            if l.startswith('>'):
                yield l.split(' ', 1)[0][1:]


@jit(nopython=True, cache=True, fastmath=True)
def cooccur_jaccard_dense(u, v):
    union = (u | v)
    intersect = (u & v)
    return intersect.sum() / union.sum()


@jit(nopython=True, cache=True, fastmath=True)
def cooccur_dice_dense(u, v):
    union = (u | v)
    intersect = (u & v)
    return 2 * intersect.sum() / np.sum(union + intersect)


@jit(nopython=True, cache=True, fastmath=True)
def cooccur_faith_dense(u, v):
    '''A variation of Faith metric.'''
    a = (u & v)
    d = (1 - u) & (1 - v)
    n = len(u)
    return np.sum(a + d / n) / n


@jit(nopython=True, cache=True, fastmath=True)
def cooccur_phi_dense(u, v):
    '''A variation of Phi metric.'''
    a = np.sum(u & v)
    d = np.sum((1 - u) & (1 - v))
    b = np.sum(u) - a
    c = np.sum(v) - a
    m = np.sqrt((a + b) * (a + c) * (b + d) * (c + d))
    return (a * d - b * c) / m if m else 0


@jit(nopython=True, cache=True, fastmath=True)
def cooccur_braycurtis_dense(u, v):
    m = u - v
    min_uv = u - np.maximum(m, 0)
    return min_uv.sum() / (u + v).sum()


@jit(nopython=True, cache=True, fastmath=True)
def cooccur_abundance_dense(u, v):
    """compute co-occurrence strength between 2 microbes.

    Parameters
    ----------
    u : 1-D array of float.
        The abundance vector for microbe.
    v : 1-D array of float.
        The abundance vector for microbe.
    Returns
    -------
    float
         co-occurrence strength of a pair of features across samples.
    """
    m = u - v
    min_uv = u - np.maximum(m, 0)
    return (np.sum(min_uv - np.abs(m) * min_uv)) / len(u)


def cooccur_abundance_sparse(u, v):
    """compute co-occurrence strength between 2 microbes.

    Parameters
    ----------
    u : 1-D sparse array of float.
        The abundance vector for microbe.
    v : 1-D sparse array of float.
        The abundance vector for microbe.
    Returns
    -------
    float
         co-occurrence strength of a pair of features across samples.
    """
    u = u.toarray()[0]
    v = v.toarray()[0]
    return cooccur_abundance_dense(u, v) / len(u)

def gen_even_pairs(n, n_packs):
    '''
    Examples
    --------
    >>> n = 3
    >>> x = np.zeros(shape=(n, n), dtype='int16')
    >>> for beg, end, i, t in gen_even_pairs(n, 2):
    ...     for row in range(beg, end):
    ...         for col in range(row):
    ...             x[row, col] = i
    >>> x
    array([[0, 0, 0],
           [1, 0, 0],
           [2, 2, 0]], dtype=int16)
    '''
    total = n * (n-1) // 2
    pack = total // n_packs
    beg = end = 0
    for i in range(1, n_packs+1):
        end = (1 + sqrt(1 - 4 * (beg - beg**2 - 2*pack))) / 2
        end = ceil(end)
        if end > n:
            end = n
        yield beg, end, i, (beg+end-1) * (end-beg) // 2
        beg = end


def build_cooccur_matrix(X, metric, cpus=1):
    """Build a feature co-occurrence list for the given abundance table.

    Parameters
    ----------
    table : numpy.ndarray
        The abundance table for microbes. Samples are in row and features in column.
    formula : str, optional
        "russellrao" : calculate co-occurrence, only considering russellrao presence and absence.
        "abundance" (default) : calculate co-occurrence using abundance. (ASV1, ASV2)
    cpus : int
        The maximum number of concurrently running jobs, such as the number of Python worker processes.
        If -1 all CPUs are used. If 1 is given, no parallel computing code is used at all, which is useful for debugging.

    Returns
    -------
    list
        Each element (representing a co-occurrence pair) is of the form
        (feature_idx1, feature_idx2, co-occurrence) where `feature_idx1` the colindex of the self._feature_embeddings
        and `feature_idx2` the colindex of the self._feature_symmetry_embeddings
    """
    logger.debug(f'Use {cpus} processes.')
    logger.debug(f'X data type: {type(X)}')
    n = X.shape[0]
    # disable tqdm progress bar if it is above DEBUG logging level.
    disable = logger.getEffectiveLevel() > 10

    def _fill_cooccur_mat(beg, end, i, pairs):
        '''s is slice'''
        logger.debug(f'Begin row {beg} and end row {end} have {pairs} pairs')
        # only store non-zero cooccurrence
        #print('X id:', id(X))
        # dok_array is efficient for incremental construction.
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
        chunks = parallel(delayed(_fill_cooccur_mat)(beg, end, i, pairs) for beg, end, i, pairs in gen_even_pairs(n, cpus))

        cooccur = np.sum(chunks)

    return cooccur


def cooccur_workflow(biom_file, normalize=False, percentile=False, dense=True,
                     metric='russellrao', cooccur_file=None, cpus=1):
    """"compute co-occurrence workflow

    biom_file : str or file object
        file path to the biom table.
    normalize : bool, optional
        True : percentage of one feature in all features in one sample
        False : Do not do any processing
    percentile : bool, optional
        True : Correct for batch effects using percentile-normalization
        False : Do not do any processing
    dense : bool, optional
        True : Convert to dense array
        False : Return sparse array
    metric : str, optional
        "russellrao" : calculate co-occurrence, only considering russellrao presence and absence.
        "abundance" (default) : calculate co-occurrence using abundance. (ASV1, ASV2)
    cpus : int
        The maximum number of concurrently running jobs, such as the number of Python worker processes.
        If -1 all CPUs are used. If 1 is given, no parallel computing code is used at all, which is useful for debugging.
    """
    logger.info(f'Set normalize={normalize}, percentile={percentile}, dense={dense} for biom table.')
    if metric == 'russellrao':
        if normalize or percentile:
            raise ValueError('You do NOT need to normalize or percentile transform the biom table!')
        if dense:
            logger.warning('Are you sure to convert sparse array to dense?! Sparse array is much faster for russellrao metric.')

    sid, fid, data = read_biom(biom_file, normalize=normalize, percentile=percentile, dense=dense)
    logger.info('Compute cooccur matrix using %r...' % metric)

    if metric == 'russellrao':
        # convert abundance to presence/absence
        data = data.astype(bool).astype('int')
        # divided by sample size to normalize the cooccurence across data sets
        cooccur = tril(np.dot(data, data.T), k=-1) / len(sid)
        logger.debug('Done computing cooccurence.')
    else:
        if metric == 'abundance':
            metric = cooccur_abundance_dense if dense else cooccur_abundance_sparse
        if metric == 'braycurtis':
            if not dense:
                raise ValueError(f'{metric} is not supported for sparse array.')
            metric = cooccur_braycurtis_dense
        if metric == 'jaccard':
            data = data.astype(bool)
            if not dense:
                raise ValueError(f'{metric} is not supported for sparse array.')
            metric = cooccur_jaccard_dense
        if metric == 'dice':
            data = data.astype(bool)
            if not dense:
                raise ValueError(f'{metric} is not supported for sparse array.')
            metric = cooccur_dice_dense
        if metric == 'faith':
            data = data.astype(bool)
            if not dense:
                raise ValueError(f'{metric} is not supported for sparse array.')
            metric = cooccur_faith_dense
        if metric == 'phi':
            data = data.astype(bool)
            if not dense:
                raise ValueError(f'{metric} is not supported for sparse array.')
            metric = cooccur_phi_dense
        elif callable(metric):
            pass
        else:
            # use scipy distance functions for other metrics.
            # check scipy.spatial.distance._METRICS_NAMES for all available metrics.
            metric = getattr(distance, metric)
        cooccur = build_cooccur_matrix(data, metric=metric, cpus=cpus)
    logger.debug('Cooccur matrix consumed mem:%.2fMB' % (sys.getsizeof(cooccur)/(1024*1024)))
    if cooccur_file is None:
        cooccur_file = f'{biom_file}.cooccur'

    cooccur = coo_array(cooccur)
    mmwrite(cooccur_file, cooccur)

    return cooccur


def classification(embedding_vector):
    # random tree model
    import pandas as pd
    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import GridSearchCV

    data_dir = "/home/dongbiao/word_embedding_microbiome/Classification_prediction/data/"
    # 加载数据集，并转化成相对丰度
    X_train = biom.load_table(data_dir + "X_train.biom").norm(axis='sample', inplace=True)
    x_test = biom.load_table(data_dir + "X_test.biom").norm(axis='sample', inplace=True)
    y_train = list(np.load(data_dir + "y_train.npy"))
    y_test = list(np.load(data_dir + "y_test.npy"))
    validation_healthy_11757 = biom.load_table(data_dir + "validation_healthy_11757.biom").norm(axis='sample', inplace=True)
    validation_ibd_10317 = biom.load_table(data_dir + "validation_ibd_10317.biom").norm(axis='sample', inplace=True)


    fid = pd.DataFrame(
        {'fid': X_train.ids(axis='observation').astype(int), 'id_feature': range(len(X_train.ids(axis='observation')))})

    agp_taxonomy = pd.read_csv(
        "/home/dongbiao/word_embedding_microbiome/programe_test/x-lang/tests/data/AGP-ready-preval.1_taxonomy.txt",
        sep='\t')

    hits_feature = pd.merge(agp_taxonomy, fid, left_on="greengene_id", right_on="fid")
    qual_vecs_half = embedding_vector[list(hits_feature['id'])]

    X_train = pd.DataFrame(np.dot(X_train.matrix_data.toarray().T[:, list(hits_feature['id_feature'])], qual_vecs_half))
    x_test = pd.DataFrame(np.dot(x_test.matrix_data.toarray().T[:, list(hits_feature['id_feature'])], qual_vecs_half))
    validation_healthy = pd.DataFrame(np.dot(validation_healthy_11757.matrix_data.toarray().T[:, list(hits_feature['id_feature'])], qual_vecs_half))
    validation_ibd = pd.DataFrame(np.dot(validation_ibd_10317.matrix_data.toarray().T[:, list(hits_feature['id_feature'])], qual_vecs_half))

    param_test1 = {'n_estimators': range(40, 1000, 10)}
    gsearch1 = GridSearchCV(estimator=GradientBoostingClassifier(learning_rate=0.1, min_samples_split=2,
                                                                 min_samples_leaf=1, max_depth=3, random_state=10,
                                                                 max_features="log2"),
                            param_grid=param_test1, scoring='roc_auc', cv=5, n_jobs=10)
    gsearch1.fit(X_train, y_train)
    train_auc = gsearch1.best_score_

    clf = GradientBoostingClassifier(n_estimators=gsearch1.best_params_['n_estimators'], learning_rate=0.1, min_samples_split=2, min_samples_leaf=1,
                                     max_depth=3, random_state=10, max_features="log2").fit(X_train, y_train)

    y_pred = clf.predict(x_test)

    AUC_test = roc_auc_score(y_test, y_pred)

    valid_pred = clf.predict(validation_ibd)
    valid_accuracy_ibd = np.sum(valid_pred > 0.5) / validation_ibd.shape[1]

    valid_pred = clf.predict(validation_healthy)
    valid_accuracy_healthy = np.sum(valid_pred < 0.5) / validation_healthy.shape[1]

    return train_auc, AUC_test, valid_accuracy_ibd, valid_accuracy_healthy


def train_glove_model(cooccur_file, config_file, model_file=None, embedding_file=None, feature_ref=None,
                      init=True, device=None, cpus=1):
    """Training GloVe model.

    Parameters
    ----------
    init : bool, optional
        True : init vector model
        False : load previous model
    embedding_size : int
        dimension of feature vector.
    feature_ref : str or file object
        file path to the 16S rRNA reference.
    x_max : int
        define weighting function when computing the cost for two feature pairs; see the GloVe paper for more details.
        Returns the cost associated with the given weight assignments and updates the weights by online AdaGrad in place.
    alpha : float
        define weighting function when computing the cost for two feature pairs; see the GloVe paper for more details.
        Returns the cost associated with the given weight assignments and updates the weights by online AdaGrad in place.
    epochs : int
        iteration times while moving toward a minimum of a loss function.
    batch_size : int
        how many samples per batch to load.
    learning_rate : float
        step size at each iteration while moving toward a minimum of a loss function
    model_file : str or file object
        a path to model
    Returns
    -------
    glove model
    """
    with open(config_file) as f:
        config_defaults = yaml.load(f.read(), Loader=yaml.FullLoader)['parameters']

    wandb.init(config=config_defaults)
    config = wandb.config
    embedding_size = config.embedding_size
    x_max = config.x_max
    alpha = config.alpha
    epochs = config.epochs
    batch_size = config.batch_size
    learning_rate = config.learning_rate

    logger.debug(f'Load {cooccur_file}')
    cooccur = mmread(cooccur_file)
    n = cooccur.shape[0]
    # init vector model or load model
    if init:
        logger.info("Initialize model hyper-parameter")
        model = GloVeModel(n, embedding_size, x_max, alpha)
    else:
        logger.info(f"Load previous model {model_file}")
        model = torch.load(model_file)

    # specify device type
    if device is None:
        device = torch.device('cuda:3' if torch.cuda.is_available() else 'cpu')

    model.to(device)
    wandb.watch(model, log="all")
    model.train(cooccur, epochs, device, batch_size, learning_rate)

    if embedding_file is None:
        embedding_file = f'{cooccur_file}_embeddingsize{embedding_size}-xmax{x_max}-alpha{alpha}-epochs{epochs}-batchsize{batch_size}-learningrate{learning_rate}.gz'
    np.savetxt(embedding_file, model.fetch_embedding())
    # save model for evaluation
    if model_file is None:
        model_file = f'{cooccur_file}_embeddingsize{embedding_size}-xmax{x_max}-alpha{alpha}-epochs{epochs}-batchsize{batch_size}-learningrate{learning_rate}.pt'
    torch.save(model.state_dict(), model_file)

    train_auc, test_auc, valid_accuracy_ibd, valid_accuracy_healthy = classification(model.fetch_embedding())
    wandb.log({"train_auc": train_auc})
    wandb.log({"test_auc": test_auc})
    wandb.log({"valid_accuracy_ibd": valid_accuracy_ibd})
    wandb.log({"valid_accuracy_healthy": valid_accuracy_healthy})


def train_phi_model(cooccur_file, config_file, model_file=None, embedding_file=None, feature_ref=None,
                      init=True, device=None, cpus=1):
    """Training GloVe model.

    Parameters
    ----------
    init : bool, optional
        True : init vector model
        False : load previous model
    embedding_size : int
        dimension of feature vector.
    feature_ref : str or file object
        file path to the 16S rRNA reference.
    x_max : int
        define weighting function when computing the cost for two feature pairs; see the GloVe paper for more details.
        Returns the cost associated with the given weight assignments and updates the weights by online AdaGrad in place.
    alpha : float
        define weighting function when computing the cost for two feature pairs; see the GloVe paper for more details.
        Returns the cost associated with the given weight assignments and updates the weights by online AdaGrad in place.
    epochs : int
        iteration times while moving toward a minimum of a loss function.
    batch_size : int
        how many samples per batch to load.
    learning_rate : float
        step size at each iteration while moving toward a minimum of a loss function
    model_file : str or file object
        a path to model
    Returns
    -------
    glove model
    """
    with open(config_file) as f:
        config_defaults = yaml.load(f.read(), Loader=yaml.FullLoader)['parameters']

    wandb.init(config=config_defaults)
    config = wandb.config
    embedding_size = config.embedding_size
    epochs = config.epochs
    batch_size = config.batch_size
    learning_rate = config.learning_rate
    lambda_1 = config.lambda_1
    lambda_2 = config.lambda_2

    logger.debug(f'Load {cooccur_file}')
    cooccur = mmread(cooccur_file)
    n = cooccur.shape[0]
    # init vector model or load model
    if init:
        logger.info("Initialize model hyper-parameter")
        model = GloVeModel(n, embedding_size, lambda_1, lambda_2)
    else:
        logger.info(f"Load previous model {model_file}")
        model = torch.load(model_file)

    # specify device type
    if device is None:
        device = torch.device('cuda:4' if torch.cuda.is_available() else 'cpu')

    model.to(device)
    wandb.watch(model, log="all")

    model.train(cooccur, epochs, device, batch_size, learning_rate)

    if embedding_file is None:
        embedding_file = f'{cooccur_file}_embeddingsize{embedding_size}-lambda1{lambda_1}-lambda2{lambda_2}-epochs{epochs}-batchsize{batch_size}-learningrate{learning_rate}.gz'
    np.savetxt(embedding_file, model.fetch_embedding())
    # save model for evaluation
    if model_file is None:
        model_file = f'{cooccur_file}_embeddingsize{embedding_size}-lambda1{lambda_1}-lambda2{lambda_2}-epochs{epochs}-batchsize{batch_size}-learningrate{learning_rate}.pt'
    torch.save(model.state_dict(), model_file)

    train_auc, test_auc, valid_accuracy_ibd, valid_accuracy_healthy = classification(model.fetch_embedding())
    wandb.log({"train_auc": train_auc})
    wandb.log({"test_auc": test_auc})
    wandb.log({"valid_accuracy_ibd": valid_accuracy_ibd})
    wandb.log({"valid_accuracy_healthy": valid_accuracy_healthy})
