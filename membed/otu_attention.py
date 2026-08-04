import time
import csv
import gc
import os
import biom
import random
import numpy as np
import pandas as pd
from torch import nn
import matplotlib.pyplot as plt

import torch
from torch.utils.data import Dataset, DataLoader, Subset
from sklearn.metrics import confusion_matrix, f1_score, matthews_corrcoef, accuracy_score
from sklearn.model_selection import StratifiedShuffleSplit, GroupShuffleSplit
from sklearn import metrics


def gpu(i=0):
    """Select a CUDA device by index.

    Parameters
    ----------
    i : int, optional
        Zero-based CUDA device index. Default is 0.

    Returns
    -------
    torch.device
        Device handle referring to ``cuda:i``.
    """
    return torch.device(f'cuda:{i}')


def num_gpus():
    """Count the visible CUDA devices.

    Sets ``CUDA_DEVICE_ORDER`` to ``PCI_BUS_ID`` and restricts
    ``CUDA_VISIBLE_DEVICES`` to device ``1`` before counting.

    Returns
    -------
    int
        Number of CUDA devices visible to the process.
    """
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    os.environ["CUDA_VISIBLE_DEVICES"] = "1"
    return torch.cuda.device_count()


def try_all_gpus(numb):
    """Wrap a single GPU index in the device list used by the training loops.

    Parameters
    ----------
    numb : int
        Zero-based CUDA device index.

    Returns
    -------
    list of torch.device
        Single-element list holding the selected device.
    """
    return [gpu(i) for i in [numb]]


def read_imdb(biom_table, metadata, labels_col, sample_id_col):
    """Read an OTU table and its sample metadata with flexible column mapping.

    Abundances are converted to within-sample ranks and divided by the per-sample
    maximum rank, so every entry lies in ``(0, 1]``.

    Parameters
    ----------
    biom_table : str
        Path to a BIOM-format OTU table.
    metadata : str
        Path to a tab-separated metadata file.
    labels_col : str
        Metadata column to use as labels.
    sample_id_col : str
        Metadata column holding sample IDs; must match the BIOM sample IDs.

    Returns
    -------
    normalized_table : numpy.ndarray
        Rank-normalized OTU matrix of shape ``(n_samples, n_features)``.
    feature_ids : numpy.ndarray
        OTU feature IDs, aligned with the columns of `normalized_table`.
    labels : list
        Sample labels, aligned with the rows of `normalized_table`.
    sample_ids : numpy.ndarray
        BIOM sample IDs, aligned with the rows of `normalized_table`. Needed to
        label the rows of the saved prediction files.

    Raises
    ------
    KeyError
        If `metadata` does not contain `sample_id_col`, or if it does not cover
        every sample of `biom_table`.
    ValueError
        If `sample_id_col` holds duplicate IDs, which would silently misalign
        the labels, or if `labels_col` does not already hold numeric 0/1
        labels.
    """
    table = biom.load_table(biom_table)
    feature_ids = table.ids(axis="observation")

    try:
        mapping_file = pd.read_csv(
            metadata,
            sep="\t",
            index_col=sample_id_col,  
            dtype={sample_id_col: str},  
            low_memory=False
        )
    except ValueError as e:
        raise KeyError(f"Metadata missing required column: {sample_id_col}") from e

    # A duplicated sample ID makes ``.loc`` return more rows than it was asked
    # for, which shifts the labels out of step with the OTU matrix instead of
    # raising. Catch it here rather than downstream.
    if not mapping_file.index.is_unique:
        duplicated = mapping_file.index[mapping_file.index.duplicated()].unique()
        raise ValueError(
            f"{sample_id_col!r} holds duplicate IDs in {metadata}, which would "
            f"misalign the labels; first offenders: "
            f"{list(duplicated[:5])}")

    ranked_table = table.rankdata(axis='sample', inplace=False)
    max_values = ranked_table.max(axis="sample")
    normalized_table = ranked_table.matrix_data.multiply(1 / max_values)
    normalized_table = normalized_table.toarray().T

    sample_ids = table.ids(axis='sample')
    missing = [sid for sid in sample_ids if sid not in mapping_file.index]
    if missing:
        raise KeyError(
            f"{len(missing)} of {len(sample_ids)} samples in {biom_table} have "
            f"no row in {metadata}; first offenders: {missing[:5]}")

    labels = pd.to_numeric(mapping_file.loc[sample_ids, labels_col],
                           errors='coerce')
    if labels.isna().any():
        offenders = labels.index[labels.isna()][:5].tolist()
        raise ValueError(
            f"column {labels_col!r} must already hold numeric 0/1 labels; "
            f"encode the classes before training. Offending samples: "
            f"{offenders}")
    unexpected = set(labels.unique()) - {0, 1}
    if unexpected:
        raise ValueError(
            f"column {labels_col!r} must hold binary 0/1 labels; found "
            f"{sorted(unexpected)}")
    labels = labels.astype(int).tolist()

    return normalized_table, feature_ids, labels, sample_ids


class Fid:
    """Assign a stable integer index to every OTU in the study.

    The encoder reads a sample as an unordered set of discrete OTUs, so each
    feature ID has to carry an index into the embedding table. This class holds
    that assignment and its inverse, alongside two markers the encoder needs in
    addition to the observed taxa:

    ``'<pad>'``
        Filler for samples whose observed richness falls short of the fixed
        sequence length. Padded positions are masked out of attention and of
        the pooled community representation.
    ``'<unk>'``
        Catch-all for feature IDs that this index does not cover, which is
        where taxa unique to a held-out cohort end up.

    Build the index once from the training table and reuse it for every other
    table in the study. One OTU then keeps one index across cohorts; building a
    separate index per table reassigns indices to different taxa and silently
    scrambles the embedding lookups.

    Parameters
    ----------
    tokens : iterable of str
        OTU feature IDs, e.g. the ``observation`` axis of a BIOM table.
        Duplicates collapse and the survivors are sorted, so the index depends
        only on which IDs are present, not on the order they arrive in.
    reserved_tokens : list of str, optional
        Extra markers to register alongside the taxa, e.g. ``['<pad>']``. They
        are placed *before* the taxa, ahead of ``'<unk>'``, so their indices do
        not move with the naming scheme of the feature IDs: sorting the markers
        in with the taxa would put ``'<pad>'`` after numeric IDs but before
        letter-prefixed ones. Default is None, treated as an empty list.

    Attributes
    ----------
    idx_to_token : list of str
        Feature IDs and markers in index order.
    token_to_idx : dict
        Inverse lookup, from feature ID to index.

    Examples
    --------
    >>> otu_ids = ['OTU_123', 'OTU_456']
    >>> fid = Fid(otu_ids, reserved_tokens=['<pad>'])
    >>> fid.idx_to_token
    ['<pad>', '<unk>', 'OTU_123', 'OTU_456']
    >>> fid['OTU_123']
    2
    >>> fid['OTU_999']  # not in this study, resolves to '<unk>'
    1
    >>> fid.to_tokens(2)
    'OTU_123'

    The markers keep those indices whatever the feature IDs look like:

    >>> Fid(['4479944', '1201'], reserved_tokens=['<pad>']).idx_to_token
    ['<pad>', '<unk>', '1201', '4479944']
    """

    def __init__(self, tokens, reserved_tokens=None):
        """Place the markers first, sort the feature IDs, build the inverse."""
        markers = list(reserved_tokens or []) + ['<unk>']
        seen = set()
        # Keep the markers in the order they were given, minus repeats, so
        # reserved_tokens=['<pad>'] always lands '<pad>' at index 0.
        self.idx_to_token = [t for t in markers
                             if not (t in seen or seen.add(t))]
        self.idx_to_token += sorted(set(tokens) - seen)
        self.token_to_idx = {
            token: idx
            for idx, token in enumerate(self.idx_to_token)
        }

    def __len__(self):
        """Return how many entries the index covers.

        Returns
        -------
        int
            Number of OTUs plus the ``'<pad>'`` and ``'<unk>'`` markers. This is
            the number of rows the embedding table needs.
        """
        return len(self.idx_to_token)

    def __getitem__(self, tokens):
        """Look up the index of one or more OTUs.

        Parameters
        ----------
        tokens : str or list of str or tuple of str
            Feature ID, or sequence of feature IDs, to resolve.

        Returns
        -------
        int or list of int
            Index for a single feature ID, or a list of indices for a sequence.
            IDs this index does not cover resolve to ``'<unk>'``, which is how
            taxa absent from the training table are handled.
        """
        if not isinstance(tokens, (list, tuple)):
            return self.token_to_idx.get(tokens, self.unk)
        return [self.__getitem__(token) for token in tokens]

    def to_tokens(self, indices):
        """Recover the OTU feature IDs behind one or more indices.

        Useful for reading a model back out in taxonomic terms, e.g. naming the
        OTUs that carry the highest attention weight in a sample.

        Parameters
        ----------
        indices : int or sequence of int
            Index, or sequence of indices, to resolve.

        Returns
        -------
        str or list of str
            Feature ID for a scalar index, or a list of feature IDs for a
            sequence, including a sequence of length 1. Markers come back as
            themselves.
        """
        if hasattr(indices, '__len__'):
            return [self.idx_to_token[int(index)] for index in indices]
        return self.idx_to_token[int(indices)]

    @property
    def unk(self):
        """int : Index of ``'<unk>'``, where taxa outside this index resolve."""
        return self.token_to_idx['<unk>']


class load_data_imdb(Dataset):
    """Torch dataset that turns a BIOM OTU table into padded token sequences.

    Each sample is reduced to its most abundant OTUs, encoded as vocabulary
    indices, and padded to a fixed length.

    Parameters
    ----------
    biom_table : str
        Path to a BIOM-format table file (``.biom``).
    metadata : str
        Path to the sample metadata file (TSV format).
    labels_col : str
        Metadata column holding the class labels.
    sample_id_col : str
        Metadata column holding the sample IDs.
    num_steps : int, optional
        Number of OTU positions kept per sample. Default is 500.
    fid : Fid, optional
        Vocabulary to encode this table with. Default is None, which builds a
        fresh vocabulary from this table's own feature IDs. Pass the training
        vocabulary when loading a held-out table, so both are encoded against
        the same index space; OTUs absent from that vocabulary fall back to
        ``'<unk>'``.

    Attributes
    ----------
    fid : Fid
        Vocabulary used to encode the table.
    features : torch.Tensor
        OTU indices of shape ``(n_samples, num_steps)``, dtype int64.
    abundance : torch.Tensor
        Abundance weights of shape ``(n_samples, num_steps)``.
    mask : torch.Tensor
        Attention mask of shape ``(n_samples, num_steps)``, dtype int64, where
        0 marks a padded position or an OTU outside `fid`.
    labels : torch.Tensor
        Sample labels of shape ``(n_samples,)``.
    sample_ids : numpy.ndarray
        BIOM sample IDs, aligned with the rows of `features`.
    otu : numpy.ndarray
        Rank-normalized OTU matrix of shape ``(n_samples, n_features)``.
    """

    def __init__(self, biom_table, metadata, labels_col, sample_id_col,
                 num_steps=500, fid=None):
        """Load the table and build the padded token/abundance/mask tensors."""
        otu, feature_ids, labels, sample_ids = read_imdb(
            biom_table, metadata, labels_col, sample_id_col)
        self.fid = fid if fid is not None else Fid(feature_ids,
                                                   reserved_tokens=['<pad>'])
        self.sample_ids = sample_ids
        self.features, self.abundance, self.mask = self.truncate_pad(
            otu, feature_ids, num_steps)
        self.labels = torch.tensor(labels)
        self.otu = otu

    def truncate_pad(self, otu, feature_ids, num_steps):
        """Truncate or pad each sample's OTUs to a fixed-length sequence.

        Samples with at least `num_steps` non-zero OTUs keep the `num_steps`
        most abundant ones; sparser samples keep all non-zero OTUs and are
        padded with ``'<pad>'`` at abundance 0. Encoding goes through
        ``self.fid``, so a table loaded with a borrowed vocabulary is placed in
        that vocabulary's index space.

        Parameters
        ----------
        otu : numpy.ndarray
            Rank-normalized OTU matrix of shape ``(n_samples, n_features)``.
        feature_ids : numpy.ndarray
            OTU feature IDs aligned with the columns of `otu`.
        num_steps : int
            Target sequence length.

        Returns
        -------
        features : torch.Tensor
            OTU indices of shape ``(n_samples, num_steps)``, dtype int64.
        abundance : torch.Tensor
            Abundance weights of shape ``(n_samples, num_steps)``.
        input_mask : torch.Tensor
            Attention mask of shape ``(n_samples, num_steps)``, dtype int64,
            where 1 marks an OTU this vocabulary covers and 0 marks either a
            padded position or an OTU that resolved to ``'<unk>'``.

        Notes
        -----
        ``'<unk>'`` positions are masked out alongside the padding. They all
        share one embedding row, so leaving them in would let a held-out cohort
        attend to a single meaningless token and, worse, count it in the
        denominator of the mean pooling -- shrinking the representation of
        exactly those samples that carry the most unseen taxa, while the
        training table (whose vocabulary this is) is never affected.
        """
        pad_idx = self.fid['<pad>']
        unk_idx = self.fid.unk
        features = np.full((otu.shape[0], num_steps), pad_idx, dtype=int)
        abundance = np.zeros((otu.shape[0], num_steps))
        for i in range(0, otu.shape[0]):
            nonzero_count = np.count_nonzero(otu[i,])
            if nonzero_count >= num_steps:
                sorted_indices = np.argsort(otu[i,])[::-1]
                features[i, ] = np.array([self.fid[line]
                                          for line in
                                          feature_ids[sorted_indices[:num_steps]]])
                abundance[i, ] = otu[i, sorted_indices[:num_steps]]
            else:
                nonzero_indices = np.nonzero(otu[i,])[0]
                features[i, :nonzero_count] = np.array([self.fid[line]
                                                        for line in
                                                        feature_ids[nonzero_indices]])
                features[i, nonzero_count:] = pad_idx
                abundance[i, :nonzero_count] = otu[i, nonzero_indices]

        input_mask = torch.ones(features.shape[0],
                                features.shape[1],
                                dtype=torch.long)
        # Mask padded positions and OTUs outside this vocabulary alike.
        hidden = (features == pad_idx) | (features == unk_idx)
        input_mask[torch.tensor(hidden)] = 0
        return torch.tensor(features), torch.tensor(abundance), input_mask

    def __call__(self):
        """Return the vocabulary built from the table's feature IDs.

        Returns
        -------
        Fid
            Vocabulary mapping OTU feature IDs to indices.
        """
        return self.fid

    def __getitem__(self, index):
        """Return one sample.

        Parameters
        ----------
        index : int
            Sample index.

        Returns
        -------
        features : torch.Tensor
            OTU indices of shape ``(num_steps,)``.
        abundance : torch.Tensor
            Abundance weights of shape ``(num_steps,)``.
        label : torch.Tensor
            Scalar label for the sample.
        mask : torch.Tensor
            Attention mask of shape ``(num_steps,)``.
        """
        return self.features[index, ], self.abundance[index, ], self.labels[index], self.mask[index]

    def __len__(self):
        """Return the number of samples.

        Returns
        -------
        int
            Number of samples in the dataset.
        """
        return self.otu.shape[0]


def split_train_valid(labels, valid_ratio=0.2, seed=11, groups=None):
    """Split sample indices into a training and a validation part.

    The split is stratified on `labels`, so the class balance of the training
    table is preserved in both parts. When `groups` is given, the split is made
    at the group level instead, keeping every sample of a group on the same
    side; this is what stops repeated samples from one subject appearing in
    both parts.

    Parameters
    ----------
    labels : array_like
        Class label per sample, of shape ``(n_samples,)``.
    valid_ratio : float, optional
        Fraction of samples held out for validation. Default is 0.2.
    seed : int, optional
        Seed for the split, kept separate from the model seed so the partition
        can be held fixed while the model is reseeded. Default is 11.
    groups : array_like, optional
        Group label per sample, e.g. a subject ID. Default is None, which
        splits at the sample level. Stratification is dropped when groups are
        used, since a group carries whole samples rather than single labels.
        Note that `valid_ratio` then applies to the number of *groups*, not the
        number of samples, so the two parts rarely split the samples in exactly
        that ratio.

    Returns
    -------
    train_idx : numpy.ndarray
        Indices of the training part.
    valid_idx : numpy.ndarray
        Indices of the validation part.

    Raises
    ------
    ValueError
        If the split would leave the validation part empty, or if a class is
        too rare to appear in both parts.
    """
    labels = np.asarray(labels)
    n_samples = len(labels)
    n_valid = int(round(n_samples * valid_ratio))
    if n_valid < 1 or n_valid >= n_samples:
        raise ValueError(
            f"valid_ratio={valid_ratio} yields {n_valid} validation samples "
            f"out of {n_samples}; choose a ratio that leaves both parts "
            f"non-empty")

    if groups is not None:
        # GroupShuffleSplit applies test_size to the groups, so the sample-level
        # check above says nothing about whether this split is viable.
        groups = np.asarray(groups)
        n_groups = len(np.unique(groups))
        n_valid_groups = int(round(n_groups * valid_ratio))
        if n_valid_groups < 1 or n_valid_groups >= n_groups:
            raise ValueError(
                f"valid_ratio={valid_ratio} yields {n_valid_groups} validation "
                f"groups out of {n_groups}; with `groups` the ratio applies to "
                f"groups rather than samples, so it must leave both parts "
                f"holding at least one group")
        splitter = GroupShuffleSplit(n_splits=1, test_size=valid_ratio,
                                     random_state=seed)
        train_idx, valid_idx = next(
            splitter.split(np.zeros(n_samples), labels, groups))
    else:
        _, class_counts = np.unique(labels, return_counts=True)
        if class_counts.min() < 2:
            raise ValueError(
                "every class needs at least 2 samples to be stratified across "
                "the split; the rarest class has "
                f"{int(class_counts.min())}")
        splitter = StratifiedShuffleSplit(n_splits=1, test_size=valid_ratio,
                                          random_state=seed)
        train_idx, valid_idx = next(splitter.split(np.zeros(n_samples), labels))

    return train_idx, valid_idx


class otuEmbedding:
    """Pretrained embedding vectors for OTUs, loaded from a text file.

    The file is expected in GloVe format: one OTU identifier per line followed
    by its space-separated vector components. An all-zero ``'<unk>'`` vector is
    inserted at index 0 for OTUs missing from the file.

    Parameters
    ----------
    data_dir : str
        Path to the OTU embedding file.

    Attributes
    ----------
    idx_to_otu : list of str
        Ordered OTU identifier list, starting with ``'<unk>'``.
    idx_to_vec : torch.Tensor
        Embedding matrix of shape ``(n_otus, embedding_dim)``.
    otu_to_idx : dict
        Inverse mapping from OTU identifier to row index in `idx_to_vec`.
    unknown_idx : int
        Row index used for unknown OTUs; always 0.
    """

    def __init__(self, data_dir):
        """Load the embedding file and build the OTU/vector lookups."""
        self.idx_to_otu, self.idx_to_vec = self._load_embedding(data_dir)
        self.unknown_idx = 0
        self.otu_to_idx = {
            otu: idx
            for idx, otu in enumerate(self.idx_to_otu)
        }

    def _load_embedding(self, data_dir):
        """Parse an embedding file into an OTU list and a vector matrix.

        Lines carrying a single component are skipped, which drops the
        dimension header some embedding formats put on the first line.

        Parameters
        ----------
        data_dir : str
            Path to the OTU embedding file.

        Returns
        -------
        idx_to_otu : list of str
            Ordered OTU identifier list, starting with ``'<unk>'``.
        idx_to_vec : torch.Tensor
            Embedding matrix of shape ``(n_otus, embedding_dim)``, with an
            all-zero first row for ``'<unk>'``.
        """
        idx_to_otu, idx_to_vec = ['<unk>'], []
        with open(data_dir, 'r') as f:
            for line in f:
                elems = line.rstrip().split(' ')
                otu, elems = elems[0], [float(elem) for elem in elems[1:]]
                if len(elems) > 1:
                    idx_to_otu.append(otu)
                    idx_to_vec.append(elems)
        idx_to_vec = [[0] * len(idx_to_vec[0])] + idx_to_vec
        return idx_to_otu, torch.tensor(idx_to_vec)

    def __getitem__(self, otus):
        """Look up embedding vectors for a sequence of OTUs.

        Parameters
        ----------
        otus : list of str
            OTU identifiers to look up. Identifiers absent from the file
            receive the all-zero ``'<unk>'`` vector.

        Returns
        -------
        torch.Tensor
            Embedding vectors of shape ``(len(otus), embedding_dim)``.
        """
        indices = [
            self.otu_to_idx.get(otu, self.unknown_idx) for otu in otus
        ]
        vecs = self.idx_to_vec[torch.tensor(indices)]
        return vecs

    def __len__(self):
        """Return the number of OTUs with embedding vectors.

        Returns
        -------
        int
            Number of OTUs, including ``'<unk>'``.
        """
        return len(self.idx_to_otu)


class ScaledDotProductAttention(nn.Module):
    """Scaled dot-product attention.

    Computes ``softmax(Q @ K.T / sqrt(d_k)) @ V``, with masked positions
    driven to a large negative score before the softmax.

    Parameters
    ----------
    d_k : int
        Dimension of the key vectors, used as the scaling factor.
    """

    def __init__(self, d_k):
        """Store the key dimension used for score scaling."""
        super(ScaledDotProductAttention, self).__init__()
        self.d_k = d_k

    def forward(self, q, k, v, attn_mask):
        """Apply scaled dot-product attention.

        Parameters
        ----------
        q : torch.Tensor
            Queries of shape ``(batch, n_heads, seq_len, d_k)``.
        k : torch.Tensor
            Keys of shape ``(batch, n_heads, seq_len, d_k)``.
        v : torch.Tensor
            Values of shape ``(batch, n_heads, seq_len, d_v)``.
        attn_mask : torch.Tensor
            Boolean mask broadcastable to ``(batch, n_heads, seq_len, seq_len)``
            -- normally ``(batch, 1, seq_len, seq_len)`` -- where ``True`` marks
            a position to suppress.

        Returns
        -------
        output : torch.Tensor
            Attended values of shape ``(batch, n_heads, seq_len, d_v)``.
        attn_weights : torch.Tensor
            Attention weights of shape ``(batch, n_heads, seq_len, seq_len)``.
        """
        attn_score = torch.matmul(q, k.transpose(-1, -2)) / np.sqrt(self.d_k)
        attn_score.masked_fill_(attn_mask, -1e9)
        attn_weights = nn.Softmax(dim=-1)(attn_score)
        output = torch.matmul(attn_weights, v)
        return output, attn_weights


class MultiHeadAttention(nn.Module):
    """Multi-head self-attention.

    Projects the inputs into `n_heads` subspaces of size ``d_model // n_heads``,
    attends within each, then concatenates and re-projects to `d_model`.

    Parameters
    ----------
    d_model : int
        Model dimension. Must be divisible by `n_heads`.
    n_heads : int
        Number of attention heads.

    Attributes
    ----------
    d_k, d_v : int
        Per-head key and value dimension, both ``d_model // n_heads``.
    """

    def __init__(self, d_model, n_heads):
        """Build the Q/K/V projections, the attention op, and the output layer."""
        super(MultiHeadAttention, self).__init__()
        self.n_heads = n_heads
        self.d_k = self.d_v = d_model // n_heads
        self.WQ = nn.Linear(d_model, d_model)
        self.WK = nn.Linear(d_model, d_model)
        self.WV = nn.Linear(d_model, d_model)
        self.scaled_dot_product_attn = ScaledDotProductAttention(self.d_k)
        self.linear = nn.Linear(n_heads * self.d_v, d_model)

    def forward(self, Q, K, V, attn_mask):
        """Apply multi-head attention.

        Parameters
        ----------
        Q : torch.Tensor
            Queries of shape ``(batch, seq_len, d_model)``.
        K : torch.Tensor
            Keys of shape ``(batch, seq_len, d_model)``.
        V : torch.Tensor
            Values of shape ``(batch, seq_len, d_model)``.
        attn_mask : torch.Tensor
            Boolean mask of shape ``(batch, seq_len, seq_len)``, broadcast
            across heads internally, where ``True`` marks a position to suppress.

        Returns
        -------
        output : torch.Tensor
            Attention output of shape ``(batch, seq_len, d_model)``.
        attn_weights : torch.Tensor
            Attention weights of shape ``(batch, n_heads, seq_len, seq_len)``.

        Notes
        -----
        The mask is broadcast over the head axis rather than materialized per
        head, which matters at the sequence lengths this encoder runs at: a
        ``(batch, n_heads, seq_len, seq_len)`` copy is the largest tensor in
        the layer after the scores themselves.
        """
        batch_size = Q.size(0)
        q_heads = self.WQ(Q).view(batch_size, -1, self.n_heads,
                                  self.d_k).transpose(1, 2)
        k_heads = self.WK(K).view(batch_size, -1, self.n_heads,
                                  self.d_k).transpose(1, 2)
        v_heads = self.WV(V).view(batch_size, -1, self.n_heads,
                                  self.d_v).transpose(1, 2)
        attn_mask = attn_mask.unsqueeze(1)
        attn, attn_weights = self.scaled_dot_product_attn(
            q_heads, k_heads, v_heads, attn_mask)
        attn = attn.transpose(1, 2).contiguous().view(batch_size, -1,
                                                      self.n_heads * self.d_v)
        output = self.linear(attn)
        return output, attn_weights


class EncoderLayer(nn.Module):
    """Attention-only transformer encoder layer.

    Applies multi-head self-attention with dropout and a residual connection
    followed by layer normalization. There is no position-wise feed-forward
    sub-layer, so each layer mixes information across OTUs but applies no
    per-OTU nonlinearity.

    Parameters
    ----------
    d_model : int
        Model dimension.
    n_heads : int
        Number of attention heads.
    p_drop : float
        Dropout probability.
    """

    def __init__(self, d_model, n_heads, p_drop):
        """Build the attention sub-layer, its dropout, and its norm."""
        super(EncoderLayer, self).__init__()
        self.mha = MultiHeadAttention(d_model, n_heads)
        self.dropout1 = nn.Dropout(p_drop)
        self.layernorm1 = nn.LayerNorm(d_model, eps=1e-6)

    def forward(self, inputs, attn_mask):
        """Run the attention sub-layer with a residual connection.

        Parameters
        ----------
        inputs : torch.Tensor
            Input of shape ``(batch, seq_len, d_model)``.
        attn_mask : torch.Tensor
            Boolean mask of shape ``(batch, seq_len, seq_len)``, where ``True``
            marks a position to suppress.

        Returns
        -------
        attn_outputs : torch.Tensor
            Normalized output of shape ``(batch, seq_len, d_model)``.
        attn_weights : torch.Tensor
            Attention weights of shape ``(batch, n_heads, seq_len, seq_len)``.
        """
        attn_outputs, attn_weights = self.mha(inputs, inputs, inputs,
                                              attn_mask)
        attn_outputs = self.dropout1(attn_outputs)
        attn_outputs = self.layernorm1(inputs + attn_outputs)
        return attn_outputs, attn_weights


class OtuAttentionEncoder(nn.Module):
    """Attention encoder over abundance-weighted OTU sets.

    Each sample is treated as a *set* of OTUs rather than a sequence. OTU
    embeddings are scaled by their relative abundance, refined by a stack of
    attention-only :class:`EncoderLayer`, mean-pooled over the unmasked
    positions into a sample representation, and mapped to a single output by
    a linear head.

    This departs from a standard transformer encoder in two ways: there is no
    positional encoding, so the model is invariant to the order in which the
    OTUs are listed; and the encoder layers carry no position-wise
    feed-forward sub-layer, so each layer mixes information across OTUs
    without applying a per-OTU nonlinearity.

    Parameters
    ----------
    otu_size : int
        Size of the OTU vocabulary.
    d_model : int, optional
        Model dimension. Default is 128.
    n_layers : int, optional
        Number of encoder layers. Default is 6.
    n_heads : int, optional
        Number of attention heads per layer. Default is 8.
    p_drop : float, optional
        Dropout probability. Default is 0.1.
    pad_id : int, optional
        Index of the padding token. Default is 0.

    Attributes
    ----------
    embedding : torch.nn.Embedding
        Token embedding table of shape ``(otu_size, d_model)``.
    layers : torch.nn.ModuleList
        The `n_layers` encoder layers.
    mlp : torch.nn.Sequential
        Linear head mapping the pooled embedding to a scalar.
    """

    def __init__(self,
                 otu_size,
                 d_model=128,
                 n_layers=6,
                 n_heads=8,
                 p_drop=0.1,
                 pad_id=0):
        """Build the embedding table, the encoder stack, and the output head."""
        super(OtuAttentionEncoder, self).__init__()
        self.embedding = nn.Embedding(otu_size, d_model, padding_idx=pad_id)
        self.layers = nn.ModuleList([
            EncoderLayer(d_model, n_heads, p_drop)
            for _ in range(n_layers)
        ])
        self.sigmoid = nn.Sigmoid()
        self.pad_id = pad_id
        self.dropout = nn.Dropout(p=p_drop)
        self.mlp = nn.Sequential(
                    nn.Linear(d_model, 1), 
                    )

    def forward(self, x, weight, mask, encoder=False):
        """Encode a batch of samples and produce a scalar logit per sample.

        Parameters
        ----------
        x : torch.Tensor
            Token indices of shape ``(batch, seq_len)``.
        weight : torch.Tensor
            Abundance weights of shape ``(batch, seq_len)``, multiplied into the
            token embeddings.
        mask : torch.Tensor
            Attention mask of shape ``(batch, seq_len)``, where 1 marks a
            position to attend over and to include in the mean pooling, and 0
            marks padding or an OTU outside the vocabulary.
        encoder : bool, optional
            If True, return intermediate representations instead of the
            prediction. Default is False.

        Returns
        -------
        outputs : torch.Tensor
            Raw logits of shape ``(batch, 1)``; the caller applies the sigmoid.
            Returned when ``encoder`` is False.
        attn_weights : torch.Tensor
            Attention weights from the final layer, of shape
            ``(batch, n_heads, seq_len, seq_len)``. Returned when ``encoder`` is False.
        inputs_1 : torch.Tensor
            Abundance-weighted input embeddings of shape
            ``(batch, seq_len, d_model)``. Returned when ``encoder`` is True.
        embedding : torch.Tensor
            Mean-pooled sample representation of shape ``(batch, d_model)``.
            Returned when ``encoder`` is True.

        Notes
        -----
        `mask` is the single source of truth: the same tensor drives both the
        attention mask and the pooling, so a position hidden from one is hidden
        from the other. Deriving the attention mask from `pad_id` instead would
        let ``'<unk>'`` positions -- padding's counterpart for taxa the
        vocabulary does not cover -- slip into both.
        """
        inputs = self.embedding(x)
        inputs = inputs.permute(2, 0, 1) * weight
        inputs = inputs.permute(1, 2, 0)
        inputs_1 = inputs.clone()
        attn_pad_mask = self.get_attention_padding_mask(mask)

        for layer in self.layers:
            inputs, attn_weights = layer(inputs.float(), attn_pad_mask)

        embedding_sum = inputs * mask.unsqueeze(-1)
        # A sample can in principle carry no covered OTU at all; clamping keeps
        # the pooling finite instead of handing NaNs to the loss.
        n_valid = mask.sum(1, keepdim=True).clamp(min=1)
        embedding = embedding_sum.sum(dim=1, keepdim=False) / n_valid
        outputs = self.dropout(embedding)
        outputs = self.mlp(outputs.float())

        if encoder == False:
            return outputs, attn_weights
        else:
            return inputs_1, embedding

    def get_attention_padding_mask(self, mask):
        """Build the boolean mask that hides masked-out keys from every query.

        Parameters
        ----------
        mask : torch.Tensor
            Attention mask of shape ``(batch, seq_len)``, where 1 marks a
            position to keep.

        Returns
        -------
        torch.Tensor
            Boolean mask of shape ``(batch, seq_len, seq_len)``, ``True`` where
            the key position is to be suppressed. The query axis is an expanded
            view, not a copy.
        """
        return mask.eq(0).unsqueeze(1).expand(-1, mask.size(1), -1)



def to_numpy(t):
    """Detach a tensor and move it to a NumPy array on the host.

    Parameters
    ----------
    t : torch.Tensor
        Tensor on any device, with or without a gradient.

    Returns
    -------
    numpy.ndarray
        Array holding the same values.
    """
    return t.detach().cpu().numpy()


def fit_threshold(Y, prob):
    """Choose the decision threshold that maximizes Youden's J statistic on the ROC curve.

    Call this on the **validation** set only. Fitting a threshold on the same
    data the metrics are reported for makes those metrics optimistic, which is
    why :func:`evaluate_cls` takes the threshold as a required argument rather
    than deriving it.

    Parameters
    ----------
    Y : array_like
        True binary labels of shape ``(n_samples,)``.
    prob : array_like
        Predicted probabilities of the positive class, shape ``(n_samples,)``.

    Returns
    -------
    float
        Threshold maximizing Youden's J statistic (TPR - FPR) along the ROC curve.
        Apply it as ``prob >= threshold``, which is the convention
        :func:`sklearn.metrics.roc_curve` scored the candidates
        under; testing ``prob > threshold`` instead flips the very sample the
        threshold was read off and scores below the maximum this reports.
    """
    Y = np.asarray(Y).astype('int')
    prob = np.asarray(prob)
    fpr, tpr, thresholds = metrics.roc_curve(Y, prob, pos_label=1)
    j_scores = tpr - fpr
    best = int(np.argmax(j_scores))
    return float(thresholds[best])


def evaluate_cls(Y, prob, threshold):
    """Compute binary classification metrics at a fixed decision threshold.

    `threshold` has no default on purpose: it must be supplied by the caller,
    normally the value :func:`fit_threshold` produced on the validation set.
    That keeps threshold selection out of the set being scored.

    Parameters
    ----------
    Y : array_like
        True binary labels of shape ``(n_samples,)``.
    prob : array_like
        Predicted probabilities of the positive class, shape ``(n_samples,)``.
    threshold : float
        Decision threshold; samples with ``prob >= threshold`` are called
        positive, matching the convention :func:`fit_threshold` selects under.

    Returns
    -------
    auc : float
        Area under the ROC curve. Threshold-independent.
    aupr : float
        Area under the precision-recall curve. Threshold-independent.
    cm : numpy.ndarray
        Confusion matrix of shape ``(2, 2)`` at `threshold`.
    f1 : float
        Macro-averaged F1 score at `threshold`.
    mcc : float
        Matthews correlation coefficient at `threshold`.
    acc : float
        Accuracy at `threshold`.
    """
    Y = np.asarray(Y).astype('int')
    prob = np.asarray(prob)
    fpr, tpr, _ = metrics.roc_curve(Y, prob, pos_label=1)
    auc = metrics.auc(fpr, tpr)
    precision, recall, _ = metrics.precision_recall_curve(Y, prob, pos_label=1)
    aupr = metrics.auc(recall, precision)

    y_hat = (prob >= threshold).astype(int)
    cm = confusion_matrix(Y, y_hat)
    f1 = f1_score(Y, y_hat, average='macro')
    acc = accuracy_score(Y, y_hat)
    mcc = matthews_corrcoef(Y, y_hat)

    return auc, aupr, cm, f1, mcc, acc


def save_predictions(path, sample_ids, labels, prob, threshold):
    """Write per-sample predicted probabilities to a CSV file.

    Parameters
    ----------
    path : str
        Destination file path.
    sample_ids : array_like
        Sample identifiers of shape ``(n_samples,)``.
    labels : array_like
        True labels of shape ``(n_samples,)``.
    prob : array_like
        Predicted probabilities of shape ``(n_samples,)``.
    threshold : float
        Threshold used to derive the hard label written alongside `prob`, as
        ``prob >= threshold``.

    Returns
    -------
    None
    """
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['sample_id', 'true_label', 'prob', 'pred_label'])
        for sid, y, p in zip(sample_ids, labels, prob):
            writer.writerow([sid, int(y), float(p), int(float(p) >= threshold)])


class Timer:
    """Stopwatch that accumulates a series of elapsed times.

    Attributes
    ----------
    times : list of float
        Elapsed seconds recorded by each :meth:`stop` call.
    """

    def __init__(self):
        """Initialize the record list and start the first interval."""
        self.times = []
        self.start()

    def start(self):
        """Start a new interval, discarding any interval already running."""
        self.tik = time.time()

    def stop(self):
        """End the current interval and record its duration.

        Returns
        -------
        float
            Seconds elapsed since the last :meth:`start` call.
        """
        self.times.append(time.time() - self.tik)
        return self.times[-1]

    def avg(self):
        """Return the mean recorded interval.

        Returns
        -------
        float
            Average of `times`, in seconds.
        """
        return sum(self.times) / len(self.times)

    def sum(self):
        """Return the total recorded time.

        Returns
        -------
        float
            Sum of `times`, in seconds.
        """
        return sum(self.times)

    def cumsum(self):
        """Return the running total across recorded intervals.

        Returns
        -------
        list of float
            Cumulative sum of `times`, in seconds.
        """
        return np.array(self.times).cumsum().tolist()


def train_batch_ch13(net, X, y, abundance, loss, mask, trainer, devices):
    """Run one training step on a minibatch.

    The model output is squeezed to shape ``(batch,)`` and squashed with a
    sigmoid, so `loss` must be a criterion that consumes probabilities.

    Parameters
    ----------
    net : torch.nn.Module
        Model whose forward pass returns ``(logits, attention_weights)``.
    X : torch.Tensor or list of torch.Tensor
        Token indices of shape ``(batch, seq_len)``, or a list of such tensors.
    y : torch.Tensor
        Binary target labels of shape ``(batch,)``.
    abundance : torch.Tensor
        Relative abundance weights of shape ``(batch, seq_len)``.
    loss : torch.nn.Module
        Loss criterion operating on probabilities, e.g.
        :class:`torch.nn.BCELoss` or :class:`FocalLoss`.
    mask : torch.Tensor
        Attention mask of shape ``(batch, seq_len)``, where 1 marks a valid
        position and 0 marks padding.
    trainer : torch.optim.Optimizer
        Optimizer used to apply the gradient step.
    devices : list of torch.device
        Devices to use; only ``devices[0]`` is used.

    Returns
    -------
    train_loss_sum : torch.Tensor
        Scalar loss for the minibatch.
    pred : torch.Tensor
        Predicted probabilities of shape ``(batch,)``.
    """
    if isinstance(X, list):
        X = [x.to(devices[0]) for x in X]
    else:
        X = X.to(devices[0])
    y = y.to(devices[0], dtype=torch.int64)
    abundance = abundance.to(devices[0])
    mask = mask.to(devices[0])
    net.train()
    trainer.zero_grad()
    pred, _ = net(X, abundance, mask)
    pred = torch.squeeze(pred, dim=1)
    pred = torch.sigmoid(pred)
    l = loss(pred.float(), y.float())
    l.backward()
    trainer.step()
    train_loss_sum = l.sum()
    return train_loss_sum, pred


def set_axes(axes, xlabel, ylabel, xlim, ylim, xscale, yscale, legend):
    """Apply labels, scales, limits, legend, and grid to a matplotlib axes.

    Parameters
    ----------
    axes : matplotlib.axes.Axes
        Axes to configure, modified in place.
    xlabel, ylabel : str
        Axis labels.
    xlim, ylim : tuple of float
        Axis limits, as ``(min, max)``.
    xscale, yscale : str
        Axis scales, e.g. ``'linear'`` or ``'log'``.
    legend : list of str
        Legend entries. Falsy values leave the legend off.

    Returns
    -------
    None
    """
    axes.set_xlabel(xlabel), axes.set_ylabel(ylabel)
    axes.set_xscale(xscale), axes.set_yscale(yscale)
    axes.set_xlim(xlim), axes.set_ylim(ylim)
    if legend:
        axes.legend(legend)
    axes.grid()


class Animator:
    """Incrementally built line plot that is re-saved to disk on every update.

    Parameters
    ----------
    xlabel, ylabel : str, optional
        Axis labels. Default is None.
    legend : list of str, optional
        Legend entries. Default is None, treated as an empty list.
    xlim, ylim : tuple of float, optional
        Axis limits, as ``(min, max)``. Default is None.
    xscale, yscale : str, optional
        Axis scales. Default is ``'linear'``.
    fmts : tuple of str, optional
        Matplotlib format strings, one per series. Default is
        ``('-', 'm--', 'g-.', 'r:')``.
    nrows, ncols : int, optional
        Subplot grid shape. Default is 1. Only the first axes is drawn on.
    figsize : tuple of float, optional
        Figure size in inches. Default is ``(5, 5)``.

    Attributes
    ----------
    X, Y : list of list of float
        Accumulated coordinates, one inner list per series.
    """

    def __init__(self,
                 xlabel=None,
                 ylabel=None,
                 legend=None,
                 xlim=None,
                 ylim=None,
                 xscale='linear',
                 yscale='linear',
                 fmts=('-', 'm--', 'g-.', 'r:'),
                 nrows=1,
                 ncols=1,
                 figsize=(5, 5)):
        """Create the figure and capture the axes configuration."""
        if legend is None:
            legend = []
        self.fig, self.axes = plt.subplots(nrows, ncols, figsize=figsize)
        if nrows * ncols == 1:
            self.axes = [
                self.axes,
            ]
        self.config_axes = lambda: set_axes(self.axes[0], xlabel, ylabel, xlim,
                                            ylim, xscale, yscale, legend)
        self.X, self.Y, self.fmts = None, None, fmts

    def add(self, x, y, plotfile):
        """Append one point per series, redraw the figure, and save it.

        Entries that are None are skipped, which lets a caller update a single
        series while leaving the others untouched.

        Parameters
        ----------
        x : float or sequence of float
            X coordinate per series. A scalar is broadcast to all series.
        y : float or sequence of float
            Y coordinate per series. A scalar is wrapped into a one-element
            sequence.
        plotfile : str
            Path the figure is written to, overwritten on every call.

        Returns
        -------
        None
        """
        if not hasattr(y, "__len__"):
            y = [y]
        n = len(y)
        if not hasattr(x, "__len__"):
            x = [x] * n
        if not self.X:
            self.X = [[] for _ in range(n)]
        if not self.Y:
            self.Y = [[] for _ in range(n)]
        for i, (a, b) in enumerate(zip(x, y)):
            if a is not None and b is not None:
                self.X[i].append(a)
                self.Y[i].append(b)
        self.axes[0].cla()
        for x, y, fmt in zip(self.X, self.Y, self.fmts):
            self.axes[0].plot(x, y, fmt)
        self.config_axes()
        self.fig.savefig(plotfile)
        
        
def predict_iter(net, data_iter, devices, loss=None):
    """Run the model over a loader in evaluation mode, without updating it.

    Parameters
    ----------
    net : torch.nn.Module
        Model to run. Left in evaluation mode on return.
    data_iter : torch.utils.data.DataLoader
        Loader yielding ``(features, abundance, labels, mask)`` batches.
    devices : list of torch.device
        Devices to use; only ``devices[0]`` is used.
    loss : torch.nn.Module, optional
        Criterion to accumulate over the loader. Default is None, which skips
        the loss and returns None in its place.

    Returns
    -------
    mean_loss : float or None
        Mean loss per batch, or None when `loss` was not given.
    prob : numpy.ndarray
        Predicted probabilities of shape ``(n_samples,)``.
    label : numpy.ndarray
        True labels of shape ``(n_samples,)``.
    """
    net.eval()
    total_loss = 0.0
    n_batches = 0
    probs, labels_all = [], []
    with torch.no_grad():
        for features, abundance, labels, mask in data_iter:
            X = features.to(devices[0])
            abundance = abundance.to(devices[0])
            mask = mask.to(devices[0])
            pred, _ = net(X, abundance, mask)
            pred = torch.squeeze(pred, dim=1)
            pred = torch.sigmoid(pred)
            if loss is not None:
                y = labels.to(devices[0], dtype=torch.int64)
                total_loss += float(loss(pred.float(), y.float()))
            probs.append(to_numpy(pred))
            labels_all.append(to_numpy(labels))
            n_batches += 1

    mean_loss = total_loss / n_batches if loss is not None else None
    return mean_loss, np.concatenate(probs), np.concatenate(labels_all)


def evaluate_on_test(net, ckpt_path, test_iter, threshold, devices):
    """Score the held-out test set once, using the selected model and threshold.

    This is the only function that touches the test set. It is called after
    training has finished, so nothing about the test set can influence which
    epoch or which threshold was chosen.

    Parameters
    ----------
    net : torch.nn.Module
        Model instance to load the checkpoint into.
    ckpt_path : str
        Path to the ``state_dict`` saved at the best validation epoch.
    test_iter : torch.utils.data.DataLoader
        Loader over the held-out test set.
    threshold : float
        Decision threshold fitted on the validation set, applied unchanged.
    devices : list of torch.device
        Devices to use; only ``devices[0]`` is used.

    Returns
    -------
    metrics_ : dict
        Keys ``auc``, ``aupr``, ``cm``, ``f1``, ``mcc``, ``acc``.
    prob : numpy.ndarray
        Predicted probabilities of shape ``(n_samples,)``.
    label : numpy.ndarray
        True labels of shape ``(n_samples,)``.
    """
    net = net.to(devices[0])
    net.load_state_dict(torch.load(ckpt_path, map_location=devices[0]))
    _, prob, label = predict_iter(net, test_iter, devices)
    auc, aupr, cm, f1, mcc, acc = evaluate_cls(label, prob, threshold)
    metrics_ = {'auc': auc, 'aupr': aupr, 'cm': cm,
                'f1': f1, 'mcc': mcc, 'acc': acc}
    return metrics_, prob, label


def train_cls(net, train_iter, valid_iter, loss, trainer, scheduler, num_epochs,
              devices, ckpt_path, plotfile_loss, plotfile_auc):
    """Train a binary classifier, selecting the best epoch on the validation set.

    The test set is deliberately absent from this function's signature. Every
    decision made here -- which epoch to keep and which decision threshold to
    use -- is made on `valid_iter` alone, so the held-out test set stays
    untouched until :func:`evaluate_on_test` runs afterwards.

    After each epoch the validation AUC is computed; when it improves on the
    best seen so far, the ``state_dict`` is written to `ckpt_path` and the
    optimal threshold **of that same epoch** is recorded, so the returned
    threshold always matches the returned weights.

    Parameters
    ----------
    net : torch.nn.Module
        Model to train. Its logits pass through a sigmoid before `loss`.
    train_iter : torch.utils.data.DataLoader
        Loader yielding ``(features, abundance, labels, mask)`` batches, with
        `features`, `abundance`, and `mask` of shape ``(batch, seq_len)`` and
        `labels` of shape ``(batch,)`` holding 0/1 class labels.
    valid_iter : torch.utils.data.DataLoader
        Validation loader with the same structure as `train_iter`.
    loss : torch.nn.Module
        Loss criterion operating on probabilities, e.g.
        :class:`torch.nn.BCELoss` or :class:`FocalLoss`.
    trainer : torch.optim.Optimizer
        Optimizer used for parameter updates.
    scheduler : torch.optim.lr_scheduler.LRScheduler
        Learning rate scheduler. Accepted but never stepped, so the learning
        rate stays constant.
    num_epochs : int
        Number of training epochs.
    devices : list of torch.device
        Devices to use; only ``devices[0]`` is used.
    ckpt_path : str
        Path the best epoch's ``state_dict`` is saved to.
    plotfile_loss : str
        Path for the train/validation loss curve plot.
    plotfile_auc : str
        Path for the train/validation AUC curve plot.

    Returns
    -------
    dict
        Record of the selected epoch, with keys ``epoch``, ``threshold``,
        ``valid_auc``, ``valid_prob``, ``valid_label``, ``train_loss``,
        ``valid_loss``, and ``train_auc``.

    Raises
    ------
    ValueError
        Propagated from scikit-learn if a split contains only one class, which
        makes the AUC undefined.
    """
    animator_1 = Animator(xlabel='epoch', xlim=[1, num_epochs],
                          legend=['train loss', 'valid loss'])
    animator_2 = Animator(xlabel='epoch', xlim=[1, num_epochs], ylim=[0, 1],
                          legend=['train auc', 'valid auc'])
    net = net.to(devices[0])
    best = None

    for epoch in range(num_epochs):
        train_loss = 0.0
        n_batches = 0
        train_probs, train_labels = [], []
        for features, abundance, labels, mask in train_iter:
            l, pred = train_batch_ch13(
                net, features, labels, abundance, loss, mask, trainer, devices)
            train_loss += float(l)
            train_probs.append(to_numpy(pred))
            train_labels.append(to_numpy(labels))
            n_batches += 1
        train_loss /= n_batches
        train_prob = np.concatenate(train_probs)
        train_label = np.concatenate(train_labels)

        valid_loss, valid_prob, valid_label = predict_iter(
            net, valid_iter, devices, loss)

        train_auc = metrics.roc_auc_score(train_label, train_prob)
        valid_auc = metrics.roc_auc_score(valid_label, valid_prob)
        animator_1.add(epoch + 1, (train_loss, valid_loss), plotfile_loss)
        animator_2.add(epoch + 1, (train_auc, valid_auc), plotfile_auc)

        if best is None or valid_auc > best['valid_auc']:
            best = {'epoch': epoch + 1,
                    'threshold': fit_threshold(valid_label, valid_prob),
                    'valid_auc': valid_auc,
                    'valid_prob': valid_prob,
                    'valid_label': valid_label,
                    'train_loss': train_loss,
                    'valid_loss': valid_loss,
                    'train_auc': train_auc}
            torch.save(net.state_dict(), ckpt_path)

    v_auc, v_aupr, v_cm, v_f1, v_mcc, v_acc = evaluate_cls(
        best['valid_label'], best['valid_prob'], best['threshold'])
    print(f"[valid] best epoch {best['epoch']}/{num_epochs} | "
          f"train loss {best['train_loss']:.4f}, valid loss {best['valid_loss']:.4f} | "
          f"train auc {best['train_auc']:.3f}, valid auc {v_auc:.3f} | "
          f"valid aupr {v_aupr:.3f}, f1 {v_f1:.3f}, mcc {v_mcc:.3f}, acc {v_acc:.3f} | "
          f"threshold {best['threshold']:.4f}")
    print(f"[valid] confusion_matrix\n{v_cm}")

    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return best


def init_transformer_weights(model, method='xavier', exclude='embedding'):
    """Initialize a model's weights in place, leaving selected ones untouched.

    Linear weights get the chosen initializer and their biases are zeroed;
    layer normalizations are reset to their standard affine state; every other
    module is left alone. Modules whose qualified name contains `exclude` are
    skipped along with everything below them.

    Call this **once, on the whole model**. Do not hand it to
    :meth:`torch.nn.Module.apply`: ``apply`` already walks the module tree, and
    this function walks it again from each node it is handed, so every
    parameter would be re-initialized once per ancestor. ``apply`` also strips
    the name prefixes, which silently defeats `exclude` -- the embedding table
    arrives under the bare name ``'weight'``.

    Parameters
    ----------
    model : torch.nn.Module
        Model to initialize, modified in place.
    method : {'xavier', 'kaiming'}, optional
        Initializer for linear weights. Any other value falls back to a
        standard normal. Default is ``'xavier'``.
    exclude : str, optional
        Substring marking modules to leave untouched. Default is
        ``'embedding'``.

    Returns
    -------
    None

    Notes
    -----
    Layer normalizations are deliberately kept out of the fan-based
    initializers. Their weight is one-dimensional, so a fan-based initializer
    reads its shape as a ``(1, d_model)`` matrix and replaces the unit scale
    with random values of standard deviation ``sqrt(2 / (1 + d_model))`` --
    around 0.14 at ``d_model=100`` -- which is not an initialization of a
    normalization layer so much as the removal of one.
    """
    for name, module in model.named_modules():
        if exclude and exclude in name:
            continue
        if isinstance(module, nn.LayerNorm):
            if module.elementwise_affine:
                nn.init.ones_(module.weight)
                nn.init.zeros_(module.bias)
        elif isinstance(module, nn.Linear):
            if method == 'xavier':
                nn.init.xavier_normal_(module.weight)
            elif method == 'kaiming':
                nn.init.kaiming_normal_(module.weight)
            else:
                nn.init.normal_(module.weight)
            if module.bias is not None:
                nn.init.constant_(module.bias, 0)


class FocalLoss(nn.Module):
    """Focal loss for imbalanced binary classification.

    Down-weights well-classified examples by a factor of ``(1 - p) ** gamma``
    so training focuses on the hard ones, with `alpha` balancing the positive
    against the negative class.

    Parameters
    ----------
    alpha : float, optional
        Weight of the positive class; the negative class gets ``1 - alpha``.
        Default is 0.6.
    gamma : float, optional
        Focusing exponent. Larger values suppress easy examples more
        aggressively. Default is 2.
    reduction : {'mean', 'sum'}, optional
        How the per-sample losses are aggregated. Any other value returns the
        unreduced tensor. Default is ``'mean'``.
    devices : list of torch.device, optional
        Devices used to move `alpha` and `gamma` alongside CUDA inputs. Required
        when the inputs are on GPU. Default is None.
    eps : float, optional
        Probabilities are clamped into ``[eps, 1 - eps]`` before the log, so a
        saturated sigmoid cannot hand back a non-finite loss. Default is 1e-7.
    """

    def __init__(self, alpha=0.6, gamma=2, reduction='mean', devices=None,
                 eps=1e-7):
        """Store the focal loss hyperparameters as tensors."""
        super(FocalLoss, self).__init__()
        self.alpha = torch.tensor(alpha)
        self.gamma = torch.tensor(gamma)
        self.reduction = reduction
        self.devices = devices
        self.eps = eps

    def forward(self, inputs, targets):
        """Compute the focal loss.

        Parameters
        ----------
        inputs : torch.Tensor
            Predicted probabilities of the positive class, shape ``(batch,)``,
            already squashed into ``(0, 1)``.
        targets : torch.Tensor
            Binary targets of shape ``(batch,)``.

        Returns
        -------
        torch.Tensor
            Scalar loss under ``'mean'`` or ``'sum'`` reduction, otherwise the
            per-sample losses of shape ``(batch,)``.

        Notes
        -----
        A sigmoid saturates to exactly 0 or 1 in float32 well before the model
        is confident enough for it to be harmless, and ``log(0)`` there would
        put a non-finite value into the gradient and from there into every
        parameter. The probabilities are clamped into ``[eps, 1 - eps]`` first.
        This is the guard :class:`torch.nn.BCELoss` provides on its own by
        clamping its log outputs at -100.
        """
        if inputs.is_cuda and not self.alpha.is_cuda:
            self.alpha = self.alpha.to(self.devices[0])
            self.gamma = self.gamma.to(self.devices[0])
        pt = inputs.clamp(min=self.eps, max=1.0 - self.eps)
        alpha = self.alpha
        F_loss = - alpha * (1 - pt) ** self.gamma * targets * torch.log(pt) - \
               (1 - alpha) * pt ** self.gamma * (1 - targets) * torch.log(1 - pt)
        if self.reduction == 'mean':
            F_loss = torch.mean(F_loss)
        elif self.reduction == 'sum':
            F_loss = torch.sum(F_loss)
        return F_loss


def set_seed(seed=11):
    """Seed the torch, NumPy, and Python random number generators.

    Parameters
    ----------
    seed : int, optional
        Seed applied to all three generators. Default is 11.

    Returns
    -------
    None
    """
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    

def Attention_biom(
        metadata,
        train_biom,
        test_biom,
        embedding_birnn,
        plotfile_loss,
        plotfile_auc,
        labels_col="group",
        sample_id_col="sample_id",
        num_steps=400,
        p_drop=0,
        batch_size=128,
        d_model=100,
        n_layers=2,
        n_heads=2,
        numb=1,
        lr=0.0005,
        weight_decay=0,
        num_epochs=100,
        loss="BCE_loss",
        alpha=0.6,
        glove_embedding=None,
        valid_ratio=0.2,
        split_seed=11,
        group_col=None,
        pred_out=None,):
    """Run the end-to-end training pipeline for OTU-based microbiome analysis.

    The training table is split into a training and a validation part. Training
    and every choice made during it -- the epoch to keep, the decision
    threshold -- use those two parts only. The test table is scored once at the
    end with the selected weights and the frozen threshold, so the reported test
    metrics are not inflated by having picked anything on the test set.

    Both tables are encoded against the **training** vocabulary, so their token
    indices refer to the same OTUs; test OTUs unseen in training fall back to
    ``'<unk>'``.

    Parameters
    ----------
    metadata : str
        Path to a TSV metadata file. Must contain `sample_id_col` and
        `labels_col`, and its sample IDs must match the BIOM tables.
    train_biom : str
        Path to the training BIOM file (features x samples).
    test_biom : str
        Path to the test BIOM file, in the same format as `train_biom`.
    embedding_birnn : str
        Output path for the best model's ``state_dict``.
    plotfile_loss : str
        Output path for the loss curve plot, e.g. ``'loss.png'``.
    plotfile_auc : str
        Output path for the train/validation AUC curve plot.
    labels_col : str, optional
        Metadata column holding the labels. Default is ``'group'``.
    sample_id_col : str, optional
        Metadata column holding the sample IDs. Default is ``'sample_id'``.
    num_steps : int, optional
        Number of OTU positions kept per sample. Default is 400.
    p_drop : float, optional
        Dropout probability in the transformer layers. Default is 0.
    batch_size : int, optional
        Batch size for both loaders. Default is 128.
    d_model : int, optional
        Transformer embedding dimension. Default is 100.
    n_layers : int, optional
        Number of encoder layers. Default is 2.
    n_heads : int, optional
        Number of attention heads per layer. Default is 2.
    numb : int, optional
        Zero-based CUDA device index. Default is 1.
    lr : float, optional
        Learning rate for the Adam optimizer. Default is 0.0005.
    weight_decay : float, optional
        L2 regularization strength. Default is 0.
    num_epochs : int, optional
        Number of training epochs. Default is 100.
    loss : {'BCE_loss', 'FocalLoss'}, optional
        Which criterion to build. Default is ``'BCE_loss'``.
    alpha : float, optional
        Positive-class weight, used only by ``'FocalLoss'``. Default is 0.6.
    glove_embedding : str, optional
        Path to a pretrained embedding file used to initialize the OTU
        embedding table; its vector width must equal `d_model`. Default is
        None, which fills the table with fixed random codes instead.
    valid_ratio : float, optional
        Fraction of the training table held out for validation. Default is 0.2.
    split_seed : int, optional
        Seed for the train/validation split, independent of the model seed.
        Default is 11.
    group_col : str, optional
        Metadata column holding a grouping key, e.g. a subject ID. Default is
        None. When given, the split keeps all samples of a group on one side,
        which is what you need if one subject contributed several samples.
    pred_out : str, optional
        Prefix for the prediction files, written as ``<pred_out>_valid.csv``
        and ``<pred_out>_test.csv``. Default is None, which derives the prefix
        from `embedding_birnn`.

    Returns
    -------
    valid_record : dict
        The record returned by :func:`train_cls` for the selected epoch.
    test_metrics : dict
        Test-set metrics at the frozen threshold, with keys ``auc``, ``aupr``,
        ``cm``, ``f1``, ``mcc``, ``acc``.

    Notes
    -----
    The embedding table is frozen either way, so what it is filled with decides
    what the OTU identities can contribute. Left at zero it would make every
    sample encode identically and the classifier could only learn a constant,
    so `glove_embedding=None` fills it with fixed random codes instead: the
    taxa stay distinguishable, they just carry no learned relationships to each
    other. The random draw comes from the seed set at the top of this function,
    so a rerun reproduces it.
    """
    set_seed(11)
    devices = try_all_gpus(numb)

    # if glove_embedding is not None:
    #     glove_embedding = f"{glove_embedding}_{d_model}.txt"

    full_train = load_data_imdb(train_biom,
                                metadata,
                                labels_col,
                                sample_id_col,
                                num_steps)
    fid_dict = full_train()
    # Encode the test table with the training vocabulary, so both sides share
    # one index space and the model's embedding rows mean the same OTU.
    test_data = load_data_imdb(test_biom,
                               metadata,
                               labels_col,
                               sample_id_col,
                               num_steps,
                               fid=fid_dict)

    train_labels = to_numpy(full_train.labels)
    groups = None
    if group_col is not None:
        group_map = pd.read_csv(metadata, sep="\t", index_col=sample_id_col,
                                dtype={sample_id_col: str}, low_memory=False)
        groups = group_map.loc[full_train.sample_ids, group_col].to_numpy()
    train_idx, valid_idx = split_train_valid(train_labels, valid_ratio,
                                             split_seed, groups)

    train_iter = DataLoader(Subset(full_train, train_idx),
                            batch_size=batch_size,
                            shuffle=True)
    valid_iter = DataLoader(Subset(full_train, valid_idx),
                            batch_size=batch_size,
                            shuffle=False)
    test_iter = DataLoader(test_data,
                           batch_size=batch_size,
                           shuffle=False)
    print(f"[split] train {len(train_idx)}, valid {len(valid_idx)}, "
          f"test {len(test_data)}")

    # How much of the test cohort the training vocabulary actually covers. These
    # positions are masked out of attention and pooling, so a high share here
    # means the test samples are being read from fewer taxa than the training
    # ones, and the metrics should be read in that light.
    observed = int((test_data.features != fid_dict['<pad>']).sum())
    unknown = int((test_data.features == fid_dict.unk).sum())
    print(f"[vocab] test OTU positions outside the training vocabulary: "
          f"{unknown}/{observed} "
          f"({100 * unknown / max(observed, 1):.1f}%), masked out")

    net = OtuAttentionEncoder(otu_size=len(fid_dict),
                              d_model=d_model,
                              n_layers=n_layers,
                              n_heads=n_heads,
                              p_drop=p_drop,
                              pad_id=fid_dict['<pad>'])
    init_transformer_weights(net)

    if glove_embedding is not None:
        glove_embedding = otuEmbedding(glove_embedding)
        embeds = glove_embedding[fid_dict.idx_to_token]
        if embeds.shape[1] != d_model:
            raise ValueError(
                f"the embedding file supplies {embeds.shape[1]}-dimensional "
                f"vectors but the model is built with d_model={d_model}")
        embeds = embeds.float()
    else:
        # Without pretrained vectors, give each OTU a fixed random code rather
        # than leaving the frozen table at zero -- zeros would erase the input
        # entirely and let the model do no better than predicting a constant.
        embeds = torch.randn(len(fid_dict), d_model) * (d_model ** -0.5)
    # Both markers stay at zero: '<pad>' is what padding_idx expects, and
    # '<unk>' is masked out anyway, so it should not add a spurious direction.
    embeds[fid_dict['<pad>']] = 0
    embeds[fid_dict.unk] = 0
    net.embedding.weight.data.copy_(embeds)
    net.embedding.weight.requires_grad = False

    trainer = torch.optim.Adam(
        net.parameters(), lr=lr, weight_decay=weight_decay)
    scheduler = torch.optim.lr_scheduler.ExponentialLR(trainer, gamma=0.9)
    if loss == "FocalLoss":
        loss = FocalLoss(alpha=alpha, devices=devices)
    elif loss == "BCE_loss" or loss is None:
        loss = nn.BCELoss(weight=None, reduction='mean')
    else:
        raise ValueError(f"unsupported loss {loss!r}; expected 'BCE_loss' or "
                         f"'FocalLoss'")

    valid_record = train_cls(net, train_iter, valid_iter, loss, trainer,
                             scheduler, num_epochs, devices, embedding_birnn,
                             plotfile_loss, plotfile_auc)

    threshold = valid_record['threshold']
    test_metrics, test_prob, test_label = evaluate_on_test(
        net, embedding_birnn, test_iter, threshold, devices)
    print(f"[test ] auc {test_metrics['auc']:.3f}, aupr {test_metrics['aupr']:.3f}, "
          f"f1 {test_metrics['f1']:.3f}, mcc {test_metrics['mcc']:.3f}, "
          f"acc {test_metrics['acc']:.3f} | threshold {threshold:.4f} (from valid)")
    print(f"[test ] confusion_matrix\n{test_metrics['cm']}")

    if pred_out is None:
        pred_out = os.path.splitext(embedding_birnn)[0]
    save_predictions(f"{pred_out}_valid.csv",
                     full_train.sample_ids[valid_idx],
                     valid_record['valid_label'],
                     valid_record['valid_prob'],
                     threshold)
    save_predictions(f"{pred_out}_test.csv",
                     test_data.sample_ids,
                     test_label,
                     test_prob,
                     threshold)
    print(f"[out  ] predictions written to {pred_out}_valid.csv and "
          f"{pred_out}_test.csv")

    return valid_record, test_metrics