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
from torch.utils.data import Dataset, DataLoader
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


def read_strata(metadata, sample_id_col, strata_col, sample_ids):
    """Read one categorical metadata column, aligned to a table's samples.

    Used to recover the disease each training sample belongs to, which the OTU
    table itself does not carry. The returned order follows `sample_ids`, i.e.
    the row order of the tensors :class:`load_data_imdb` built from the same
    BIOM table, so an index into one indexes the other.

    This is kept out of :func:`read_imdb` on purpose. The disease is a
    *sampling* concern, not part of the "BIOM table -> padded tensors"
    contract, and of the three tables a run loads only the training one needs
    it. Reading it here also means the default path never touches this column.

    Parameters
    ----------
    metadata : str
        Path to a tab-separated metadata file.
    sample_id_col : str
        Metadata column holding sample IDs; must match the BIOM sample IDs.
    strata_col : str
        Metadata column naming the stratum, e.g. ``'disease_name_ab'``.
    sample_ids : array_like of str
        Sample IDs to read, in the order the caller wants them back. Pass
        ``dataset.sample_ids``.

    Returns
    -------
    numpy.ndarray
        Stratum labels as strings, aligned with `sample_ids`.

    Raises
    ------
    KeyError
        If `metadata` lacks `sample_id_col` or `strata_col`, or does not cover
        every ID in `sample_ids`.
    ValueError
        If `sample_id_col` holds duplicate IDs, or if `strata_col` is empty or
        missing for any of `sample_ids`.

    Notes
    -----
    A blank stratum raises rather than falling back to an ``'<unknown>'``
    bucket. That bucket would become a disease of its own and, being small,
    could end up carrying a full vote in a macro average or a base rate of its
    own in the logit offsets -- both decided by a few empty cells. This is the
    same strictness :func:`read_imdb` applies to `labels_col`.
    """
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

    if strata_col not in mapping_file.columns:
        raise KeyError(
            f"{metadata} has no column {strata_col!r}; available columns: "
            f"{list(mapping_file.columns)}")

    if not mapping_file.index.is_unique:
        duplicated = mapping_file.index[mapping_file.index.duplicated()].unique()
        raise ValueError(
            f"{sample_id_col!r} holds duplicate IDs in {metadata}, which would "
            f"misalign the strata; first offenders: {list(duplicated[:5])}")

    sample_ids = np.asarray(sample_ids)
    missing = [sid for sid in sample_ids if sid not in mapping_file.index]
    if missing:
        raise KeyError(
            f"{len(missing)} of {len(sample_ids)} samples have no row in "
            f"{metadata}; first offenders: {missing[:5]}")

    values = mapping_file.loc[sample_ids, strata_col]
    blank = values.isna() | (values.astype(str).str.strip() == "")
    if blank.any():
        offenders = values.index[blank][:5].tolist()
        raise ValueError(
            f"column {strata_col!r} is empty for {int(blank.sum())} of "
            f"{len(sample_ids)} samples; fill it, or point `disease_col` at a "
            f"different column. Offending samples: {offenders}")

    return values.astype(str).to_numpy()


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
        # Zero unless a caller installs one. Carried by the dataset rather
        # than looked up in the training loop because the loader shuffles and
        # resamples, so the offset has to travel with its own sample.
        self.logit_offset = torch.zeros(len(labels), dtype=torch.float)

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

    def set_logit_offset(self, offsets):
        """Install a per-sample logit offset, aligned with the table's rows.

        Parameters
        ----------
        offsets : array_like
            One float per sample, in the order :attr:`sample_ids` holds. See
            :func:`logit_adjust_offsets`.

        Raises
        ------
        ValueError
            If `offsets` is not one value per sample.
        """
        offsets = np.asarray(offsets, dtype=float).ravel()
        if len(offsets) != len(self.labels):
            raise ValueError(
                f"offsets has {len(offsets)} entries but the table holds "
                f"{len(self.labels)} samples")
        self.logit_offset = torch.tensor(offsets, dtype=torch.float)

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
        logit_offset : torch.Tensor
            Scalar logit offset for the sample, 0 unless
            :meth:`set_logit_offset` installed one.
        """
        return (self.features[index, ], self.abundance[index, ],
                self.labels[index], self.mask[index],
                self.logit_offset[index])

    def __len__(self):
        """Return the number of samples.

        Returns
        -------
        int
            Number of samples in the dataset.
        """
        return self.otu.shape[0]


def logit_adjust_offsets(strata, labels, tau=1.0, verbose=True):
    """Per-sample logit offsets from each disease's own case/control base rate.

    Logit adjustment (Menon et al., ICLR 2021) trains on ``f(x) + tau * log
    pi_y`` and predicts with ``f(x)`` alone, so the base rate lives in the
    offset instead of in the learned scores. For a single-logit binary model,
    where ``z = f_1 - f_0``, the pair of class offsets collapses to one number
    per sample::

        z_adjusted = z + tau * log(pi_1 / pi_0)

    Taking the long-tailed label to be the joint ``(disease, class)`` cell
    makes ``pi_1 / pi_0`` conditional on the disease, and the disease's own
    size cancels out of the ratio::

        (n_d1 / N) / (n_d0 / N) = n_d1 / n_d0

    which is worth being explicit about: this corrects each disease's
    **base rate**, not the fact that one disease has sixty times the samples
    of another. That size imbalance survives untouched in the gradient
    weighting; the offsets address a different thing.

    What it buys here is that ``z`` stops carrying "cohorts like this one are
    usually cases". Under a left-out-disease split the test cohort's base rate
    is unknown and unlike any in training, so a scorer that has that term
    baked in transfers a bias it cannot justify -- and a threshold fitted on
    the training diseases misses on the test one, which is what a high AUC
    beside a poor accuracy looks like.

    Parameters
    ----------
    strata : array_like
        Disease of each training sample. See :func:`read_strata`.
    labels : array_like
        0/1 class label of each training sample, same order.
    tau : float, optional
        Strength of the correction. Default is 1.0, which removes the base
        rate exactly; 0 disables it and reproduces plain training. Values
        between are the usual sweep.
    verbose : bool, optional
        Print the per-disease table. Default is True.

    Returns
    -------
    offsets : numpy.ndarray
        Float offset per sample, aligned with `strata`, added to the logit
        during training only.
    per_disease : dict
        ``{disease: offset}``, for the record.

    Raises
    ------
    ValueError
        If `strata` and `labels` differ in length, or `tau` is negative or not
        finite.

    Notes
    -----
    A disease holding one class has an infinite log-ratio. Its offset is
    clamped to the largest finite one seen, rather than left infinite -- an
    infinite offset saturates the loss for every sample of that disease and
    stops their gradients entirely, which removes the disease from training
    instead of correcting it.

    Examples
    --------
    >>> import numpy as np
    >>> s = np.array(['A', 'A', 'A', 'B', 'B', 'B'])
    >>> y = np.array([1, 1, 0, 1, 0, 0])
    >>> off, per = logit_adjust_offsets(s, y, verbose=False)
    >>> round(per['A'], 4), round(per['B'], 4)
    (0.6931, -0.6931)
    """
    strata = np.asarray([str(s) for s in strata])
    labels = np.asarray(labels).astype(int).ravel()
    if len(strata) != len(labels):
        raise ValueError(
            f"strata has {len(strata)} entries but labels has {len(labels)}; "
            f"they must be aligned with the same dataset")
    if not np.isfinite(tau):
        raise ValueError(f"tau must be finite, got {tau!r}")
    if tau < 0:
        raise ValueError(
            f"tau must be >= 0, got {tau!r}; a negative value would amplify "
            f"each disease's base rate instead of removing it")

    counts, raw = {}, {}
    for d in sorted(set(strata)):
        m = strata == d
        n_pos = int((labels[m] == 1).sum())
        n_neg = int((labels[m] == 0).sum())
        counts[d] = (n_pos, n_neg)
        raw[d] = (np.inf if n_neg == 0 else
                  -np.inf if n_pos == 0 else
                  float(np.log(n_pos / n_neg)))

    finite = [v for v in raw.values() if np.isfinite(v)]
    cap = max((abs(v) for v in finite), default=1.0) or 1.0
    per_disease, clamped = {}, []
    for d, v in raw.items():
        if not np.isfinite(v):
            v = cap if v > 0 else -cap
            clamped.append(d)
        per_disease[d] = float(tau * v)

    offsets = np.array([per_disease[d] for d in strata], dtype=float)

    if verbose:
        print(f"[logitadj] tau {tau:g}, per-disease case/control base rate "
              f"removed from the training logits ({len(per_disease)} diseases)")
        order = sorted(per_disease, key=lambda d: -abs(per_disease[d]))
        width = max(len(str(d)) for d in order)
        for d in order[:12]:
            n_pos, n_neg = counts[d]
            star = " *" if d in clamped else ""
            print(f"[logitadj] {str(d):<{width}} case {n_pos:5d} / ctrl "
                  f"{n_neg:5d} -> offset {per_disease[d]:+.3f}{star}")
        if len(order) > 12:
            print(f"[logitadj] ... and {len(order) - 12} more")
        if clamped:
            print(f"[logitadj] * {', '.join(clamped)} hold one class only; "
                  f"their offset is clamped to +/-{tau * cap:.3f} rather than "
                  f"infinite, which would zero their gradients")
        print(f"[logitadj] validation and test see raw logits -- the offset is "
              f"a training-time device, so the scores it selects and reports "
              f"on are already free of the base rate")
    return offsets, per_disease


class LogitAdjustedLoss(nn.Module):
    """Binary cross-entropy on logits shifted by a per-sample prior offset.

    The offset is added inside the loss rather than by the model, so the
    network's output stays the prior-free score at every point it is read --
    validation, the threshold, the test predictions. Only the gradient sees
    the shifted value. See :func:`logit_adjust_offsets`.

    Parameters
    ----------
    pos_weight : torch.Tensor, optional
        Passed through to the underlying criterion. Default is None. It is
        rarely wanted alongside this loss: both correct a class prior, and the
        offset already does it per disease.
    reduction : {'mean', 'sum'}, optional
        How the per-sample losses are aggregated. Default is ``'mean'``.
    """

    def __init__(self, pos_weight=None, reduction='mean'):
        """Store the criterion configuration."""
        super().__init__()
        self.pos_weight = pos_weight
        self.reduction = reduction

    def forward(self, logits, targets, offset=None):
        """Loss on ``logits + offset``; `offset` of None means plain BCE."""
        if offset is not None:
            logits = logits + offset.to(logits.device, dtype=logits.dtype)
        return nn.functional.binary_cross_entropy_with_logits(
            logits.float(), targets.float(), pos_weight=self.pos_weight,
            reduction=self.reduction)


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
    """Transformer encoder layer: self-attention followed by a feed-forward net.

    Two residual sub-layers, each in post-norm form. The first is multi-head
    self-attention, which mixes information *across* OTUs. The second is a
    position-wise feed-forward network, which applies a nonlinearity *within*
    each OTU position independently.

    Both sub-layers are needed for the layer to be more than a linear map.
    Attention alone reweights and averages value vectors; every operation it
    applies to the representation is linear apart from the softmax over the
    weights, so a stack of attention-only layers followed by mean pooling and a
    linear head collapses to something very close to a linear model on the
    input embeddings. The feed-forward sub-layer is what gives each layer its
    per-OTU nonlinearity.

    Parameters
    ----------
    d_model : int
        Model dimension.
    n_heads : int
        Number of attention heads.
    p_drop : float
        Dropout probability, applied to both sub-layers and inside the
        feed-forward network.
    d_ff : int, optional
        Width of the feed-forward hidden layer. Default is None, which uses
        ``4 * d_model``, the usual transformer ratio. Pass 0 to drop the
        feed-forward sub-layer entirely, which reduces the layer to
        attention-only and is there to run as the ablation.

    Attributes
    ----------
    mha : MultiHeadAttention
        The self-attention sub-layer.
    ffn : torch.nn.Sequential or None
        The position-wise feed-forward sub-layer,
        ``Linear -> GELU -> Dropout -> Linear``, or None when ``d_ff`` is 0.
    """

    def __init__(self, d_model, n_heads, p_drop, d_ff=None):
        """Build the attention and feed-forward sub-layers with their norms."""
        super(EncoderLayer, self).__init__()
        d_ff = 4 * d_model if d_ff is None else d_ff
        self.d_ff = d_ff
        self.mha = MultiHeadAttention(d_model, n_heads)
        self.dropout1 = nn.Dropout(p_drop)
        self.layernorm1 = nn.LayerNorm(d_model, eps=1e-6)
        if d_ff:
            self.ffn = nn.Sequential(
                nn.Linear(d_model, d_ff),
                nn.GELU(),
                nn.Dropout(p_drop),
                nn.Linear(d_ff, d_model),
            )
            self.dropout2 = nn.Dropout(p_drop)
            self.layernorm2 = nn.LayerNorm(d_model, eps=1e-6)
        else:
            self.ffn = None

    def forward(self, inputs, attn_mask):
        """Run the attention and feed-forward sub-layers, each with a residual.

        Parameters
        ----------
        inputs : torch.Tensor
            Input of shape ``(batch, seq_len, d_model)``.
        attn_mask : torch.Tensor
            Boolean mask of shape ``(batch, seq_len, seq_len)``, where ``True``
            marks a position to suppress.

        Returns
        -------
        outputs : torch.Tensor
            Normalized output of shape ``(batch, seq_len, d_model)``.
        attn_weights : torch.Tensor
            Attention weights of shape ``(batch, n_heads, seq_len, seq_len)``.

        Notes
        -----
        The feed-forward sub-layer runs on masked-out positions too, since it
        acts on each position independently and cannot leak anything between
        them. Those positions are dropped later, at the pooling step, so what
        the network computes there never reaches the output.
        """
        attn_outputs, attn_weights = self.mha(inputs, inputs, inputs,
                                              attn_mask)
        attn_outputs = self.dropout1(attn_outputs)
        attn_outputs = self.layernorm1(inputs + attn_outputs)

        if self.ffn is None:
            return attn_outputs, attn_weights

        ffn_outputs = self.ffn(attn_outputs)
        ffn_outputs = self.dropout2(ffn_outputs)
        outputs = self.layernorm2(attn_outputs + ffn_outputs)
        return outputs, attn_weights


class OtuAttentionEncoder(nn.Module):
    """Attention encoder over abundance-weighted OTU sets.

    Each sample is treated as a *set* of OTUs rather than a sequence. OTU
    embeddings are scaled by their relative abundance, refined by a stack of
    :class:`EncoderLayer`, mean-pooled over the unmasked positions into a
    sample representation, and mapped to a single output by an MLP head.

    This departs from a standard transformer encoder in one way: there is no
    positional encoding, so the model is invariant to the order in which the
    OTUs are listed. The encoder layers themselves are standard, each pairing
    self-attention with a position-wise feed-forward sub-layer.

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
    d_ff : int, optional
        Width of the feed-forward hidden layer inside each encoder layer.
        Default is None, which uses ``4 * d_model``.
    head_hidden : int, optional
        Width of the hidden layer in the output head. Default is None, which
        uses ``d_model // 2``. Pass 0 to fall back to a single linear layer.
    abund_mode : {'multiply', 'none'}, optional
        How abundance enters the input representation. ``'multiply'``, the
        default, scales each OTU embedding by its abundance. ``'none'``
        ignores abundance entirely, reducing the input to the presence or
        absence of each covered taxon.

    Attributes
    ----------
    embedding : torch.nn.Embedding
        Token embedding table of shape ``(otu_size, d_model)``.
    layers : torch.nn.ModuleList
        The `n_layers` encoder layers.
    mlp : torch.nn.Sequential
        Output head mapping the pooled embedding to a scalar logit, with one
        hidden layer unless `head_hidden` is 0.

    Notes
    -----
    The head reads a single pooled vector, so a bare linear layer there makes
    the last thing the model does a linear function of the community
    representation. One hidden layer with a nonlinearity lets the head
    respond to combinations of the pooled dimensions -- which is where the
    interactions between taxa that the encoder assembled actually get used.

    Abundance is applied as a scale rather than an added vector, and the
    scale does not survive the first :class:`torch.nn.LayerNorm`, which
    normalizes each position over `d_model`: ``LayerNorm(c * x)`` equals
    ``LayerNorm(x)`` exactly, for any positive ``c``. That is deliberate.
    Under the rank-over-max normalization the abundance at a position is
    close to a function of the position's index and the sample's richness, so
    a channel that carried it faithfully would mostly carry richness -- a
    study-level batch effect. Dropping the magnitude and keeping the
    direction leaves the pretrained OTU codes, whose information is entirely
    directional, uncontaminated. Abundance still reaches the model through
    the attention logits, which scale with it.
    """

    def __init__(self,
                 otu_size,
                 d_model=128,
                 n_layers=6,
                 n_heads=8,
                 p_drop=0.1,
                 pad_id=0,
                 d_ff=None,
                 head_hidden=None,
                 abund_mode='multiply',
                 linear_branch=False):
        """Build the embedding table, the encoder stack, and the output head."""
        super(OtuAttentionEncoder, self).__init__()
        if abund_mode not in ('multiply', 'none'):
            raise ValueError(f"unsupported abund_mode {abund_mode!r}; expected "
                             f"'multiply' or 'none'")
        self.abund_mode = abund_mode

        self.embedding = nn.Embedding(otu_size, d_model, padding_idx=pad_id)
        self.layers = nn.ModuleList([
            EncoderLayer(d_model, n_heads, p_drop, d_ff)
            for _ in range(n_layers)
        ])
        self.sigmoid = nn.Sigmoid()
        self.pad_id = pad_id
        self.dropout = nn.Dropout(p=p_drop)
        head_hidden = max(d_model // 2, 1) if head_hidden is None else head_hidden
        if head_hidden:
            self.mlp = nn.Sequential(
                nn.Linear(d_model, head_hidden),
                nn.GELU(),
                nn.Dropout(p_drop),
                nn.Linear(head_hidden, 1),
            )
        else:
            self.mlp = nn.Sequential(
                nn.Linear(d_model, 1),
            )

        # Linear residual branch: one scalar coefficient per OTU, applied to
        # the abundance directly. This is the pathway the pooled encoder cannot
        # express -- mean pooling collapses the per-OTU values into a single
        # d_model vector, so a coefficient attached to an individual taxon has
        # nowhere to live. A linear probe on the raw abundances reaches 0.767
        # on IBD where the pooled representation reaches 0.733, and that gap is
        # what this branch exists to close.
        #
        # A plain Parameter rather than nn.Embedding: Embedding draws its rows
        # from the RNG on construction, which would shift the stream for every
        # layer initialized afterwards and stop linear_branch=True/False from
        # being a clean ablation. torch.zeros consumes no random numbers.
        self.linear_branch = linear_branch
        if linear_branch:
            self.lin_w = nn.Parameter(torch.zeros(otu_size))
            self.lin_b = nn.Parameter(torch.zeros(1))

    def forward(self, x, weight, mask, encoder=False):
        """Encode a batch of samples and produce a scalar logit per sample.

        Parameters
        ----------
        x : torch.Tensor
            Token indices of shape ``(batch, seq_len)``.
        weight : torch.Tensor
            Abundance weights of shape ``(batch, seq_len)``, multiplied into
            the token embeddings unless `abund_mode` is ``'none'``.
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

        Under ``abund_mode='none'`` the scaling is simply skipped, which is
        equivalent to an abundance of 1 everywhere. Padding and ``'<unk>'``
        positions stay at zero regardless, since their embedding rows are
        zero.
        """
        # The abundance arrives as float64 from the dataset. Cast it once here
        # rather than letting it promote the embeddings to float64 and
        # materializing the whole batch at double width.
        weight = weight.float()
        mask_f = mask.unsqueeze(-1).float()

        inputs = self.embedding(x)
        if self.abund_mode == 'multiply':
            inputs = inputs * weight.unsqueeze(-1)
        inputs_1 = inputs.clone()
        attn_pad_mask = self.get_attention_padding_mask(mask)

        for layer in self.layers:
            inputs, attn_weights = layer(inputs, attn_pad_mask)

        embedding_sum = inputs * mask_f
        # A sample can in principle carry no covered OTU at all; clamping keeps
        # the pooling finite instead of handing NaNs to the loss.
        n_valid = mask.sum(1, keepdim=True).clamp(min=1)
        embedding = embedding_sum.sum(dim=1, keepdim=False) / n_valid
        outputs = self.dropout(embedding)
        outputs = self.mlp(outputs)

        if self.linear_branch:
            # Abundance enters here at full magnitude whatever abund_mode says:
            # this branch exists precisely to carry the per-taxon abundance the
            # LayerNorm strips out of the attention path. Set
            # linear_branch=False to ablate it.
            #
            # The positions are summed, matching a logistic regression on the
            # abundance vector. If training turns out unstable, divide by
            # `n_valid` instead -- that makes the term a mean and so invariant
            # to sample richness, which is a study-level batch effect.
            z = self.lin_w[x] * weight
            lin = (z * mask.float()).sum(1, keepdim=True)
            outputs = outputs + lin + self.lin_b

        if encoder:
            return inputs_1, embedding
        return outputs, attn_weights

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


def macro_auc(labels, prob, strata, min_group_n=0):
    """Mean of the per-group AUCs, one group per disease.

    A pooled AUC over a validation set spanning several diseases is not the
    quantity this task is after. AUC counts case/control pairs, and pooling
    counts the *cross-disease* pairs too -- whether a CRC case outranks a T2DM
    control. Nothing asks that question: every test cohort here is scored
    within one disease. Worse, those pairs are the majority. With K diseases of
    similar size only about ``1 / K`` of the pooled pairs are within-disease,
    so at a dozen diseases the pooled number is mostly reporting whether the
    model can tell the diseases apart -- a between-disease offset in the scores
    -- rather than whether it separates cases from controls inside any of them.
    A model with no within-disease signal at all can score well on it.

    Averaging the per-disease AUCs discards the cross-disease pairs entirely
    and gives every disease one vote, whatever its sample count.

    `min_group_n` exists because one vote each equalises the *level* but not
    the *variance*. An AUC's sampling error goes as ``1 / sqrt(n)``, so a
    disease with 24 validation samples swings by around 0.13 from epoch to
    epoch on noise alone while one with 300 swings by 0.04. In an unweighted
    mean the small group therefore supplies most of the epoch-to-epoch
    movement, and whichever epoch it happens to like is the epoch the average
    picks -- the size bias comes back, pointing the other way. Dropping the
    groups too small to estimate an AUC on is the blunt fix, and the one whose
    effect can be read off the output.

    Parameters
    ----------
    labels : array_like
        True binary labels of shape ``(n_samples,)``.
    prob : array_like
        Predicted probabilities of the positive class, shape ``(n_samples,)``.
    strata : array_like
        Group of each sample -- normally the disease -- same length and order
        as `labels`. See :func:`read_strata`.
    min_group_n : int, optional
        Groups with fewer than this many samples are left out, their AUC being
        too noisy to select an epoch on. Default is 0, which keeps every group
        that holds both classes. Note the sampling error is really driven by
        the *smaller* class, so a group of 100 that is 97 cases and 3 controls
        is as noisy as a balanced group of 6 and passes a threshold on the
        total; the per-group class counts are reported alongside so that case
        is visible.

    Returns
    -------
    value : float
        Unweighted mean of the per-group AUCs, NaN when no group could be
        scored.
    per_group : dict
        ``{group: auc}`` for the groups that could be scored, in sorted order.
    skipped : dict
        ``{group: reason}`` for the groups that were left out -- either they
        hold one class, or they fall under `min_group_n`. A group holding one
        class has no case/control pair to rank and no AUC to contribute; it is
        left out rather than counted as 0.5, which would be a claim about a
        model that was never tested there.

    Examples
    --------
    >>> import numpy as np
    >>> y = np.array([0, 1, 0, 1])
    >>> p = np.array([0.1, 0.9, 0.8, 0.2])
    >>> s = np.array(['A', 'A', 'B', 'B'])
    >>> value, per, skipped = macro_auc(y, p, s)
    >>> per['A'], per['B'], value
    (1.0, 0.0, 0.5)
    """
    labels = np.asarray(labels)
    prob = np.asarray(prob)
    strata = np.asarray([str(s) for s in strata])
    if not (len(labels) == len(prob) == len(strata)):
        raise ValueError(
            f"labels ({len(labels)}), prob ({len(prob)}) and strata "
            f"({len(strata)}) must be the same length and in the same order")

    per_group, skipped = {}, {}
    for g in sorted(set(strata)):
        m = strata == g
        y = labels[m]
        n = int(m.sum())
        n_pos, n_neg = int((y == 1).sum()), int((y == 0).sum())
        classes = np.unique(y)
        if len(classes) < 2:
            only = "cases" if classes[0] == 1 else "controls"
            skipped[g] = f"n={n}, all {only}"
            continue
        if n < min_group_n:
            skipped[g] = (f"n={n} ({n_pos}/{n_neg}) < min_group_n "
                          f"{min_group_n}, too noisy to select on")
            continue
        per_group[g] = float(metrics.roc_auc_score(y, prob[m]))

    value = float(np.mean(list(per_group.values()))) if per_group else float('nan')
    return value, per_group, skipped


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


def criterion_value(loss, logits, targets, offset=None):
    """Evaluate a criterion on whichever of logits or probabilities it expects.

    :class:`torch.nn.BCEWithLogitsLoss` folds the sigmoid into the loss and
    must be handed raw logits; :class:`torch.nn.BCELoss` and :class:`FocalLoss`
    both consume probabilities. This picks the right one so the training and
    evaluation paths do not each need to know which criterion is in play.

    Parameters
    ----------
    loss : torch.nn.Module
        The criterion.
    logits : torch.Tensor
        Raw model outputs of shape ``(batch,)``.
    targets : torch.Tensor
        Binary targets of shape ``(batch,)``, as floats.
    offset : torch.Tensor, optional
        Per-sample logit offset of shape ``(batch,)``, used only by
        :class:`LogitAdjustedLoss` and ignored by every other criterion.
        Default is None.

    Returns
    -------
    torch.Tensor
        Scalar loss.

    Notes
    -----
    Handing logits to ``BCEWithLogitsLoss`` is not merely tidier than applying
    the sigmoid first. The fused form evaluates the loss through
    ``log(1 + exp(-|x|))``, which stays finite at any logit, whereas a
    separate sigmoid saturates to exactly 0 or 1 in float32 and leaves
    ``BCELoss`` to clamp its logarithm at -100 -- a clamp that returns a
    gradient of zero for the samples the model is most wrong about.
    """
    if isinstance(loss, LogitAdjustedLoss):
        return loss(logits, targets, offset)
    if isinstance(loss, nn.BCEWithLogitsLoss):
        return loss(logits.float(), targets)
    return loss(torch.sigmoid(logits).float(), targets)


def train_batch_ch13(net, X, y, abundance, loss, mask, trainer, devices,
                     offset=None):
    """Run one training step on a minibatch.

    The model output is squeezed to shape ``(batch,)``. Whether the criterion
    sees those logits directly or their sigmoid is decided by
    :func:`criterion_value`; the probabilities are returned either way, since
    the caller needs them for the running AUC.

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
        Loss criterion, on either logits or probabilities.
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
    prob : torch.Tensor
        Predicted probabilities of shape ``(batch,)``, detached.
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
    logits, _ = net(X, abundance, mask)
    logits = torch.squeeze(logits, dim=1)
    l = criterion_value(loss, logits, y.float(),
                        offset=None if offset is None
                        else offset.to(devices[0]))
    l.backward()
    trainer.step()
    train_loss_sum = l.sum().detach()
    return train_loss_sum, torch.sigmoid(logits).detach()


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
        for features, abundance, labels, mask, offset in data_iter:
            X = features.to(devices[0])
            abundance = abundance.to(devices[0])
            mask = mask.to(devices[0])
            pred, _ = net(X, abundance, mask)
            logits = torch.squeeze(pred, dim=1)
            prob = torch.sigmoid(logits)
            if loss is not None:
                y = labels.to(devices[0], dtype=torch.int64)
                # No offset here on purpose: the held-out cohorts are scored
                # on the prior-free logits, which is the whole point of
                # adjusting during training only.
                total_loss += float(criterion_value(loss, logits, y.float()))
            probs.append(to_numpy(prob))
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


def train_cls(net, train_iter, valid_iter, loss, trainer, num_epochs,
              devices, ckpt_path, plotfile_loss, plotfile_auc,
              patience=10, min_delta=0.001, select_by='loss',
              valid_strata=None, min_group_n=0):
    """Train a binary classifier, selecting the best epoch by validation loss.

    The test set is deliberately absent from this function's signature. Every
    decision made here -- which epoch to keep and which decision threshold to
    use -- is made on `valid_iter` alone, so the held-out test set stays
    untouched until :func:`evaluate_on_test` runs afterwards.

    After each epoch the validation loss is computed; when it drops below the
    best seen so far by at least `min_delta`, the ``state_dict`` is written to
    `ckpt_path` and the optimal threshold **of that same epoch** is recorded,
    so the returned threshold always matches the returned weights.

    Parameters
    ----------
    net : torch.nn.Module
        Model to train.
    train_iter : torch.utils.data.DataLoader
        Loader yielding ``(features, abundance, labels, mask)`` batches, with
        `features`, `abundance`, and `mask` of shape ``(batch, seq_len)`` and
        `labels` of shape ``(batch,)`` holding 0/1 class labels.
    valid_iter : torch.utils.data.DataLoader
        Validation loader with the same structure as `train_iter`.
    loss : torch.nn.Module
        Loss criterion, on either logits or probabilities; see
        :func:`criterion_value`.
    trainer : torch.optim.Optimizer
        Optimizer used for parameter updates.
    num_epochs : int
        Maximum number of training epochs.
    devices : list of torch.device
        Devices to use; only ``devices[0]`` is used.
    ckpt_path : str
        Path the selected ``state_dict`` is saved to.
    plotfile_loss : str
        Path for the train/validation loss curve plot.
    plotfile_auc : str
        Path for the train/validation AUC curve plot.
    patience : int, optional
        Stop once the validation loss has not improved by at least `min_delta`
        for this many epochs. Default is 10.
    min_delta : float, optional
        Minimum absolute decrease in validation loss that counts as an
        improvement. Default is 0.001.
    valid_strata : array_like, optional
        Group of each validation sample -- normally its disease -- aligned with
        the order `valid_iter` yields, which is the dataset's own order since
        that loader does not shuffle. Default is None, which scores validation
        with a single pooled AUC. When given, ``valid_auc`` becomes the mean of
        the per-group AUCs instead, so every disease weighs the same and the
        cross-disease pairs a pooled AUC would count are dropped; see
        :func:`macro_auc`. Only the AUC changes -- the loss, and the threshold
        :func:`fit_threshold` picks, stay pooled.

        This changes which epoch is kept **only under** ``select_by='auc'``,
        which is the branch that compares ``valid_auc``. Under
        ``select_by='loss'`` the epoch is chosen on the pooled validation loss
        and the macro AUC is computed and plotted but decides nothing -- and
        the pooled loss is a per-sample mean, so it weights each disease by its
        sample count, which is the imbalance the macro average removes. Pair
        `valid_strata` with ``select_by='auc'``.
    min_group_n : int, optional
        Under `valid_strata`, drop diseases with fewer than this many
        validation samples from the average. Default is 0, which keeps every
        disease holding both classes. One vote each equalises the level but
        not the variance, and a small disease's noise can end up choosing the
        epoch; see :func:`macro_auc`.

    Returns
    -------
    dict
        Record of the selected model, with keys ``epoch``, ``threshold``,
        ``valid_auc``, ``valid_prob``, ``valid_label``, ``train_loss``,
        ``valid_loss``, ``train_auc``, and ``valid_auc_per_disease`` -- the
        last being ``{disease: auc}`` when `valid_strata` was given and None
        otherwise. Note ``train_auc`` stays pooled either way: the training
        part's disease mix is not the validation part's, so a macro number
        there would not be comparable to the one beside it on the plot.

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
    stale = 0

    if select_by not in ('loss', 'auc', 'last'):
        raise ValueError(f"unsupported select_by {select_by!r}; "
                         f"expected 'loss', 'auc', or 'last'")
    if select_by == 'last':
        print(f"[train] select_by='last': training the full {num_epochs} "
              f"epochs and keeping the last, with no early stopping and no "
              f"selection on validation")
        
    for epoch in range(num_epochs):
        train_loss = 0.0
        n_batches = 0
        train_probs, train_labels = [], []
        for features, abundance, labels, mask, offset in train_iter:
            l, pred = train_batch_ch13(
                net, features, labels, abundance, loss, mask, trainer, devices,
                offset=offset)
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
        if valid_strata is None:
            valid_auc = metrics.roc_auc_score(valid_label, valid_prob)
            per_disease = None
        else:
            # One AUC per disease, averaged. The validation set does not move
            # between epochs, so which diseases can be scored is settled on
            # the first epoch; report it once there rather than on all
            # hundred.
            valid_auc, per_disease, skipped = macro_auc(
                valid_label, valid_prob, valid_strata,
                min_group_n=min_group_n)
            if epoch == 0:
                print(f"[macro] validation AUC averaged over "
                      f"{len(per_disease)} disease(s), not pooled over "
                      f"{len(valid_label)} samples")
                # Sizes of the groups that do count, so the threshold can be
                # tuned from one run's output: an AUC's noise goes as
                # 1/sqrt(n), and it is the smaller class that binds.
                sizes = []
                for g in sorted(per_disease):
                    m = np.asarray(valid_strata) == g
                    yg = valid_label[m]
                    sizes.append(f"{g} n={int(m.sum())}"
                                 f"({int((yg == 1).sum())}/"
                                 f"{int((yg == 0).sum())})")
                print("[macro] scored: " + "  ".join(sizes))
                if skipped:
                    print(f"[macro] {len(skipped)} disease(s) left out: "
                          + "; ".join(f"{g} ({r})"
                                      for g, r in sorted(skipped.items())))
                if not per_disease:
                    raise ValueError(
                        f"no disease in the validation set can be scored, so "
                        f"a macro AUC cannot be formed; lower min_group_n "
                        f"(currently {min_group_n}), use valid_auc='pooled', "
                        f"or rebuild the split")
        animator_1.add(epoch + 1, (train_loss, valid_loss), plotfile_loss)
        animator_2.add(epoch + 1, (train_auc, valid_auc), plotfile_auc)

        # Select the epoch with the lowest validation loss. An improvement is
        # counted only when the drop exceeds min_delta, so tiny fluctuations do
        # not reset the early-stopping counter.
        if select_by == 'last':
            # Keep every epoch, so `best` ends up holding the final one. The
            # stale counter never advances either, which turns early stopping
            # off -- deliberately the same switch, because stopping early is
            # only meaningful when the metric being watched is worth watching.
            # When validation AUC does not predict test AUC, an early peak in
            # it is a lucky draw, and letting that draw end training after two
            # epochs is worse than not selecting at all.
            improved = True
        elif select_by == 'loss':
            improved = (best is None
                        or (best['valid_loss'] - valid_loss) > min_delta)
        else:
            improved = (best is None
                        or (valid_auc - best['valid_auc']) > min_delta)

        if improved:
            best = {'epoch': epoch + 1,
                    'threshold': fit_threshold(valid_label, valid_prob),
                    'valid_auc': valid_auc,
                    'valid_prob': valid_prob,
                    'valid_label': valid_label,
                    'train_loss': train_loss,
                    'valid_loss': valid_loss,
                    'train_auc': train_auc,
                    'valid_auc_per_disease': per_disease}
            torch.save(net.state_dict(), ckpt_path)
            stale = 0
        else:
            stale += 1

        # Early stopping: halt when validation loss has not improved by at
        # least min_delta for patience consecutive epochs.
        if patience is not None and stale >= patience:
            print(f"[early] no valid loss improvement >= {min_delta} for "
                  f"{patience} epochs; stopping at epoch {epoch + 1}/{num_epochs} "
                  f"(best was epoch {best['epoch']})")
            break

    v_auc, v_aupr, v_cm, v_f1, v_mcc, v_acc = evaluate_cls(
        best['valid_label'], best['valid_prob'], best['threshold'])
    print(f"[valid] best epoch {best['epoch']}/{num_epochs} | "
          f"train loss {best['train_loss']:.4f}, valid loss {best['valid_loss']:.4f} | "
          f"train auc {best['train_auc']:.3f}, valid auc {v_auc:.3f} | "
          f"valid aupr {v_aupr:.3f}, f1 {v_f1:.3f}, mcc {v_mcc:.3f}, acc {v_acc:.3f} | "
          f"threshold {best['threshold']:.4f}")
    print(f"[valid] confusion_matrix\n{v_cm}")
    if best['valid_auc_per_disease'] is not None:
        # v_auc above is the pooled number, which is not what selected this
        # epoch. Print the macro one that did, and the per-disease AUCs behind
        # it -- the point of the mode is to see whether any disease was left
        # behind, which the average alone hides.
        per = best['valid_auc_per_disease']
        print(f"[macro] valid auc {best['valid_auc']:.3f} "
              f"(mean over {len(per)} diseases; the pooled {v_auc:.3f} above "
              f"is not what selected this epoch)")
        print("[macro] " + "  ".join(f"{g} {a:.3f}"
                                     for g, a in sorted(per.items())))
        worst = min(per, key=per.get)
        print(f"[macro] weakest disease {worst} {per[worst]:.3f}, "
              f"spread {max(per.values()) - min(per.values()):.3f}")

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
        valid_biom,
        test_biom,
        embedding_birnn,
        plotfile_loss,
        plotfile_auc,
        labels_col="group",
        sample_id_col="sample_id",
        num_steps=600,
        p_drop=0,
        batch_size=64,
        d_model=100,
        n_layers=1,
        n_heads=1,
        numb=1,
        lr=0.001,
        weight_decay=0,
        num_epochs=100,
        loss="BCE_loss",
        alpha=0.6,
        glove_embedding=None,
        pred_out=None,
        d_ff=None,
        head_hidden=None,
        abund_mode='multiply',
        model_seed=11,
        patience=10,
        min_delta=0.001,
        pos_weight=None,
        select_by='loss',
        linear_branch=True,
        lin_weight_decay=1e-3,
        disease_col='disease_name_ab',
        valid_auc='pooled',
        min_group_n=0,
        logit_adjust_tau=1.0):
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
        Path to the training BIOM file (features x samples). Used in full --
        there is no internal split; the caller decides the partition.
    valid_biom : str
        Path to the validation BIOM file, in the same format. This cohort
        drives early stopping, checkpoint selection and the decision
        threshold; it is never trained on.
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
    loss : {'BCE_loss', 'BCEWithLogits', 'FocalLoss'}, optional
        Which criterion to build. Default is ``'BCE_loss'``.
        ``'BCEWithLogits'`` is the same objective computed in its numerically
        stable fused form, and is the only one that accepts `pos_weight`.
    alpha : float, optional
        Positive-class weight, used only by ``'FocalLoss'``. Default is 0.6.
    glove_embedding : str, optional
        Path to a pretrained embedding file used to initialize the OTU
        embedding table; its vector width must equal `d_model`. Default is
        None, which fills the table with fixed random codes instead.
    pred_out : str, optional
        Prefix for the prediction files, written as ``<pred_out>_valid.csv``
        and ``<pred_out>_test.csv``. Default is None, which derives the prefix
        from `embedding_birnn`.
    d_ff : int, optional
        Width of the feed-forward hidden layer inside each encoder layer.
        Default is None, which uses ``4 * d_model``. This is the widest tensor
        in the model, so it is also the first knob to turn down if the encoder
        overfits or runs out of memory.
    head_hidden : int, optional
        Width of the hidden layer in the output head. Default is None, which
        uses ``d_model // 2``. Pass 0 for the single-linear-layer head.
    abund_mode : {'multiply', 'none'}, optional
        How abundance enters the input representation. Default is
        ``'multiply'``, which scales each OTU embedding by its abundance.
        ``'none'`` ignores abundance, leaving the model only the presence or
        absence of each covered taxon; comparing the two measures what
        abundance is contributing under the current normalization.
    model_seed : int, optional
        Seed for weight initialization, dropout, batch shuffling, and the
        random embedding draw. Default is 11. Since the partition is now fixed
        by the input files rather than by a seed, varying `model_seed` alone
        re-runs the same three cohorts with a different model -- which is what
        averaging or ensembling over seeds needs.
    patience : int, optional
        Stop once the validation loss has not improved by at least `min_delta`
        for this many epochs. Default is 10.
    min_delta : float, optional
        Minimum absolute decrease in validation loss that counts as an
        improvement for model selection and early stopping. Default is 0.001.
    pos_weight : float or {'auto'}, optional
        Weight applied to the positive class, valid only with
        ``loss='BCEWithLogits'``. ``'auto'`` uses the ratio of negatives to
        positives in the training part. Default is None, which weights the
        classes equally. Note this shifts where the probabilities sit rather
        than how the samples are ordered, so it moves the thresholded metrics
        far more than it moves the AUC.
    disease_col : str, optional
        Metadata column naming each sample's disease, read only under
        ``valid_auc='macro'`` and ``loss='LogitAdjusted'``. Default is
        ``'disease_name_ab'``.
    valid_auc : {'pooled', 'macro'}, optional
        How the validation AUC that selects the epoch is formed. ``'pooled'``
        (default) is one AUC over every validation sample at once -- unchanged
        from before this parameter existed. ``'macro'`` scores each disease
        (`disease_col`) separately and averages, so a disease with many samples
        no longer dominates the choice of epoch, and the cross-disease pairs a
        pooled AUC counts -- which are most of them once several diseases are
        present -- are dropped. See :func:`macro_auc`. Diseases whose
        validation samples are all one class cannot be scored and are left out,
        reported once at the first epoch. Only the AUC changes: the loss and
        the decision threshold stay pooled, since under a left-out-disease
        split the test disease is absent from validation and a per-disease
        threshold would have nothing to transfer from.
    min_group_n : int, optional
        Under ``valid_auc='macro'``, leave out diseases with fewer than this
        many validation samples. Default is 0, which keeps every disease that
        holds both classes. Equal votes equalise the level but not the
        variance -- a disease with a couple of dozen validation samples swings
        by around 0.13 on sampling noise alone, enough to decide which epoch
        the average prefers -- so the smallest cohorts are worth excluding
        from the *selection* even though the model is still trained on them.
        Read the ``[macro] scored:`` line of a run to pick a value.
    logit_adjust_tau : float, optional
        Strength of the logit adjustment applied when ``loss='LogitAdjusted'``,
        ignored otherwise. Default is 1.0, which removes each disease's
        case/control base rate exactly (Menon et al., ICLR 2021); 0 leaves the
        logits alone and reduces the loss to plain ``BCEWithLogits``. Requires
        `disease_col`. See :func:`logit_adjust_offsets`.

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
    other. The random draw comes from `model_seed`, so a rerun reproduces it.
    """
    set_seed(model_seed)
    devices = try_all_gpus(numb)

    # if glove_embedding is not None:
    #     glove_embedding = f"{glove_embedding}_{d_model}.txt"

    full_train = load_data_imdb(train_biom,
                                metadata,
                                labels_col,
                                sample_id_col,
                                num_steps)
    fid_dict = full_train()
    # Encode the held-out tables with the training vocabulary, so all three
    # share one index space and the model's embedding rows mean the same OTU.
    valid_data = load_data_imdb(valid_biom,
                                metadata,
                                labels_col,
                                sample_id_col,
                                num_steps,
                                fid=fid_dict)
    test_data = load_data_imdb(test_biom,
                               metadata,
                               labels_col,
                               sample_id_col,
                               num_steps,
                               fid=fid_dict)

    train_labels = to_numpy(full_train.labels)

    train_iter = DataLoader(full_train,
                            batch_size=batch_size,
                            shuffle=True)
    valid_iter = DataLoader(valid_data,
                            batch_size=batch_size,
                            shuffle=False)
    test_iter = DataLoader(test_data,
                           batch_size=batch_size,
                           shuffle=False)
    print(f"[split] train {len(full_train)}, valid {len(valid_data)}, "
          f"test {len(test_data)} (all three from separate tables)")

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
                              pad_id=fid_dict['<pad>'],
                              d_ff=d_ff,
                              head_hidden=head_hidden,
                              abund_mode=abund_mode,
                              linear_branch=linear_branch)
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

    # Counted after the embedding is frozen, so this is the number of weights
    # the optimizer actually moves.
    n_trainable = sum(p.numel() for p in net.parameters() if p.requires_grad)
    print(f"[model] d_model {d_model}, layers {n_layers}, heads {n_heads}, "
          f"d_ff {net.layers[0].d_ff}, dropout {p_drop} | "
          f"abundance {abund_mode}, model_seed {model_seed} | "
          f"{n_trainable} trainable parameters")

    if linear_branch:
        # The branch holds one free coefficient per OTU -- more parameters than
        # the rest of the encoder combined -- so it needs its own, much
        # stronger decay. Sharing `weight_decay` with the attention path would
        # either leave the branch unregularized or over-regularize the encoder.
        lin_p = [p for n, p in net.named_parameters()
                 if n.startswith('lin_') and p.requires_grad]
        oth_p = [p for n, p in net.named_parameters()
                 if not n.startswith('lin_') and p.requires_grad]
        print(f"[model] linear residual branch on | {len(lin_p)} branch "
              f"tensors, lin_weight_decay {lin_weight_decay}")
        trainer = torch.optim.Adam(
            [{'params': oth_p, 'weight_decay': weight_decay},
             {'params': lin_p, 'weight_decay': lin_weight_decay}], lr=lr)
    else:
        trainer = torch.optim.Adam(
            net.parameters(), lr=lr, weight_decay=weight_decay)
    if pos_weight is not None and loss != "BCEWithLogits":
        raise ValueError(f"pos_weight applies only to loss='BCEWithLogits', "
                         f"not {loss!r}"
                         + (" -- LogitAdjusted already corrects a class prior, "
                            "per disease, so stacking pos_weight on top would "
                            "correct it twice"
                            if loss == "LogitAdjusted" else ""))
    if loss == "FocalLoss":
        loss = FocalLoss(alpha=alpha, devices=devices)
    elif loss == "BCEWithLogits":
        pw = None
        if pos_weight is not None:
            if pos_weight == 'auto':
                # 'auto' describes the class balance of what the loader
                # actually yields, which is the whole training part.
                n_pos = int((train_labels == 1).sum())
                n_neg = int((train_labels == 0).sum())
                if n_pos == 0:
                    raise ValueError("pos_weight='auto' needs at least one "
                                     "positive sample in the training part")
                pw_value = n_neg / n_pos
            else:
                pw_value = float(pos_weight)
            print(f"[loss ] BCEWithLogits, pos_weight {pw_value:.3f}")
            pw = torch.tensor(pw_value, device=devices[0])
        loss = nn.BCEWithLogitsLoss(pos_weight=pw)
    elif loss == "LogitAdjusted":
        # The offsets ride on the training dataset, so they follow each sample
        # through the loader's shuffling.
        la_strata = read_strata(metadata, sample_id_col, disease_col,
                               full_train.sample_ids)
        offsets, _ = logit_adjust_offsets(la_strata, train_labels,
                                          tau=logit_adjust_tau)
        full_train.set_logit_offset(offsets)
        loss = LogitAdjustedLoss(pos_weight=None, reduction='mean')
    elif loss == "BCE_loss" or loss is None:
        loss = nn.BCELoss(weight=None, reduction='mean')
    else:
        raise ValueError(f"unsupported loss {loss!r}; expected 'BCE_loss', "
                         f"'BCEWithLogits', 'FocalLoss', or 'LogitAdjusted'")

    print(f"[train] lr {lr}, epochs {num_epochs}, patience {patience}, "
          f"min_delta {min_delta}, select_by {select_by}")

    # valid_iter does not shuffle, so predict_iter returns probabilities in the
    # dataset's own row order and the strata line up index for index.
    if valid_auc == 'pooled':
        valid_strata = None
    elif valid_auc == 'macro':
        valid_strata = read_strata(metadata, sample_id_col, disease_col,
                                   valid_data.sample_ids)
    else:
        raise ValueError(
            f"unsupported valid_auc {valid_auc!r}; expected 'pooled' or 'macro'")

    valid_record = train_cls(net, train_iter, valid_iter, loss, trainer,
                             num_epochs, devices, embedding_birnn,
                             plotfile_loss, plotfile_auc,
                             patience=patience, min_delta=min_delta,
                             select_by=select_by,
                             valid_strata=valid_strata,
                             min_group_n=min_group_n)

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
                     valid_data.sample_ids,
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