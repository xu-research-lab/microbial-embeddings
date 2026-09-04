#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Per-taxon attributions for the LOSO attention ensemble.

This replaces ``run_shap_python_crc.py``, which explained the previous model
(``membed.Attention_embedding.TransformerEncoder``: a CLS token, one checkpoint
per study, and its own hand-written tokenizer). Three things about the current
model make that script wrong rather than merely dated:

* :class:`~membed.otu_attention.OtuAttentionEncoder` has no CLS token, masks
  ``'<unk>'`` alongside ``'<pad>'``, and returns a raw logit rather than a
  probability.
* A fold is not one model. Under ``--inner-split loso`` a fold trains one
  member per training cohort, and the number this run reports is
  ``--report-combiner prob`` -- the mean of the members' probabilities. What
  gets explained here is therefore that mean, not any single member.
* The model carries a linear residual branch, so its logit is
  ``encoder(pooled) + sum_i lin_w[otu_i] * abundance_i + lin_b``. The second
  term is exactly attributable per taxon, which this script uses as a hard
  check on the gradient path rather than as a separate result.

What is explained
-----------------
For each fold, the samples of the held-out part -- the ones the fold actually
predicts -- against a background of the per-taxon median over the fold's
training and test tables together. The feature space is the union of the two
tables' observation axes, so a taxon the test part lacks but the background has
still gets an attribution, which is the point of a held-out split.

What a fold is comes from the task, and nothing here hard-codes it: ``disease``
holds out one study of one disease (``{disease}_{study}``) and ``lodo`` holds
out a whole disease (``{left_out_disease}``), so the run-TSV columns are read
off each task's own fold template. Two consequences worth knowing under lodo:

* Its training table is every other disease at once, so the background median
  is computed blocked (:func:`union_median`) rather than by stacking the two
  tables dense, which for the pooled table would be tens of gigabytes.
* Each disease is a single fold, so the cross-cohort consistency table the
  ``disease`` task produces has nothing to be consistent across. It is written
  as a plain ranking there, and the columns that would read as agreement are
  left out.

Three readouts
--------------
``attn``
    What the attention layer looks at. ``attn_weights`` is already a
    distribution over the taxa a sample carries -- the softmax runs over the
    key axis with masked keys driven to -1e9 -- so the readout is the mean of
    those rows over the unmasked queries, one forward pass per sample per
    member. No integration, no backward pass, two to three orders of magnitude
    cheaper than the two methods below, and it is the direct answer to "which
    microbes does the model attend to". Runs on every sample of every fold.

    Read it as **enrichment**, ``attention * n_valid``, which is what the files
    hold: 1.0 is "attended as much as an average taxon of this sample". And
    read the **case-minus-control differential** rather than the level. The
    attention logits are ``q.k/sqrt(d_k)`` with both ``q`` and ``k`` linear in
    ``embedding * abundance``, so the scores scale with abundance squared and
    the raw ranking is substantially a ranking of abundance -- the run reports
    ``attn_abundance_spearman`` per fold so that is a number rather than a
    surprise. The shared part cancels out of a case-minus-control difference.

    Attention is not attribution: it says what the encoder reads, not what the
    output responds to. The two below say the latter, and
    ``attn_vs_grad_{task}.csv`` reports how far apart they land.

``grad``
    Integrated gradients on the abundance vector, with the sample's own
    tokenization held fixed and the abundance moved from the background's value
    to the sample's. Cheap enough to run on every fold and every held-out
    sample: 2 * ``--ig-steps`` forward/backward passes per sample, against
    ``2 * n_features + 1`` forward passes for the permutation explainer.

    The approximation is the fixed tokenization: a real background sample would
    select a different set of top-``num_steps`` taxa, and this holds the
    sample's own selection while sliding the abundances. Positions are masked
    exactly as in training, so ``'<pad>'`` and ``'<unk>'`` contribute nothing.

``perm``
    ``shap.Explainer(f, background)`` -- the permutation explainer, the same
    call the old script made, on the same background and the same feature
    space. Run on a couple of folds only, to measure what the fixed-tokenization
    approximation costs: the two methods are compared on the per-taxon ordering
    of ``mean|attribution|``.

Both explain the same difference, ``f(x) - f(background)``, so their values are
on one scale and worth correlating.

Reproducing the runs being explained
------------------------------------
Neither ``results_with_SNEs`` run passed ``--keep-ckpt``, so both deleted their
weights. Retrain what you want to explain first, keeping the checkpoints, under
a run name of its own so the original results are left alone. The two runs
differ in more than the task -- lodo was trained with ``--no-linear-branch``,
``--inner-split disease_loso`` and 32 members -- so each retrain has to repeat
its own run's flags exactly:

    # disease: only the CRC and IBD folds, so cut the fold list first
    python - <<'PY'
    import pandas as pd
    t = pd.read_csv("run_leave_one_study_out_each_diease_list.tsv", sep="\\t")
    t[t.disease.isin(["CRC", "IBD"])].to_csv(
        "run_leave_one_study_out_explain.tsv", sep="\\t", index=False)
    PY

    python run_attention_biom_with_SNEs.py --tasks disease \\
        --gpus 0 1 2 3 4 5 6 7 --inner-split loso --linear-branch \\
        --report-combiner prob \\
        --run-tsv run_leave_one_study_out_explain.tsv \\
        --run-name results_with_SNEs_ckpt --keep-ckpt

    # lodo: every left-out disease
    python run_attention_biom_with_SNEs.py --tasks lodo \\
        --gpus 0 1 2 3 4 5 6 7 --inner-split disease_loso --valid-auc macro \\
        --n-estimators 32 --report-combiner prob --no-linear-branch \\
        --run-name results_with_SNEs_ckpt --keep-ckpt

Cutting the fold list is safe: a job's ``model_seed`` is ``--ensemble-seed0``
plus its index among its own fold's members, its ``member_seed`` is
``--member-seed`` plus the same index, and the members are enumerated from that
fold's own training table, so the surviving jobs are bit-for-bit the ones the
full run produced. Check that by comparing each new results CSV against the
original one before trusting anything below.

A note on reading the lodo attributions. With ``--no-linear-branch`` the whole
model is the pooled encoder, and the encoder's first LayerNorm discards the
abundance magnitude by construction -- see the note in
``OtuAttentionEncoder``'s docstring, which is deliberate. So the lodo numbers
say *which taxa are in the sample's top-num_steps*, far more than *how abundant
they are*. That is a property of the model rather than of this method, but it
belongs next to the numbers.

Usage
-----
    python explain_attention.py --task disease \\
        --run-name results_with_SNEs_ckpt \\
        --run-tsv run_leave_one_study_out_explain.tsv \\
        --diseases CRC IBD --gpus 0 1 2 3 4 5 6 7

    python explain_attention.py --task lodo \\
        --run-name results_with_SNEs_ckpt --gpus 0 1 2 3 4 5 6 7 \\
        --perm-per-disease 0

The permutation explainer is the one thing worth running on its own: it is a
single fold on a single GPU while everything else is spread over all of them,
so leaving it in the main run means seven idle cards waiting for one.

    python explain_attention.py --task disease \\
        --run-name results_with_SNEs_ckpt \\
        --run-tsv run_leave_one_study_out_explain.tsv --diseases CRC IBD \\
        --gpus 0 --perm-per-disease 1 --perm-max-folds 1

Is --max-samples 300 enough? It is an empirical question about this data, so
measure it rather than trusting the default: run one mid-sized fold twice,
once capped and once not, and compare that fold's rows of the two
``shap_summary_{task}.csv`` files on ``mean_abs`` -- Spearman and the size of
the top-30 intersection. Below about 0.9 and 25/30, raise the cap.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import time
import traceback
from contextlib import contextmanager
from glob import glob
from string import Formatter

import biom
import numpy as np
import pandas as pd

# The task table, the column names and the training CFG all come from the run
# script, so a fold path or a hyper-parameter can only be defined in one place.
from run_attention_biom_with_SNEs import CFG, LABELS_COL, SAMPLE_ID_COL, TASKS

# torch and membed are imported inside the workers. At module level they would
# load before the pool initializer sets CUDA_VISIBLE_DEVICES, which breaks the
# per-process GPU binding -- the same reason the run script defers them.


# =====================================================================
# 1. Fold and member discovery (parent process, never touches CUDA)
# =====================================================================
def fold_fields(task):
    """The run-TSV columns a task's fold name is built from.

    ``disease`` names its folds ``{disease}_{study}`` and ``lodo`` names them
    ``{left_out_disease}``, so nothing can read a fixed pair of column names.
    Both, though, put the disease first, which is what ``--diseases`` filters
    on and what the summary groups by.
    """
    return [f for _, f, _, _ in Formatter().parse(TASKS[task]["fold"]) if f]


def build_jobs(args):
    """One job per fold, carrying its tables and its members' checkpoints."""
    t = TASKS[args.task]
    run_tsv = args.run_tsv or t["run_tsv"]
    run = pd.read_csv(run_tsv, sep="\t")
    artifact = f"{t['artifact']}{args.run_name}"
    fields = fold_fields(args.task)

    jobs, skipped = [], []
    for _, r in run.iterrows():
        row = r.to_dict()
        fold = t["fold"].format(**row)
        disease = str(row[fields[0]])
        # The cohort under test, for the tasks that hold one out. Empty under
        # lodo, whose fold is a whole disease rather than one of its studies.
        study = str(row[fields[1]]) if len(fields) > 1 else ""
        if args.diseases and disease not in args.diseases:
            continue

        fold_dir = os.path.join(artifact, fold)
        members = sorted(
            d for d in glob(os.path.join(fold_dir, "members", "*"))
            if os.path.isfile(os.path.join(d, "model.pth")))
        if not members:
            skipped.append(fold)
            continue

        jobs.append(dict(
            task=args.task, fold=fold, disease=disease, study=study,
            train=t["train"].format(**row), test=t["test"].format(**row),
            meta=t["meta"].format(**row),
            members=members, out_dir=args.out_dir,
            num_steps=args.num_steps, n_heads=args.n_heads,
            abund_mode=args.abund_mode,
            ig_steps=args.ig_steps, ig_cal_steps=args.ig_cal_steps,
            ig_batch=args.ig_batch,
            ig_power=args.ig_power, ig_tol=args.ig_tol, tf32=args.tf32,
            max_samples=args.max_samples, sample_seed=args.sample_seed,
            target=args.target, batch_size=args.batch_size,
            per_member=not args.no_per_member,
            dump_pooled=args.dump_pooled, skip_ig=args.skip_ig,
            self_check=not args.no_self_check,
            check_tol=args.check_tol, check_background=args.check_background,
            median_block=args.median_block,
            perm=False, perm_max_evals=args.perm_max_evals,
            perm_samples=args.perm_samples, perm_chunk=args.perm_chunk,
            perm_seed=args.perm_seed))

    if not jobs and not skipped:
        raise SystemExit(
            f"no fold of {run_tsv!r} is left to explain. --diseases "
            f"{args.diseases} filters on the {fields[0]!r} column, which holds "
            f"{sorted(set(run[fields[0]].astype(str)))}")
    if skipped:
        raise SystemExit(
            f"{len(skipped)} fold(s) under {artifact!r} have no member with a "
            f"model.pth: {', '.join(skipped[:8])}"
            + (" ..." if len(skipped) > 8 else "")
            + "\nThe run being explained deleted its checkpoints. Retrain "
              "those folds with --keep-ckpt (see this script's docstring) "
              "before explaining them.")
    return jobs


def choose_perm_folds(jobs, args):
    """Mark the folds that also get the permutation explainer.

    Explicit names win. Otherwise the smallest held-out cohort of each disease
    is taken, since the permutation explainer costs ``2 * n_features + 1``
    forward passes *per sample* and the check does not need a big cohort to be
    informative. The size is read from a member's ``pred_test.csv`` rather than
    by loading the BIOM table, which the parent has no other reason to touch.

    ``--perm-max-folds`` then caps the total. Per-disease selection alone is
    not a cap under lodo, where each disease *is* one fold and "one per
    disease" therefore means every fold in the run.
    """
    if args.perm_folds:
        wanted = set(args.perm_folds)
        for j in jobs:
            j["perm"] = j["fold"] in wanted
        missing = wanted - {j["fold"] for j in jobs}
        if missing:
            raise SystemExit(f"--perm-folds names folds that are not in this "
                             f"run: {', '.join(sorted(missing))}")
        return

    if args.perm_per_disease <= 0:
        return

    def n_test(job):
        p = os.path.join(job["members"][0], "pred_test.csv")
        try:
            with open(p) as f:
                return sum(1 for _ in f) - 1
        except OSError:
            return np.inf

    by_disease = {}
    for j in jobs:
        by_disease.setdefault(j["disease"], []).append(j)
    picked = []
    for _, group in sorted(by_disease.items()):
        picked += sorted(group, key=n_test)[:args.perm_per_disease]
    for j in sorted(picked, key=n_test)[:args.perm_max_folds]:
        j["perm"] = True


# =====================================================================
# 2. Loading a fold: vocabulary, tables, members
# =====================================================================
def load_fold(job):
    """Tokenizer, dense feature space, background row, and the member models.

    The vocabulary is rebuilt from the fold's own ``train_loo.biom`` rather
    than from each member's training part. Those two are the same index:
    :class:`~membed.otu_attention.Fid` is built from a table's observation axis
    alone, and ``member_split`` filters the sample axis only, so every member of
    a fold shares one vocabulary and it is the fold's. The checkpoint's
    embedding table is the right width by construction, and that is asserted
    below rather than assumed.
    """
    import torch
    from membed.otu_attention import load_data_imdb

    num_steps = job["num_steps"]
    train_ds = load_data_imdb(job["train"], job["meta"], LABELS_COL,
                              SAMPLE_ID_COL, num_steps)
    fid = train_ds()
    test_ds = load_data_imdb(job["test"], job["meta"], LABELS_COL,
                             SAMPLE_ID_COL, num_steps, fid=fid)

    train_cols = np.asarray(biom.load_table(job["train"]).ids(axis="observation"))
    # The feature space is the held-out table's own observation axis, in its
    # own order, and nothing may be added to it, dropped from it or reordered
    # within it.
    #
    # This is not a stylistic choice. `truncate_pad` selects a sample's taxa
    # with `np.argsort(otu[i, ])[::-1]`, and after the rank-over-max
    # normalisation a sample's values are massively tied -- a fixture row had 5
    # distinct values across 20 non-zero taxa. `np.argsort` defaults to
    # quicksort, which is not stable, so which of a set of tied taxa lands in
    # the top-`num_steps` depends on the length and layout of the array being
    # sorted. Widen the space with all-zero columns for taxa only the training
    # table declares, and the tie-breaking moves: measured on the lodo fixture,
    # the rebuilt predictions drifted from the run's own by 3.9e-02, which the
    # forward check refuses.
    #
    # What this costs is the taxa the held-out table does not declare. They are
    # zero in every explained sample, so the gradient method could never
    # attribute to them anyway -- it only visits positions in a sample's own
    # tokenization -- and the permutation explainer loses the "the background
    # has this and the sample does not" direction for them alone. When the two
    # tables are sample-wise splits of one table, which is how this pipeline
    # builds them, the two axes are identical and nothing is lost at all.
    columns = np.asarray(biom.load_table(job["test"]).ids(axis="observation"))
    col_of = {c: i for i, c in enumerate(columns.tolist())}
    if not np.array_equal(train_cols, columns):
        print(f"[warn] {job['fold']}: the training and held-out tables declare "
              f"different observation axes ({len(train_cols)} vs "
              f"{len(columns)}); the {len(set(train_cols.tolist()) - set(columns.tolist()))} "
              f"taxon(s) only the training table has are outside the explained "
              f"feature space, and contribute to the background only where the "
              f"held-out table also declares them")

    background = union_median(train_ds.otu, train_cols, test_ds.otu, columns,
                              col_of, len(columns), job)
    if job["check_background"]:
        naive = naive_union_median(train_ds.otu, train_cols, test_ds.otu,
                                   columns, col_of, len(columns))
        if not np.array_equal(background, naive):
            bad = int((background != naive).sum())
            raise ValueError(
                f"{job['fold']}: the blocked background median differs from "
                f"the definition in {bad}/{len(columns)} column(s), worst by "
                f"{np.abs(background - naive).max():.3e}. The zero shortcut or "
                f"the blocking is wrong; every attribution is measured against "
                f"this row")
    # The training table's dense matrix is the largest object a fold loads --
    # under lodo it is every disease at once -- and the background above was
    # the last thing that needed it. What the dataset is kept for afterwards is
    # its vocabulary and its `truncate_pad`, and that method reads nothing but
    # `self.fid`.
    train_ds.otu = None

    # The held-out table's matrix as the library built it, not a copy placed
    # into some other space. See the note above.
    X_test = test_ds.otu

    # Token id -> column, precomputed. The alternative is a dict lookup per
    # position per sample, and this is called once per integration step. Only
    # IDs the vocabulary actually covers are entered, so '<pad>', '<unk>' and
    # every test-only taxon that collapses onto '<unk>' stay at -1 and are
    # dropped from the scatter -- which is right, since they are masked out of
    # the model too.
    col_of_token = np.full(len(fid), -1, dtype=np.int64)
    for c in columns.tolist():
        if c in fid.token_to_idx:
            col_of_token[fid.token_to_idx[c]] = col_of[c]

    nets = [_build_member(d, fid, job, torch) for d in job["members"]]
    return dict(fid=fid, train_ds=train_ds, test_ds=test_ds, columns=columns,
                col_of=col_of, col_of_token=col_of_token, X=X_test,
                background=background, nets=nets,
                labels=test_ds.labels.numpy(),
                sample_ids=np.asarray(test_ds.sample_ids))


def union_median(train_otu, train_cols, test_otu, test_cols, col_of, width,
                 job):
    """Per-taxon median over both tables, without ever stacking them dense.

    Over the held-out table's feature space: a taxon only the training table
    declares is not a column here at all, and a taxon only the held-out table
    declares takes zeros from the training side, which is what "not observed"
    means. `test_cols` is that space and `col_of` indexes it.

    The straightforward version -- widen both matrices, ``vstack``,
    ``np.median(axis=0)`` -- is fine for a single-disease table and is tens of
    gigabytes for lodo's pooled one, whose training half is every other disease
    at once. Two exact reductions replace it:

    * **The zeros settle most columns.** Abundances are non-negative, so a
      column whose non-zero count `k` satisfies ``2 * k < n - 1`` has zeros at
      both middle order statistics and a median of exactly 0. That is most of a
      sparse OTU table, and counting non-zeros costs one pass with no widening.
    * **The rest are done in column blocks**, so the peak allocation is
      ``n_rows x --median-block`` rather than ``n_rows x n_features``.

    Neither is an approximation. ``--check-background`` asserts the result
    against the naive computation element for element, which is worth doing on
    a small fold whenever this function is touched.
    """
    n_rows = train_otu.shape[0] + test_otu.shape[0]
    # Training columns outside the held-out table's axis are not columns of
    # this space and are dropped; the ones inside it keep their own position in
    # `train_otu`, which the blocked pass reads them from.
    tr_keep = np.array([i for i, c in enumerate(train_cols.tolist())
                        if c in col_of], dtype=np.int64)
    tr_idx = np.array([col_of[train_cols[i]] for i in tr_keep], dtype=np.int64)
    te_idx = np.array([col_of[c] for c in test_cols.tolist()], dtype=np.int64)

    nnz = np.zeros(width, dtype=np.int64)
    if len(tr_keep):
        np.add.at(nnz, tr_idx, np.count_nonzero(train_otu[:, tr_keep], axis=0))
    np.add.at(nnz, te_idx, np.count_nonzero(test_otu, axis=0))
    # The boundary is left to the blocked path rather than reasoned about: at
    # 2k == n - 1 or 2k == n one of the two middle values is the smallest
    # non-zero, and the median is not 0.
    todo = np.flatnonzero(2 * nnz >= n_rows - 1)

    background = np.zeros(width, dtype=np.float64)
    if len(todo) == 0:
        return background

    # Column of this space -> column of each source matrix, -1 where absent.
    src_train = np.full(width, -1, dtype=np.int64)
    src_train[tr_idx] = tr_keep
    src_test = np.full(width, -1, dtype=np.int64)
    src_test[te_idx] = np.arange(len(te_idx))

    block = job["median_block"]
    buf = np.zeros((n_rows, min(block, len(todo))), dtype=np.float64)
    n_tr = train_otu.shape[0]
    for i in range(0, len(todo), block):
        cols = todo[i:i + block]
        b = buf[:, :len(cols)]
        b[:] = 0.0
        st, se = src_train[cols], src_test[cols]
        b[:n_tr, st >= 0] = train_otu[:, st[st >= 0]]
        b[n_tr:, se >= 0] = test_otu[:, se[se >= 0]]
        background[cols] = np.median(b, axis=0)
    return background


def naive_union_median(train_otu, train_cols, test_otu, test_cols, col_of,
                       width):
    """The definition `union_median` is optimised away from. Testing only."""
    return np.median(np.vstack([
        _widen(train_otu, train_cols, col_of, width),
        _widen(test_otu, test_cols, col_of, width)]), axis=0)


def _widen(otu, cols, col_of, width):
    """A table's matrix in the explained feature space, zero-filled.

    Columns the space does not contain are dropped. Only the background
    computation and its test use this; the matrix that reaches the tokenizer is
    never widened, since that would change which tied taxa it selects.
    """
    src = [i for i, c in enumerate(cols.tolist()) if c in col_of]
    out = np.zeros((otu.shape[0], width), dtype=np.float64)
    out[:, [col_of[cols[i]] for i in src]] = otu[:, src]
    return out


def _build_member(member_dir, fid, job, torch):
    """Rebuild one member from its checkpoint.

    Every shape is read off the ``state_dict``, so a member trained under a
    different ``--set`` than this script's defaults still loads. ``n_heads`` is
    the exception -- the projections are ``d_model x d_model`` whatever it was
    -- so it comes from CFG and can be overridden with ``--n-heads``.
    """
    from membed.otu_attention import OtuAttentionEncoder

    path = os.path.join(member_dir, "model.pth")
    sd = torch.load(path, map_location="cpu")

    otu_size, d_model = sd["embedding.weight"].shape
    if otu_size != len(fid):
        raise ValueError(
            f"{path}: the checkpoint's embedding table has {otu_size} rows but "
            f"the vocabulary rebuilt from {job['train']} has {len(fid)}. The "
            f"member was trained against a different table than this fold's "
            f"train_loo.biom, so its embedding rows do not mean the taxa this "
            f"script would feed it")
    n_layers = 1 + max(int(k.split(".")[1]) for k in sd if k.startswith("layers."))
    # 0 rather than None when the layer has no feed-forward sub-layer at all:
    # None means "4 * d_model" to EncoderLayer, which is a different model.
    d_ff = (sd["layers.0.ffn.0.weight"].shape[0]
            if "layers.0.ffn.0.weight" in sd else 0)
    head_hidden = 0 if sd["mlp.0.weight"].shape[0] == 1 else sd["mlp.0.weight"].shape[0]
    linear_branch = "lin_w" in sd

    net = OtuAttentionEncoder(
        otu_size=otu_size, d_model=d_model, n_layers=n_layers,
        n_heads=job["n_heads"], p_drop=0.0, pad_id=fid["<pad>"],
        d_ff=d_ff, head_hidden=head_hidden, abund_mode=job["abund_mode"],
        linear_branch=linear_branch)
    net.load_state_dict(sd)
    net.eval()
    return net


# =====================================================================
# 3. The prediction the explanation is of
# =====================================================================
def tokenize(fold, X):
    """Dense matrix -> the exact tensors the training pipeline would build.

    ``truncate_pad`` is reused rather than reimplemented. It is a bound method
    but touches nothing on the instance except ``self.fid``, so calling it with
    a perturbed matrix puts that matrix through the very code path training and
    inference use. The old script kept a second copy of this logic and the copy
    was wrong for the current model: it prepended a CLS column and left
    ``'<unk>'`` unmasked.
    """
    return fold["train_ds"].truncate_pad(np.asarray(X, dtype=np.float64),
                                         fold["columns"],
                                         fold["train_ds"].features.shape[1])


def member_probs(fold, feats, abund, mask, device, batch_size=64):
    """Every member's probability for every row: ``(n_members, n_samples)``."""
    import torch
    out = []
    with torch.no_grad():
        for net in fold["nets"]:
            net.to(device)
            ps = []
            for i in range(0, feats.shape[0], batch_size):
                sl = slice(i, i + batch_size)
                logit, _ = net(feats[sl].to(device), abund[sl].to(device),
                               mask[sl].to(device))
                ps.append(torch.sigmoid(logit.squeeze(1)).cpu().numpy())
            out.append(np.concatenate(ps))
    return np.vstack(out)


def predict(fold, X, device, batch_size=64):
    """``f(X)``: the mean of the members' probabilities, as the run reports it."""
    feats, abund, mask = tokenize(fold, X)
    return member_probs(fold, feats, abund, mask, device, batch_size).mean(axis=0)


# =====================================================================
# 4. Self-checks
# =====================================================================
def check_forward(fold, job, device):
    """Does this script's pipeline reproduce the run's own predictions?

    The single most important check here. Each member wrote a ``pred_test.csv``
    during the run; feeding the held-out table through ``tokenize`` and the
    reloaded checkpoint has to return those same probabilities. If it does,
    the explanation is of the model that was scored, and not of a lookalike
    assembled from the wrong vocabulary, the wrong feature order, or the wrong
    head width.
    """
    feats, abund, mask = tokenize(fold, fold["X"])
    P = member_probs(fold, feats, abund, mask, device, job["batch_size"])

    worst, checked = 0.0, 0
    for k, d in enumerate(job["members"]):
        p = os.path.join(d, "pred_test.csv")
        if not os.path.isfile(p):
            continue
        ref = pd.read_csv(p)
        if not np.array_equal(ref["sample_id"].astype(str).to_numpy(),
                              fold["sample_ids"].astype(str)):
            raise ValueError(
                f"{p}: sample order differs from the test table's, so the "
                f"comparison would be between different samples")
        worst = max(worst, float(np.abs(P[k] - ref["prob"].to_numpy()).max()))
        checked += 1

    if not checked:
        raise ValueError(f"{job['fold']}: no member has a pred_test.csv to "
                         f"check the forward pass against")
    if worst > job["check_tol"]:
        raise ValueError(
            f"{job['fold']}: rebuilt predictions differ from the run's "
            f"pred_test.csv by up to {worst:.2e} (tolerance "
            f"{job['check_tol']:.0e}) over {checked} member(s). The "
            f"explanation would be of a different model than the one scored; "
            f"fix this before reading any attribution")
    return dict(max_abs_prob_diff=worst, members_checked=checked)


def check_linear_branch(fold, job, device, n_steps):
    """The linear branch's attribution has a closed form; the code must hit it.

    Integrated gradients are exact for a linear term, and the branch is exactly
    ``sum_i lin_w[token_i] * abundance_i * mask_i``. Running the first member
    twice -- branch on, branch off -- and subtracting isolates what the code
    attributed to it, which must equal ``lin_w[token_i] * (x_i - base_i)``
    position by position. Nothing else tests the gradient path this sharply:
    a wrong baseline, a misaligned scatter or a dropped mask all show up here.

    Done in logit space, where the branch is additive. In probability space it
    is not separable and the identity does not hold.
    """
    net = fold["nets"][0]
    if not getattr(net, "linear_branch", False):
        # Under --no-linear-branch, as the lodo run was trained, there is no
        # branch to check. The caller records that this check did not run
        # rather than letting its absence read as a pass.
        return None

    feats, abund, mask, base = _ig_inputs(fold, job)
    kw = dict(device=device, n_steps=n_steps, batch_rows=job["ig_batch"],
              target="logit", power=job["ig_power"])
    with_branch, _, _, _ = integrated_gradients([net], feats, abund, mask,
                                                base, **kw)
    net.linear_branch = False
    try:
        without, _, _, _ = integrated_gradients([net], feats, abund, mask,
                                                base, **kw)
    finally:
        net.linear_branch = True

    lin_w = net.lin_w.detach().cpu().numpy()
    exact = lin_w[feats.numpy()] * (abund.numpy() - base) * mask.numpy()
    got = with_branch - without
    scale = max(float(np.abs(exact).max()), 1e-12)
    err = float(np.abs(got - exact).max()) / scale
    if err > 1e-3:
        raise ValueError(
            f"{job['fold']}: the integrated gradients attribute "
            f"{err:.2e} (relative) more or less to the linear branch than its "
            f"closed form gives. The gradient path is wrong -- check the "
            f"baseline, the mask and the position->taxon scatter")
    return dict(linear_branch_rel_err=err)


# =====================================================================
# 5. Integrated gradients
# =====================================================================
def _ig_inputs(fold, job):
    """Tokenized held-out samples plus the background abundance per position.

    The baseline is the background's value *for the taxon that sits at that
    position*, so the path moves each taxon's abundance from typical to
    observed while the set of taxa stays the sample's own. Positions holding
    ``'<pad>'`` or ``'<unk>'`` get 0: they are masked out of attention, of the
    pooling and of the linear branch, so nothing they are given can reach the
    output.
    """
    feats, abund, mask = tokenize(fold, fold["X"])
    cols = fold["col_of_token"][feats.numpy()]      # (n_samples, num_steps)

    base = np.zeros(feats.shape, dtype=np.float64)
    hit = cols >= 0
    base[hit] = fold["background"][cols[hit]]
    base *= mask.numpy()
    return feats, abund.double(), mask, base


def integrated_gradients(nets, feats, abund, mask, base, device, n_steps,
                         batch_rows, target="prob", power=3.0):
    """Attribution per position, plus ``f`` at the sample and at the baseline.

    The path is the straight line from `base` to `abund` in abundance space
    with the tokenization held fixed. Rows of different samples are batched
    together freely: each row's output depends on that row alone, so one
    backward pass over their weighted sum yields every row's own gradient.

    The path points are **not** spaced evenly. Under ``power=p`` the integral
    is substituted as ``alpha = u**p``, sampling `u` at the midpoints of
    `n_steps` equal intervals and weighting each point by the Jacobian
    ``p * u**(p-1)``. This is not a refinement; without it the method does not
    work on this model at any affordable step count.

    The reason is in the architecture. Most positions of the baseline are 0 --
    the background median of a sparse taxon is 0 -- so near ``alpha = 0`` the
    whole community shrinks towards the zero vector, where the encoder's
    ``LayerNorm(..., eps=1e-6)`` is dominated by its epsilon and ``f`` swings
    through most of its range inside a boundary layer narrower than
    ``alpha = 0.01``. Measured on a fixture: uniform spacing leaves a
    completeness error of 0.48 that barely moves between 16 and 1024 steps,
    while ``power=3`` at 256 steps closes it to 0.0000. An error that refuses
    to shrink with step count looks like a broken gradient and is not one, so
    it is worth knowing which of the two you are looking at.

    Every member is attributed separately and the ensemble is their mean. That
    is not an extra pass: ``f`` is the mean of the members, integrated
    gradients are linear in ``f``, so the mean of the members' attributions is
    the ensemble's exactly -- and backward over one summed output already
    differentiates through all of their graphs, so one backward per member
    costs what a single backward over the sum did. Per-member attributions,
    which under ``--n-estimators 32`` used to mean a 33-fold bill, are
    therefore free.

    Returns ``(attr, f_x, f_base, per_member)``. `attr` has shape
    ``(n_samples, num_steps)`` and `per_member` ``(n_members, n_samples,
    num_steps)``. ``attr.sum(1) == f_x - f_base`` is the completeness identity,
    which the caller checks.
    """
    import torch

    n, steps = feats.shape
    m = len(nets)
    x = abund.numpy().astype(np.float64)
    dx = x - base
    u = (np.arange(n_steps) + 0.5) / n_steps
    alphas = u ** power
    jac = power * u ** (power - 1.0)

    rows = np.repeat(np.arange(n), n_steps)
    al = np.tile(alphas, n)
    wt = np.tile(jac, n)
    grad_sum = np.zeros((m, n, steps), dtype=np.float64)

    for i in range(0, len(rows), batch_rows):
        js = rows[i:i + batch_rows]
        jt = torch.as_tensor(js, dtype=torch.long)
        a_np = base[js] + al[i:i + batch_rows, None] * dx[js]
        ft = feats[jt].to(device)
        mk = mask[jt].to(device)
        w = torch.tensor(wt[i:i + batch_rows], dtype=torch.float32,
                         device=device)

        for k, net in enumerate(nets):
            net.to(device)
            a = torch.tensor(a_np, dtype=torch.float32, device=device,
                             requires_grad=True)
            logit, _ = net(ft, a, mk)
            v = logit.squeeze(1)
            if target == "prob":
                v = torch.sigmoid(v)
            (v * w).sum().backward()
            # np.add.at, not `+=`: consecutive rows of a batch are path points
            # of the same sample, and fancy-indexed assignment would keep only
            # the last of them.
            np.add.at(grad_sum[k], js,
                      a.grad.detach().cpu().numpy().astype(np.float64))

    per_member = dx[None, :, :] * grad_sum / n_steps
    attr = per_member.mean(axis=0)

    def _f(vals):
        out = []
        with torch.no_grad():
            for i in range(0, n, 64):
                sl = slice(i, i + 64)
                a = torch.tensor(vals[sl], dtype=torch.float32, device=device)
                acc = None
                for net in nets:
                    logit, _ = net(feats[sl].to(device), a, mask[sl].to(device))
                    v = logit.squeeze(1)
                    if target == "prob":
                        v = torch.sigmoid(v)
                    acc = v if acc is None else acc + v
                out.append((acc / len(nets)).cpu().numpy())
        return np.concatenate(out)

    return attr, _f(x), _f(base), per_member


#: Step counts the calibration walks, cheapest first.
IG_STEP_LADDER = (64, 128, 256, 512, 1024)


@contextmanager
def tf32(enabled):
    """Allow TF32 matmuls inside this block, on Ampere and later.

    Worth two to three times on the integration, which is nothing but small
    matmuls, at a relative error around 1e-3. That error is far below the
    completeness tolerance and the check stays self-consistent under it -- the
    endpoints and the path gradients are evaluated at the same precision.

    It is *not* below the 1e-5 the forward check reproduces the run's own
    ``pred_test.csv`` to, and it is pointless on the attention readout, which
    is one pass either way. Both of those run outside this block, in fp32, so
    turning it on cannot weaken the checks that decide whether the explanation
    is of the right model. On pre-Ampere hardware the flags are inert.
    """
    import torch
    if not enabled:
        yield
        return
    mm, cudnn = (torch.backends.cuda.matmul.allow_tf32,
                 torch.backends.cudnn.allow_tf32)
    torch.backends.cuda.matmul.allow_tf32 = True
    torch.backends.cudnn.allow_tf32 = True
    try:
        yield
    finally:
        torch.backends.cuda.matmul.allow_tf32 = mm
        torch.backends.cudnn.allow_tf32 = cudnn


def completeness(attr, f_x, f_base):
    """``(relative, absolute)`` shortfall of the attributions of a sample."""
    gap = np.abs(attr.sum(axis=1) - (f_x - f_base))
    span = max(float(np.abs(f_x - f_base).max()), 1e-12)
    return float(gap.max() / span), float(gap.max())


def calibrate_steps(fold, job, device, feats, abund, mask, base):
    """Smallest step count on the ladder whose completeness clears the tolerance.

    Measured on a handful of samples, which costs about a percent of the fold
    and decides a cost that scales with every sample and every member. 256 was
    the number a fixture needed; a fold that is happy with 64 should not pay
    four times that because of it, and one that needs 1024 should not quietly
    produce attributions that do not add up.

    The bar here is half of ``--ig-tol``, not ``--ig-tol``. The calibration
    sees a few samples and the gate afterwards sees all of them, so the worst
    sample of the fold is almost certainly not in the calibration set; leaving
    no margin means routinely choosing a step count that the full cohort then
    fails on, and paying for the calibration twice.

    Returns ``(n_steps, calibration_rows)``.
    """
    n = feats.shape[0]
    idx = np.unique(np.linspace(0, n - 1,
                                num=min(job["ig_cal_steps"], n)).astype(int))
    sub = (feats[idx], abund[idx], mask[idx], base[idx])

    trace = []
    for n_steps in IG_STEP_LADDER:
        attr, f_x, f_b, _ = integrated_gradients(
            fold["nets"], *sub, device=device, n_steps=n_steps,
            batch_rows=job["ig_batch"], target=job["target"],
            power=job["ig_power"])
        rel, _ = completeness(attr, f_x, f_b)
        trace.append((n_steps, rel))
        if rel <= job["ig_tol"] / 2:
            return n_steps, trace
    return IG_STEP_LADDER[-1], trace


# =====================================================================
# 5a. Attention: what the encoder looks at
# =====================================================================
def attention_profile(fold, feats, abund, mask, device, batch_size=32):
    """Each member's attention over the taxa of each sample.

    The question this answers -- which microbes the attention layer attends to
    -- does not need an attribution method, and gets one forward pass per
    sample per member instead of the integration's hundreds of forward *and*
    backward passes.

    ``attn_weights`` is ``(batch, n_heads, seq, seq)`` with the softmax over
    the **key** axis and masked keys driven to -1e9 before it, so each query
    row is already a distribution over the taxa the sample actually carries.
    The readout is the mean of those rows over the **unmasked queries**:

        A[b, j] = sum_i mask[b, i] * attn[b, ., i, j] / sum_i mask[b, i]

    which is how much the community as a whole attends to taxon j, and sums to
    1 over j. Masked query rows are dropped rather than averaged in: they are
    padding, they still emit a full distribution, and the pooling discards
    their output, so counting them would dilute every real row by however much
    padding the sample happens to carry -- that is, by its richness, which is a
    cohort-level batch effect.

    Averaging the query axis uniformly is not a simplification here. The model
    pools the encoder output as a plain mean over unmasked positions
    (no CLS token), so every query position contributes to the prediction with
    the same weight ``1 / n_valid``, and a weighted rollout would be that same
    uniform weight written out.

    Returns ``(per_member, n_valid)`` with `per_member` of shape
    ``(n_members, n_samples, seq)`` and `n_valid` the unmasked count per sample.
    """
    import torch

    n, steps = feats.shape
    out = np.zeros((len(fold["nets"]), n, steps), dtype=np.float64)
    n_valid = mask.numpy().sum(axis=1).astype(np.float64)

    with torch.no_grad():
        for k, net in enumerate(fold["nets"]):
            net.to(device)
            for i in range(0, n, batch_size):
                sl = slice(i, i + batch_size)
                mk = mask[sl].to(device)
                _, attn = net(feats[sl].to(device), abund[sl].to(device), mk)
                # Mean over heads, then over the unmasked query rows. Reduced
                # on the device: the (batch, heads, 600, 600) scores are the
                # largest tensor in the pass and only the (batch, 600) column
                # profile is worth moving back.
                q = mk.unsqueeze(1).unsqueeze(-1).to(attn.dtype)
                prof = (attn * q).sum(dim=2).mean(dim=1)
                prof = prof / mk.sum(dim=1, keepdim=True).clamp(min=1).to(prof.dtype)
                out[k, sl] = prof.cpu().numpy()
    return out, n_valid


def pooled_embedding(fold, feats, abund, mask, device, batch_size=32):
    """The sample representation the output head reads, meaned over members.

    ``forward(..., encoder=True)`` returns the mean, over the unmasked
    positions, of the encoder stack's output -- so after the attention *and*
    the feed-forward sub-layer, and before the head. That ``d_model``-wide
    vector is what the classifier actually sees of a sample, and it is what a
    samples-by-dimensions heatmap is a picture of.

    One forward per sample per member, no backward, like the attention
    readout; nothing here needs the integration.
    """
    import torch

    n = feats.shape[0]
    acc = None
    with torch.no_grad():
        for net in fold["nets"]:
            net.to(device)
            out = []
            for i in range(0, n, batch_size):
                sl = slice(i, i + batch_size)
                _, emb = net(feats[sl].to(device), abund[sl].to(device),
                             mask[sl].to(device), encoder=True)
                out.append(emb.cpu().numpy().astype(np.float64))
            e = np.concatenate(out)
            acc = e if acc is None else acc + e
    return acc / len(fold["nets"])


def check_attention(A, n_valid, X_positions, tol=1e-5):
    """The two identities the readout has to satisfy.

    A distribution over the taxa a sample carries: each row sums to 1, and
    nothing lands on a position the sample masked out. The second is what
    catches a mask used the wrong way round -- the rows would still sum to 1.
    """
    s = A.sum(axis=-1)
    worst_sum = float(np.abs(s - 1.0).max())
    leak = float(np.abs(A[..., ~X_positions]).max()) if (~X_positions).any() else 0.0
    if worst_sum > tol or leak > tol:
        raise ValueError(
            f"attention rows sum to 1 +/- {worst_sum:.2e} and put {leak:.2e} "
            f"on masked positions (tolerance {tol:.0e}); the query mask or the "
            f"softmax axis is not what this readout assumes")
    return dict(attn_row_sum_err=worst_sum, attn_mask_leak=leak,
                attn_n_valid_med=float(np.median(n_valid)))


def _spearman(a, b):
    """Spearman rho, indexed rather than attribute-accessed.

    ``.statistic`` only exists from SciPy 1.9; before that the result is a
    plain named tuple with ``.correlation``. The first element is the
    coefficient in every version, and this runs on whatever SciPy the cluster
    happens to have.
    """
    from scipy.stats import spearmanr
    return float(spearmanr(a, b)[0])


def write_attention(job, fold, A_mem, n_valid, feats, mask):
    """The attention outputs, and the two numbers needed to read them honestly.

    Written per fold:

    ``attn_{task}_{fold}.csv``
        ``n_samples x n_taxa`` **enrichment**, ``A * n_valid``, averaged over
        members. 1.0 means "attended exactly as much as an average taxon of
        this sample". The raw share is not what goes in the file: it is
        divided by however many taxa the sample carries, so ranking on it
        ranks sample richness -- a cohort-level batch effect -- as much as it
        ranks taxa.
    ``attn_{task}_{fold}_members.npz``
        The same, per member, for the spread across the ensemble.

    And into the run info:

    ``attn_abundance_spearman``
        How much of the attention ranking is just abundance. The attention
        logits are ``q.k/sqrt(d_k)`` with ``q`` and ``k`` both linear in
        ``embedding * abundance``, so the scores scale with abundance
        *squared* and the top of the ranking may be nothing but the sample's
        most abundant taxa. This puts a number on that rather than leaving it
        to be discovered by a reviewer. When it is high, the differential
        below is the only part worth reporting.
    """
    A = A_mem.mean(axis=0)
    enrich = scatter_to_taxa(fold, A * n_valid[:, None], feats, mask)

    _write_matrix(job, "attn", enrich, fold)
    if job["per_member"]:
        np.savez_compressed(
            os.path.join(job["out_dir"],
                         f"attn_{job['task']}_{job['fold']}_members.npz"),
            columns=fold["columns"], sample_ids=fold["sample_ids"],
            n_valid=n_valid,
            **{os.path.basename(d):
               scatter_to_taxa(fold, A_mem[k] * n_valid[:, None], feats,
                               mask).astype(np.float32)
               for k, d in enumerate(job["members"])})

    mean_attn = enrich.mean(axis=0)
    mean_abund = fold["X"].mean(axis=0)
    live = mean_abund > 0
    rho = (_spearman(mean_attn[live], mean_abund[live])
           if live.sum() > 2 else np.nan)
    return dict(attn_abundance_spearman=rho)


def scatter_to_taxa(fold, attr, feats, mask):
    """Positions -> the dense taxon space, zeros where a taxon is not present."""
    cols = fold["col_of_token"][feats.numpy()]
    keep = (cols >= 0) & (mask.numpy() > 0)

    out = np.zeros((attr.shape[0], len(fold["columns"])), dtype=np.float64)
    rows = np.repeat(np.arange(attr.shape[0]), attr.shape[1]).reshape(attr.shape)
    np.add.at(out, (rows[keep], cols[keep]), attr[keep])
    return out


# =====================================================================
# 6. Permutation SHAP, on the pruned feature space
# =====================================================================
def prune_columns(X, background):
    """Columns that can carry a non-zero attribution against this background.

    With a single background row, a feature equal to the background in a sample
    contributes exactly zero to that sample. A column equal to the background
    in *every* explained sample therefore has an attribution of zero
    everywhere, and dropping it is lossless rather than an approximation. On a
    sparse OTU table this is most of the table, and it is what decides whether
    the permutation explainer -- ``2 * n_features + 1`` forward passes per
    sample -- can be run at all.
    """
    keep = np.flatnonzero((X != background[None, :]).any(axis=0))
    return keep


def stratified_rows(labels, n, seed):
    """`n` row indices drawn per class in proportion, in the table's own order.

    Not the first `n` rows: a BIOM table's samples usually arrive grouped by
    cohort and often by class, so the head of a held-out cohort can be all
    controls, and an agreement statistic computed on one class only is not the
    one being claimed.
    """
    labels = np.asarray(labels)
    if n <= 0 or n >= len(labels):
        return np.arange(len(labels))
    rng = np.random.default_rng(seed)
    take = []
    for c in np.unique(labels):
        idx = np.flatnonzero(labels == c)
        k = max(1, int(round(n * len(idx) / len(labels))))
        take.append(rng.choice(idx, size=min(k, len(idx)), replace=False))
    return np.sort(np.concatenate(take))


def permutation_shap(fold, job, device, keep, candidates):
    """The old script's explainer, on this model and this feature space.

    Draws from `candidates` -- the samples the gradient method was run on --
    rather than from the whole cohort, so the two methods are compared on the
    same samples. Returns positions *within* `candidates`, which is what
    indexes the gradient attributions.
    """
    import shap

    X = fold["X"]
    bg = fold["background"]
    sel = stratified_rows(fold["labels"][candidates], job["perm_samples"],
                          job["perm_seed"])
    rows = candidates[sel]
    Xp = pd.DataFrame(X[rows][:, keep], columns=fold["columns"][keep])
    bgp = pd.DataFrame(bg[keep].reshape(1, -1), columns=fold["columns"][keep])

    def f(vals):
        # Chunked: the explainer hands over up to max_evals rows at once, and
        # widening those back to the full taxon space is an
        # ``n_rows x n_features`` float64 array -- gigabytes on a real table,
        # for a matrix that is then consumed a few hundred rows at a time
        # anyway.
        vals = np.asarray(vals)
        out = []
        for i in range(0, len(vals), job["perm_chunk"]):
            part = vals[i:i + job["perm_chunk"]]
            full = np.repeat(bg[None, :], len(part), axis=0)
            full[:, keep] = part
            out.append(predict(fold, full, device, job["batch_size"]))
        return np.concatenate(out)

    # Seeded: the permutation explainer draws its orderings at random, so two
    # runs of the same fold disagree by a little and the agreement statistic
    # this feeds would carry that noise on top of the difference between the
    # two methods, which is the thing being measured. `seed` is not in older
    # shap, and an unseeded run is worth far more than no run at all after the
    # hours this step costs -- so it degrades rather than raising.
    max_evals = job["perm_max_evals"] or (2 * len(keep) + 1)
    try:
        ex = shap.Explainer(f, bgp, seed=job["perm_seed"])
    except TypeError:
        print(f"[warn] {job['fold']}: this shap does not take seed=; the "
              f"permutation values will move a little between runs")
        ex = shap.Explainer(f, bgp)
    sv = ex(Xp, max_evals=max_evals)

    out = np.zeros((len(rows), len(fold["columns"])), dtype=np.float64)
    out[:, keep] = sv.values
    return out, sel, rows, max_evals


# =====================================================================
# 7. One fold, end to end
# =====================================================================
_GPU = None


def _init_worker(gpu_queue):
    """Bind this worker to one GPU for its lifetime, before torch is imported."""
    global _GPU
    _GPU = gpu_queue.get()
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    os.environ["CUDA_VISIBLE_DEVICES"] = str(_GPU)
    os.environ["MPLBACKEND"] = "Agg"
    os.environ.setdefault("OMP_NUM_THREADS", "4")


def explain_fold(job):
    """Returns ``(fold, summary | None, error)``."""
    try:
        import torch
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

        fold = load_fold(job)
        os.makedirs(job["out_dir"], exist_ok=True)
        info = dict(fold=job["fold"], disease=job["disease"],
                    study=job["study"], gpu=_GPU,
                    n_members=len(fold["nets"]),
                    n_samples=int(fold["X"].shape[0]),
                    n_features=int(len(fold["columns"])))

        if job["self_check"]:
            info.update(check_forward(fold, job, device))

        feats, abund, mask, base = _ig_inputs(fold, job)

        # --- what the attention layer looks at, over every held-out sample ---
        # First, and in fp32, and never subsampled: one forward per sample per
        # member is cheap enough that there is nothing to trade away here.
        A_mem, n_valid = attention_profile(fold, feats, abund, mask, device,
                                           job["batch_size"])
        if job["self_check"]:
            info.update(check_attention(A_mem, n_valid, mask.numpy() > 0))
        info.update(write_attention(job, fold, A_mem, n_valid, feats, mask))
        del A_mem

        if job["dump_pooled"]:
            P = pooled_embedding(fold, feats, abund, mask, device,
                                 job["batch_size"])
            # The label travels with the matrix. Everything downstream of this
            # file is a picture of samples, and a picture of samples needs to
            # know which are cases without having to find the metadata again.
            out = pd.DataFrame(P, index=fold["sample_ids"],
                               columns=[f"dim_{k}" for k in range(P.shape[1])])
            out.insert(0, "label", fold["labels"])
            path = os.path.join(job["out_dir"],
                                f"pool_{job['task']}_{job['fold']}.csv")
            out.to_csv(path)
            info["pooled_dims"] = int(P.shape[1])
            print(f"  [out ] {path}", flush=True)

        if job["skip_ig"]:
            # The cheap read-outs only. Everything below is the integration and
            # what compares against it.
            return job["fold"], info, None

        # --- and what moves the prediction, over a stratified subsample ---
        # The taxon ranking is a mean of |attribution| over samples, so it
        # converges long before the cohort is exhausted, and this is the one
        # cost here that scales with every sample *and* every member *and*
        # every path point. Verified rather than assumed: see --max-samples.
        rows = stratified_rows(fold["labels"], job["max_samples"],
                               job["sample_seed"])
        info["ig_samples"] = int(len(rows))
        ig_in = (feats[rows], abund[rows], mask[rows], base[rows])

        if job["ig_steps"]:
            n_steps, trace = job["ig_steps"], []
        else:
            n_steps, trace = calibrate_steps(fold, job, device, *ig_in)
        info["ig_steps"] = int(n_steps)
        info["ig_calibration"] = [[int(s), r] for s, r in trace]

        with tf32(job["tf32"]):
            attr, f_x, f_base, per_member = integrated_gradients(
                fold["nets"], *ig_in, device=device,
                n_steps=n_steps, batch_rows=job["ig_batch"],
                target=job["target"], power=job["ig_power"])

        # Completeness, now over every sample rather than the calibration's
        # handful. The attributions of a sample have to add up to the change in
        # the prediction they are attributions of; when they do not, some of
        # the model's response has been left unassigned and the per-taxon
        # numbers are a fraction of the story without saying which fraction.
        rel, absolute = completeness(attr, f_x, f_base)
        info["ig_completeness_rel_err"] = rel
        info["ig_completeness_abs_err"] = absolute
        if job["self_check"] and rel > job["ig_tol"]:
            raise ValueError(
                f"{job['fold']}: the attributions of the worst sample miss the "
                f"change they explain by {absolute:.4f} ({rel:.1%} of the "
                f"largest such change), over the tolerance {job['ig_tol']:.1%}"
                + (f". The calibration settled on {n_steps} steps from "
                   f"{len(trace)} samples and the full cohort disagrees, so "
                   f"raise --ig-cal-samples"
                   if trace else
                   f". Drop --ig-steps to let the step count be calibrated, or "
                   f"raise --ig-power (now {job['ig_power']:g}) -- the boundary "
                   f"layer this model has near the baseline is narrow rather "
                   f"than wide"))

        # Only where there is a branch to check. Recorded either way: a missing
        # check must not read as a passing one.
        lin = (check_linear_branch(fold, job, device, n_steps)
               if job["self_check"] else None)
        info["linear_branch"] = bool(
            getattr(fold["nets"][0], "linear_branch", False))
        if lin:
            info.update(lin)

        feats, mask = ig_in[0], ig_in[2]
        dense = scatter_to_taxa(fold, attr, feats, mask)
        _write_matrix(job, "shap_grad", dense, fold, rows=rows)

        # The members' disagreement, on record: the reported number is their
        # mean, and a taxon only one member leans on is a different claim from
        # one they all do. These came out of the same pass as the mean above.
        if job["per_member"]:
            np.savez_compressed(
                os.path.join(job["out_dir"],
                             f"shap_grad_{job['task']}_{job['fold']}"
                             f"_members.npz"),
                columns=fold["columns"], sample_ids=fold["sample_ids"],
                **{os.path.basename(d):
                   scatter_to_taxa(fold, per_member[k], feats,
                                   mask).astype(np.float32)
                   for k, d in enumerate(job["members"])})

        keep = prune_columns(fold["X"], fold["background"])
        info["n_features_informative"] = int(len(keep))

        if job["perm"]:
            t0 = time.time()
            perm, sel, perm_rows, max_evals = permutation_shap(
                fold, job, device, keep, rows)
            _write_matrix(job, "shap_perm", perm, fold, rows=perm_rows)
            # dense is indexed by the gradient subsample, and `sel` are
            # positions within it -- the two methods are being compared on the
            # same samples, not on two draws that happen to be the same size.
            info.update(_agreement(dense[sel], perm, fold["columns"]))
            info.update(perm_samples=len(perm_rows), perm_max_evals=max_evals,
                        perm_minutes=round((time.time() - t0) / 60, 1))

        return job["fold"], info, None

    except Exception:
        return job["fold"], None, traceback.format_exc()


def _write_matrix(job, tag, values, fold, rows=None):
    """`tag` is the whole prefix, so the attention files are not called shap_*."""
    path = os.path.join(job["out_dir"],
                        f"{tag}_{job['task']}_{job['fold']}.csv")
    idx = fold["sample_ids"] if rows is None else fold["sample_ids"][rows]
    pd.DataFrame(values, index=idx, columns=fold["columns"]).to_csv(path)


def _agreement(grad, perm, columns, top_k=30):
    """How far the cheap method and the reference method actually agree."""
    g = np.abs(grad).mean(axis=0)
    p = np.abs(perm).mean(axis=0)
    live = (g != 0) | (p != 0)
    rho = _spearman(g[live], p[live]) if live.sum() > 2 else np.nan
    gt = set(columns[np.argsort(-g)[:top_k]])
    pt = set(columns[np.argsort(-p)[:top_k]])
    return dict(spearman_grad_vs_perm=rho,
                top30_overlap=len(gt & pt),
                pearson_signed=float(np.corrcoef(grad.ravel(), perm.ravel())[0, 1]))


# =====================================================================
# 8. Aggregation across folds
# =====================================================================
def attention_summary(jobs, args):
    """Per-fold attention rankings, and the case-minus-control differential.

    The differential is the one to read. Mean attention is dominated by
    abundance -- the scores scale with it squared -- and by whatever taxa are
    simply common, both of which are the same in cases and controls and cancel
    out of a difference. "The model attends to this taxon more when the sample
    is a case" is a statement about the disease; "the model attends to this
    taxon" is largely a statement about sequencing depth.
    """
    rows = []
    for j in jobs:
        path = os.path.join(j["out_dir"],
                            f"attn_{j['task']}_{j['fold']}.csv")
        if not os.path.isfile(path):
            continue
        df = pd.read_csv(path, index_col=0)
        # Labels from a member's own prediction file, reindexed onto the
        # attention matrix rather than assumed to be in the same order.
        pred = os.path.join(j["members"][0], "pred_test.csv")
        if not os.path.isfile(pred):
            continue
        lab = pd.read_csv(pred).set_index("sample_id")["true_label"]
        y = lab.reindex(df.index).to_numpy()
        if np.isnan(y).any() or len(np.unique(y)) < 2:
            continue
        mean_all = df.mean(axis=0)
        live = mean_all[mean_all > 0].index
        case = df.loc[y == 1, live].mean(axis=0)
        ctrl = df.loc[y == 0, live].mean(axis=0)
        diff = (case - ctrl)
        order = mean_all[live].sort_values(ascending=False)
        rows.append(pd.DataFrame(dict(
            fold=j["fold"], disease=j["disease"], study=j["study"],
            otu=order.index,
            enrichment=order.values,
            enrichment_case=case[order.index].values,
            enrichment_control=ctrl[order.index].values,
            case_minus_control=diff[order.index].values,
            rank=np.arange(1, len(order) + 1),
            rank_diff=(-diff.abs()).rank(method="first")[order.index].values)))
    if not rows:
        return None

    summary = pd.concat(rows, ignore_index=True)
    out = os.path.join(args.out_dir, f"attn_summary_{args.task}.csv")
    summary.to_csv(out, index=False)
    print(f"\n[out] per-fold attention rankings -> {out}")

    for disease, part in summary.groupby("disease"):
        n_folds = part.fold.nunique()
        top = part[part["rank_diff"] <= args.top_k]
        cons = (top.groupby("otu")
                   .agg(n_folds_in_top=("fold", "nunique"),
                        mean_rank_diff=("rank_diff", "mean"),
                        case_minus_control=("case_minus_control", "mean"),
                        enrichment=("enrichment", "mean")))
        # Ranked on the size of the shift, not its direction: a taxon the model
        # attends to far less in cases is as much of a finding as one it
        # attends to far more.
        cons = (cons.assign(_abs=cons.case_minus_control.abs())
                    .sort_values(["n_folds_in_top", "_abs"],
                                 ascending=[False, False])
                    .drop(columns="_abs"))
        if n_folds > 1:
            cons["frac_folds"] = cons.n_folds_in_top / n_folds
        else:
            cons = cons.drop(columns=["n_folds_in_top"])
        p = os.path.join(args.out_dir,
                         f"attn_consistency_{args.task}_{disease}.csv")
        cons.to_csv(p)
        print(f"[out] {disease} attention: {n_folds} fold(s) -> {p}")
        print(cons.head(15).to_string(float_format=lambda v: f"{v:.4f}"))
    return summary


def attention_vs_gradient(jobs, args, top_k=30):
    """Does the model look at the taxa that move its prediction?

    Two different questions -- attention is what the encoder reads, the
    gradient is what the output responds to -- and how far apart they land is
    worth a number. Compared on the samples the gradient was run on, which are
    a subset of the attention's.
    """
    out = []
    for j in jobs:
        pa = os.path.join(j["out_dir"], f"attn_{j['task']}_{j['fold']}.csv")
        pg = os.path.join(j["out_dir"], f"shap_grad_{j['task']}_{j['fold']}.csv")
        if not (os.path.isfile(pa) and os.path.isfile(pg)):
            continue
        a = pd.read_csv(pa, index_col=0)
        g = pd.read_csv(pg, index_col=0)
        a = a.loc[g.index]                       # the gradient's subsample
        av, gv = a.mean(axis=0).to_numpy(), g.abs().mean(axis=0).to_numpy()
        live = (av > 0) | (gv > 0)
        if live.sum() < 3:
            continue
        cols = a.columns.to_numpy()
        ta = set(cols[np.argsort(-av)[:top_k]])
        tg = set(cols[np.argsort(-gv)[:top_k]])
        out.append(dict(fold=j["fold"], disease=j["disease"],
                        spearman_attn_vs_grad=_spearman(av[live], gv[live]),
                        top30_overlap=len(ta & tg)))
    if not out:
        return None
    df = pd.DataFrame(out)
    p = os.path.join(args.out_dir, f"attn_vs_grad_{args.task}.csv")
    df.to_csv(p, index=False)
    print(f"\n[out] attention vs gradient rankings -> {p}")
    print(f"   spearman: median {df.spearman_attn_vs_grad.median():.3f}, "
          f"range {df.spearman_attn_vs_grad.min():.3f}.."
          f"{df.spearman_attn_vs_grad.max():.3f}; "
          f"top30 overlap median {df.top30_overlap.median():.0f}/30")
    print("   Low is not a bug: what an encoder reads and what its output "
          "responds to are different questions. Report both, and say which "
          "figure is which.")
    return df


def summarise(jobs, infos, args):
    """Per-fold taxon rankings, and how many cohorts each taxon survives."""
    rows = []
    for j in jobs:
        path = os.path.join(j["out_dir"],
                            f"shap_grad_{j['task']}_{j['fold']}.csv")
        if not os.path.isfile(path):
            continue
        df = pd.read_csv(path, index_col=0)
        mean_abs = df.abs().mean(axis=0)
        live = mean_abs[mean_abs > 0]
        order = live.sort_values(ascending=False)
        rows.append(pd.DataFrame(dict(
            fold=j["fold"], disease=j["disease"], study=j["study"],
            otu=order.index, mean_abs=order.values,
            mean_signed=df[order.index].mean(axis=0).values,
            rank=np.arange(1, len(order) + 1))))
    if not rows:
        return None

    summary = pd.concat(rows, ignore_index=True)
    out = os.path.join(args.out_dir, f"shap_summary_{args.task}.csv")
    summary.to_csv(out, index=False)
    print(f"\n[out] per-fold taxon rankings -> {out}")

    # Cross-cohort consistency. One fold's ranking is one cohort's opinion;
    # what is worth reporting is the taxa that several held-out cohorts
    # independently put near the top.
    for disease, part in summary.groupby("disease"):
        n_folds = part.fold.nunique()
        top = part[part["rank"] <= args.top_k]
        cons = (top.groupby("otu")
                   .agg(n_folds_in_top=("fold", "nunique"),
                        mean_rank=("rank", "mean"),
                        mean_abs=("mean_abs", "mean"),
                        mean_signed=("mean_signed", "mean"))
                   .sort_values(["n_folds_in_top", "mean_abs"],
                                ascending=[False, False]))
        if n_folds > 1:
            cons["frac_folds"] = cons.n_folds_in_top / n_folds
        else:
            # Under lodo a disease is one fold, so n_folds_in_top would be 1
            # and frac_folds 1.0 for every taxon in the table -- columns that
            # look like agreement across cohorts and are arithmetic. Drop them
            # and say what the ranking actually is.
            cons = cons.drop(columns=["n_folds_in_top"])
        p = os.path.join(args.out_dir,
                         f"shap_consistency_{args.task}_{disease}.csv")
        cons.to_csv(p)
        print(f"[out] {disease}: {n_folds} fold(s) -> {p}"
              + ("" if n_folds > 1 else
                 "  (one fold, so this is that fold's ranking and not a "
                 "cross-cohort agreement)"))
        print(cons.head(15).to_string(float_format=lambda v: f"{v:.4f}"))

    agree = [i for i in infos if "spearman_grad_vs_perm" in i]
    if agree:
        p = os.path.join(args.out_dir,
                         f"shap_method_agreement_{args.task}.csv")
        pd.DataFrame(agree).to_csv(p, index=False)
        print(f"\n[out] gradient vs permutation -> {p}")
        for i in agree:
            print(f"   {i['fold']}: spearman {i['spearman_grad_vs_perm']:.3f}, "
                  f"top30 overlap {i['top30_overlap']}/30, "
                  f"{i['perm_minutes']} min on {i['perm_samples']} samples")
        worst = min(i["spearman_grad_vs_perm"] for i in agree)
        if worst < 0.6:
            print(f"   WARNING: the weakest agreement is {worst:.3f}. The "
                  f"fixed-tokenization approximation is not holding on this "
                  f"data; the gradient rankings should not be reported on "
                  f"their own.")
    return summary


# =====================================================================
# 9. Main
# =====================================================================
def main():
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--task", default="disease", choices=list(TASKS))
    ap.add_argument("--run-name", default="results_with_SNEs_ckpt",
                    help="the --run-name of the run being explained; its "
                         "member directories must still hold model.pth, which "
                         "means it was run with --keep-ckpt")
    ap.add_argument("--run-tsv", default=None,
                    help="fold list; defaults to the task's own")
    ap.add_argument("--diseases", nargs="*", default=None,
                    help="explain only these diseases. Which run-TSV column "
                         "that is comes from the task's own fold template: "
                         "'disease' under --task disease, 'left_out_disease' "
                         "under lodo. Default: every fold in the TSV")
    ap.add_argument("--out-dir", default="Data/biomark")
    ap.add_argument("--gpus", type=int, nargs="+", default=[0])
    ap.add_argument("--workers-per-gpu", type=int, default=1)

    ap.add_argument("--target", default="prob", choices=["prob", "logit"],
                    help="what is attributed. 'prob' (default) is the mean of "
                         "the members' probabilities, which is what "
                         "--report-combiner prob reports, so the attributions "
                         "sum to a difference in that number. 'logit' is "
                         "additive in the model's own output and is what the "
                         "linear-branch check runs in")
    ap.add_argument("--ig-steps", type=int, default=0,
                    help="points on the integration path. 0 (default) "
                         f"calibrates it per fold over {list(IG_STEP_LADDER)}, "
                         "taking the first that clears --ig-tol on a few "
                         "samples. The cost of this scales with every sample "
                         "and every member, so a fold that is happy with 64 "
                         "should not pay for 256; one that needs 1024 should "
                         "not quietly report attributions that do not add up")
    ap.add_argument("--ig-cal-samples", type=int, default=8,
                    dest="ig_cal_steps",
                    help="samples the step calibration measures on")
    ap.add_argument("--ig-power", type=float, default=3.0,
                    help="clusters the path points near the baseline, as "
                         "alpha = u**power. 1 spaces them evenly and does not "
                         "work here: most of the model's response to abundance "
                         "happens inside alpha < 0.01, where the near-zero "
                         "baseline puts the encoder's LayerNorm into its "
                         "epsilon, and an even spacing leaves a completeness "
                         "error that no step count removes. 3 (default) closed "
                         "it to zero on the fixture")
    ap.add_argument("--ig-batch", type=int, default=256,
                    help="path points evaluated per forward/backward pass. 64 "
                         "rows of a 600-position sequence is a small kernel "
                         "for a one-layer d_model=100 encoder, and the GPU "
                         "spends most of its time on launch overhead; raise "
                         "this as far as memory allows")
    ap.add_argument("--max-samples", type=int, default=300,
                    help="stratified cap on the held-out samples the "
                         "*integration* explains, per fold. The taxon ranking "
                         "is a mean of |attribution| over samples and settles "
                         "well before a cohort is exhausted, while the cost "
                         "scales with samples x members x path points. 0 for "
                         "no cap. The attention readout ignores this and "
                         "always covers every sample -- it is one forward pass "
                         "per sample and has nothing to trade away")
    ap.add_argument("--sample-seed", type=int, default=0,
                    help="seeds the --max-samples draw")
    ap.add_argument("--tf32", dest="tf32", action="store_true", default=True,
                    help="allow TF32 matmuls inside the integration on Ampere "
                         "and later (default). The forward check, the "
                         "attention readout and the prediction function stay "
                         "in fp32, so this cannot weaken the checks that "
                         "decide whether the right model is being explained")
    ap.add_argument("--no-tf32", dest="tf32", action="store_false")
    ap.add_argument("--ig-tol", type=float, default=0.02,
                    help="how far the attributions of a sample may fall short "
                         "of the change they explain, as a fraction of the "
                         "largest such change in the fold, before the fold is "
                         "failed rather than written out")
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--num-steps", type=int, default=CFG["num_steps"],
                    help="OTU positions per sample; must match the run")
    ap.add_argument("--n-heads", type=int, default=CFG["n_heads"],
                    help="attention heads. The only shape the checkpoint does "
                         "not pin down, since the projections are d_model x "
                         "d_model whatever it was")
    ap.add_argument("--abund-mode", default="multiply",
                    choices=["multiply", "none"],
                    help="must match the run being explained. Not recorded in "
                         "the checkpoint, and getting it wrong changes what "
                         "the abundance does without changing any shape")
    ap.add_argument("--dump-pooled", action="store_true",
                    help="also write pool_{task}_{fold}.csv: the pooled "
                         "encoder output per held-out sample, meaned over the "
                         "members, with the sample's label in the first "
                         "column. That is the d_model-wide vector the output "
                         "head reads -- after the attention and the "
                         "feed-forward sub-layer -- and what a "
                         "samples-by-dimensions heatmap is drawn from")
    ap.add_argument("--skip-ig", action="store_true",
                    help="stop after the attention readout and --dump-pooled, "
                         "skipping the integration and everything that "
                         "compares against it. Both of those are one forward "
                         "pass per sample per member, so this turns a run of "
                         "hours into one of minutes when the integration's "
                         "output is already on disk")
    ap.add_argument("--median-block", type=int, default=2048,
                    help="columns per block when the background median is "
                         "computed. Bounds that step's peak memory at "
                         "n_samples x this, which is what keeps lodo's pooled "
                         "table from being stacked dense")
    ap.add_argument("--check-background", action="store_true",
                    help="assert the blocked background median against the "
                         "definition, element for element. Worth a run on a "
                         "small fold whenever union_median is touched; it "
                         "materialises the dense stack the blocking exists to "
                         "avoid, so not on a lodo fold")
    ap.add_argument("--no-per-member", action="store_true",
                    help="skip writing the per-member attributions. They come "
                         "out of the same pass as the ensemble's and cost "
                         "nothing extra to compute, so this only saves the "
                         "npz")

    ap.add_argument("--perm-per-disease", type=int, default=1,
                    help="also run the permutation explainer on this many "
                         "folds per disease, smallest held-out cohort first, "
                         "to measure what the gradient method's fixed "
                         "tokenization costs. 0 disables it")
    ap.add_argument("--perm-max-folds", type=int, default=2,
                    help="cap on how many folds get the permutation explainer "
                         "in total, smallest held-out cohort first. Applied "
                         "after --perm-per-disease, which is not itself a cap "
                         "under lodo: a disease there *is* one fold, so 'one "
                         "per disease' would be every fold in the run")
    ap.add_argument("--perm-folds", nargs="*", default=None,
                    help="run the permutation explainer on these folds by "
                         "name instead of choosing them by size")
    ap.add_argument("--perm-samples", type=int, default=0,
                    help="explain a stratified draw of N held-out samples with "
                         "the "
                         "permutation explainer; 0 (default) does all of them")
    ap.add_argument("--perm-seed", type=int, default=0,
                    help="seeds the permutation explainer's orderings, so the "
                         "agreement it is being used to measure does not move "
                         "between runs")
    ap.add_argument("--perm-chunk", type=int, default=256,
                    help="rows the permutation explainer's masked inputs are "
                         "widened back to the full taxon space in. Bounds the "
                         "peak memory of that step, which is otherwise "
                         "max_evals x n_features")
    ap.add_argument("--perm-max-evals", type=int, default=0,
                    help="0 (default) uses 2 * n_informative_features + 1, "
                         "which is what the explainer needs to be exact. A "
                         "smaller number makes it an approximation and it will "
                         "say so")

    ap.add_argument("--no-self-check", action="store_true",
                    help="skip reproducing the run's own pred_test.csv and the "
                         "linear-branch identity. Only for debugging: those "
                         "two are what make the attributions trustworthy")
    ap.add_argument("--check-tol", type=float, default=1e-5,
                    help="tolerance on reproducing pred_test.csv")
    ap.add_argument("--top-k", type=int, default=30,
                    help="rank cut-off used by the cross-cohort consistency "
                         "table")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if args.skip_ig and (args.perm_folds or args.perm_per_disease > 0):
        # The permutation explainer is scored against the integration's own
        # attributions, and --skip-ig means there are none to score against.
        print("[warn] --skip-ig: the permutation explainer is turned off with "
              "it, since the agreement it computes is against the integration")
        args.perm_per_disease, args.perm_folds = 0, None

    jobs = build_jobs(args)
    choose_perm_folds(jobs, args)
    os.makedirs(args.out_dir, exist_ok=True)

    n_perm = sum(1 for j in jobs if j["perm"])
    print(f"[cfg] task={args.task} run_name={args.run_name} "
          f"target={args.target} "
          f"ig_steps={args.ig_steps or 'calibrated per fold'}")
    print(f"[cfg] folds named by {'/'.join(fold_fields(args.task))}; "
          f"--diseases filters on {fold_fields(args.task)[0]!r}")
    print(f"[cfg] {len(jobs)} fold(s), "
          f"{sum(len(j['members']) for j in jobs)} member(s) total; "
          f"{n_perm} fold(s) also get the permutation explainer")
    for j in jobs:
        print(f"   {j['fold']:<28} {len(j['members'])} members"
              + ("   <- permutation" if j["perm"] else ""))
    if args.dry_run:
        return

    slots = [g for g in args.gpus for _ in range(args.workers_per_gpu)]
    ctx = mp.get_context("spawn")
    q = ctx.Queue()
    for s in slots:
        q.put(s)

    infos, failed, t0 = [], [], time.time()
    with ctx.Pool(processes=len(slots), initializer=_init_worker,
                  initargs=(q,)) as pool:
        for i, (fold, info, err) in enumerate(
                pool.imap_unordered(explain_fold, jobs), 1):
            el = (time.time() - t0) / 60
            if err is None:
                infos.append(info)
                extra = ""
                if "spearman_grad_vs_perm" in info:
                    extra = f" perm_rho={info['spearman_grad_vs_perm']:.3f}"
                pd_ = info.get("max_abs_prob_diff")
                comp = info.get("ig_completeness_rel_err")
                print(f"[{i}/{len(jobs)}] {fold}  "
                      + (f"pred_diff={pd_:.1e} " if pd_ is not None
                         else "pred_diff=unchecked ")
                      + (f"steps={info['ig_steps']} "
                         f"completeness={comp:.1e}"
                         + ("" if info.get("linear_branch")
                            else " lin_branch=off(check n/a)")
                         if comp is not None else "attn+pool only")
                      + f"{extra}  | {el:.0f}m", flush=True)
            else:
                failed.append(fold)
                print(f"[{i}/{len(jobs)}] {fold}  FAILED", flush=True)
                print("   ", err.strip().splitlines()[-1], flush=True)

    if infos:
        p = os.path.join(args.out_dir, f"shap_run_info_{args.task}.json")
        with open(p, "w") as f:
            json.dump(infos, f, indent=2, default=float)
    if failed:
        print(f"\n{len(failed)} fold(s) failed: {', '.join(failed)}")

    ok = [j for j in jobs if j["fold"] not in failed]
    attention_summary(ok, args)
    summarise(ok, infos, args)
    attention_vs_gradient(ok, args)

    conf = [i.get("attn_abundance_spearman") for i in infos]
    conf = [c for c in conf if c is not None and not np.isnan(c)]
    if conf:
        print(f"\n[attn] attention vs abundance, spearman over folds: "
              f"median {np.median(conf):.3f}, max {max(conf):.3f}")
        if np.median(conf) > 0.7:
            print("   The attention ranking largely tracks abundance -- the "
                  "scores scale with it squared. Report the case-minus-control "
                  "differential rather than the raw ranking; the shared part "
                  "cancels out of it.")


if __name__ == "__main__":
    main()
