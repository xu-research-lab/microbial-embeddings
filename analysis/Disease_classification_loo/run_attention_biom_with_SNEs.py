#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Bagged ensemble of Attention_biom across all LOSO tasks, on multiple GPUs.

One protocol, the random-forest one
-----------------------------------
Each fold trains --n-estimators members. A member differs from its siblings in
two ways, and both are the ways a random forest's trees differ from each other:

  bagging          The member trains on a random --bag-frac of the fold's
                   training table and validates on the rest. The draw is plain
                   random, taking no account of study or label, which is what a
                   forest's bootstrap does; the samples a member did not draw
                   are exactly the samples it early-stops on -- its out-of-bag
                   set, used here to choose an epoch and a threshold rather
                   than to estimate error.

  feature subsets  The member sees a random --feat-frac of the taxa that are
                   non-zero in its own training part; taxa that are zero
                   throughout it are kept whole, since they are not features it
                   could learn from and dropping one only costs the vocabulary
                   an SNE vector the test cohort can use. This needs nothing
                   from the model: the vocabulary is built from the training
                   table it is handed, so a dropped taxon is simply absent, and
                   the held-out samples resolve it to '<unk>' and mask it out
                   -- which is what a
                   tree that never saw a feature does with it.

Both draws are made here rather than in the library, which takes the training
and validation cohorts as two separate BIOM tables and no longer partitions
anything itself. Each member therefore writes its own pair of tables into its
job directory, and they are deleted once it succeeds unless --keep-tables.

The fold's result is the ensemble of its members' test predictions, not any
single member. A member trained on 63% of the data and 80% of the taxa is
weaker than one trained on everything; --n-estimators is what buys that back,
which is why a forest grows hundreds of trees rather than one.

--bag-frac defaults to 1 - 1/e = 0.6321, the fraction of *distinct* samples a
bootstrap of the same size contains. Drawing without replacement at that rate
matches a bootstrap's diversity while avoiding what with-replacement would cost
here: BIOM sample ids must be unique, so duplicates would need renaming, and a
renamed sample no longer finds its metadata row.

Combining the members
---------------------
AUC reads only the ordering of the scores, and the members are not on a common
scale -- each one's probabilities come from a threshold fitted on its own
out-of-bag set. Five combiners are therefore scored side by side, all recomputed
from the same pred_test.csv files at no cost:

  ens_logit   mean of log(p/(1-p)); the geometric mean of the odds, and how
              independent evidence adds. Most exposed to one confident member.
  ens_prob    mean of p.
  ens_zlogit  mean of the logits after standardising each member over the test
              cohort, which removes that member's mean and spread but keeps
              the distances inside it.
  ens_trim    mean of the logits with the top and bottom 20% of members cut
              per sample. Between the mean and the median.
  ens_median  median of the logits; a degenerate member cannot move it.
  ens_rank    mean of each member's within-cohort ranks. Wholly immune to the
              members disagreeing about scale.
  ens_normal  ranks through the normal quantile (van der Waerden scores):
              scale-free like rank, but with the tails spaced apart.
  ens_rankmed median of the ranks -- scale-free and outlier-resistant at once.
  ens_wauc    logits weighted by valid_auc - 0.5. Left NaN under
              --inner-split loso, where each member was scored on a different
              cohort and the weights would not be comparable.

The middle group is there for what this data actually does: members calibrated
on different cohorts produce probabilities on different scales, and a mean of
raw logits trusts those scales while a mean of ranks discards them.

Parallelism
-----------
One job = one Attention_biom call, and every member of every fold and task is
queued together rather than fold by fold: at the end of a fold the pool would
otherwise idle waiting for its last member. Each worker binds one GPU through
CUDA_VISIBLE_DEVICES and uses numb=0 inside, so memory is isolated per process
and one job running out of it cannot disturb the others.

Re-running: by default a job whose output directory exists has that directory
wiped and is run again, so a second run of the same command produces fresh
results rather than echoing back what is on disk. Pass --resume to keep
finished jobs and run only the missing ones -- what you want after a crash, or
to recompute the combiners from predictions that are already there.

Usage
-----
  python run_attention_biom_with_SNEs.py --dry-run
  python run_attention_biom_with_SNEs.py --tasks disease --n-estimators 100 \
      --gpus 0 1 2 3 4 5 6 7 --run-name rf_bag
  python run_attention_biom_with_SNEs.py --tasks all --n-estimators 100 \
      --bag-frac 0.8 --feat-frac 0.8 --run-name rf_bag80
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import multiprocessing as mp
import os
import shutil
import time
import traceback
from collections import Counter

import biom
import numpy as np
import pandas as pd
from biom.util import biom_open

# torch and membed are imported lazily inside the worker functions. Importing
# them at module level would let spawned children load torch before the pool
# initializer sets CUDA_VISIBLE_DEVICES, which breaks the per-process GPU
# binding.


# =====================================================================
# 0. Task table: the four parts differ only in path templates
# =====================================================================
TASKS = {
    "disease": dict(
        run_tsv="run_leave_one_study_out_each_diease_list.tsv",
        train="Data/disease_data/{disease}/{study}/train_loo.biom",
        test="Data/disease_data/{disease}/{study}/test_loo.biom",
        meta="Data/disease_data/{disease}/metadata.tsv",
        fold="{disease}_{study}",
        artifact="Data/disease_data/",
        results="disease.csv",
    ),
    "ibd_subtype": dict(
        run_tsv="run_leave_one_study_out_IBD_subtype_list.tsv",
        train="Data/IBD_subtype_data/{disease}/{study}/train_loo.biom",
        test="Data/IBD_subtype_data/{disease}/{study}/test_loo.biom",
        meta="Data/IBD_subtype_data/{disease}/metadata.tsv",
        fold="{disease}_{study}",
        artifact="Data/IBD_subtype_data/",
        results="ibd_subtype.csv",
    ),
    "lodo": dict(
        run_tsv="run_leave_one_disease_out.tsv",
        train="Data/loo_all_diseases/data/{left_out_disease}/train_loo.biom",
        test="Data/loo_all_diseases/data/{left_out_disease}/test_loo.biom",
        meta="Data/loo_all_diseases/data/metadata.tsv",
        fold="{left_out_disease}",
        artifact="Data/loo_all_diseases/",
        results="lodo.csv",
    ),
    "loso_all": dict(
        run_tsv="run_leave_one_study_out.tsv",
        train="Data/loo_all_studies/data/{left_out_study}/train_loo.biom",
        test="Data/loo_all_studies/data/{left_out_study}/test_loo.biom",
        meta="Data/loo_all_studies/data/metadata.tsv",
        fold="{left_out_study}",
        artifact="Data/loo_all_studies/",
        results="loso_all.csv",
    ),
}

SAMPLE_ID_COL = "sample"
LABELS_COL = "group"

CFG = dict(
    num_steps=600,
    p_drop=0.4,
    batch_size=64,
    d_model=100,
    n_layers=1,
    n_heads=1,
    d_ff=200,
    lr=1e-3,
    # A single int: the reverted library builds one hidden layer,
    # Linear(d_model, head_hidden) -> GELU -> Dropout -> Linear(head_hidden, 1).
    # A list here reaches nn.Linear as out_features and raises.
    head_hidden=64,
    weight_decay=1e-4,
    num_epochs=100,
    loss="BCEWithLogits",
)
GLOVE_DEFAULT = "../../data/social_niche_embedding_removing_disease_samples_100.txt"


def parse_set(pairs):
    """``['p_drop=0.3', 'd_ff=400']`` -> ``{'p_drop': 0.3, 'd_ff': 400}``.

    Only keys already in `CFG` are accepted, and each value is read with the
    type of the default it replaces. A sweep is the one place a silent typo is
    most expensive -- it produces a full run that looks like a result and is
    the baseline -- so both halves are refused rather than guessed at.
    """
    out = {}
    for item in pairs:
        if "=" not in item:
            raise SystemExit(f"--set wants KEY=VALUE, got {item!r}")
        key, _, raw = item.partition("=")
        key, raw = key.strip(), raw.strip()
        if key not in CFG:
            near = [k for k in CFG if k.startswith(key[:3])]
            raise SystemExit(
                f"--set {key!r} is not a CFG entry"
                + (f"; did you mean {' or '.join(near)}?" if near else "")
                + f"\navailable: {', '.join(sorted(CFG))}")
        default = CFG[key]
        try:
            if isinstance(default, bool):
                if raw.lower() not in ("true", "false", "1", "0"):
                    raise ValueError(raw)
                out[key] = raw.lower() in ("true", "1")
            elif isinstance(default, int):
                out[key] = int(raw)
            elif isinstance(default, float):
                out[key] = float(raw)
            else:
                out[key] = raw
        except ValueError:
            raise SystemExit(
                f"--set {key}={raw!r}: the default is {default!r}, so this "
                f"has to read as {type(default).__name__}")
    return out

#: Fraction of distinct samples a bootstrap of the same size contains.
BAG_FRAC_DEFAULT = 1 - 1 / math.e

#: How many training samples a feature subset may empty before the member is
#: refused. A few emptied rows are the cost of random subspaces on a sparse
#: table; past this share the subset is too small to be training on.
EMPTY_SAMPLE_LIMIT = 0.05


# =====================================================================
# 1. Materialising a member's bag as two BIOM tables
# =====================================================================
def _write_table(tbl, path):
    """Write a BIOM table atomically.

    Via a temporary name and os.replace, so a job killed mid-write cannot leave
    a truncated table that the next run would happily load as if it were
    complete. The temporary name carries the pid, so two members writing at the
    same moment cannot collide.
    """
    tmp = f"{path}.tmp.{os.getpid()}"
    with biom_open(tmp, "w") as f:
        tbl.to_hdf5(f, "run_attention_biom")
    os.replace(tmp, path)


def survey_training_table(train_biom, meta_path, study_col):
    """The cohorts present in a fold's training table."""
    sids = biom.load_table(train_biom).ids(axis="sample")
    md = pd.read_csv(meta_path, sep="\t", index_col=SAMPLE_ID_COL,
                     dtype={SAMPLE_ID_COL: str}, low_memory=False)
    return md.loc[sids, study_col].unique()


def split_train_val_with_joint_stratification(pool, sample_to_study,
                                              sample_to_label, test_size,
                                              random_state):
    """Split with a (study, label) stratified split, degrading as needed.

    Three attempts in order -- the joint ``study||label`` stratum, then the
    label alone, then a plain random split. Each fallback is taken only when a
    stratum would carry a single sample, which is where ``train_test_split``
    cannot place one on each side.

    Returns ``(train_ids, val_ids, strategy)``. The strategy is worth
    recording: a fold that fell through to ``'random'`` has a validation set
    whose class balance is not controlled, and one that fell to
    ``'label_only'`` has cohorts spread across the split by chance.
    """
    from sklearn.model_selection import train_test_split
    joint, labels = [], []
    for sid in pool:
        study = sample_to_study.get(sid, "UNKNOWN_STUDY")
        lab = sample_to_label.get(sid)
        lab = "UNKNOWN_LABEL" if lab is None or str(lab).strip() == "" else str(lab)
        labels.append(lab)
        joint.append(f"{study}||{lab}")

    for strata, name in ((joint, "joint_study_label"), (labels, "label_only")):
        if min(Counter(strata).values()) < 2:
            logging.warning("a %s stratum has fewer than two samples; "
                            "degrading", name)
            continue
        try:
            tr, va = train_test_split(pool, test_size=test_size, shuffle=True,
                                      random_state=random_state,
                                      stratify=strata)
            return tr, va, name
        except ValueError as exc:
            logging.warning("%s stratification failed: %s", name, exc)
    tr, va = train_test_split(pool, test_size=test_size, shuffle=True,
                              random_state=random_state)
    return tr, va, "random"


def _valid_mask(job, sids, lab, md):
    """Boolean mask over `sids` (True = validation), plus a strategy label."""
    mode = job["inner_mode"]

    if mode == "bag":
        # A plain random draw, which is what a forest's bootstrap is: it takes
        # no account of study or label. The samples not drawn are the member's
        # out-of-bag set and become its validation table.
        rng = np.random.default_rng(job["bag_seed"])
        n_bag = max(1, int(round(len(sids) * job["bag_frac"])))
        mask = np.ones(len(sids), dtype=bool)
        mask[rng.permutation(len(sids))[:n_bag]] = False
        return mask, "random_bag"

    if mode == "joint":
        # One split, stratified on (study, label), shared by every member --
        # they differ only in model_seed. The validation samples come from the
        # same cohorts as the training ones, so this selects for generalization
        # within the training distribution rather than across cohorts.
        pool = [str(x) for x in sids]
        _, val_ids, strategy = split_train_val_with_joint_stratification(
            pool,
            md.loc[sids, job["study_col"]].astype(str).to_dict(),
            md.loc[sids, LABELS_COL].to_dict(),
            job["val_ratio"], job["split_seed"])
        return np.isin(sids, np.asarray(val_ids, dtype=object)), strategy

    if mode == "same_disease":
        return _same_disease_mask(job, sids, lab, md)

    if mode == "disease_loso":
        # One study held out per disease, redrawn per member. Plain 'loso'
        # holds out a single cohort, so its validation set speaks for whichever
        # disease that cohort happens to carry -- on the pooled tasks the epoch
        # is then chosen on one disease out of a dozen. Here every disease
        # contributes one held-out cohort, so the selected epoch is the one
        # that transfers across studies *for all diseases at once*, which is
        # what the pooled test set actually asks for.
        return _disease_loso_mask(job, sids, md)

    # leave-one-study-out: this member holds out one or more whole training
    # cohorts (job["valid_studies"]), so the epoch it selects is the epoch
    # that transfers across studies rather than within one.
    g = md.loc[sids, job["study_col"]].astype(str).to_numpy()
    valid_studies = [str(s) for s in job["valid_studies"]]
    return np.isin(g, valid_studies), "loso"


def _column(md, sids, col, what):
    """One metadata column over `sids` as strings, refusing blanks.

    A blank would survive ``astype(str)`` as the literal ``'nan'`` and become a
    bucket of its own -- a disease that does not exist, drawing a held-out
    study of its own. Cheaper to refuse here than to explain the extra stratum
    later.
    """
    if col not in md.columns:
        raise KeyError(
            f"metadata has no {what} column {col!r}; available columns: "
            f"{list(md.columns)}")
    v = md.loc[sids, col]
    blank = v.isna() | (v.astype(str).str.strip() == "")
    if blank.any():
        raise ValueError(
            f"{what} column {col!r} is empty for {int(blank.sum())} of "
            f"{len(sids)} training samples; first offenders: "
            f"{v.index[blank][:5].tolist()}")
    return v.astype(str).to_numpy()


def _same_disease_mask(job, sids, lab, md):
    """Hold out one study of the test study's disease; keep the rest training.

    Leave-one-study-out confined to the disease under test. The fold holds out
    study S of disease D and asks whether the model reaches it; this member
    holds out one of D's *other* studies, so the epoch it keeps is the one that
    crossed studies within D -- the same step the test cohort asks for. D's
    remaining studies stay in training, so the model still learns the disease.

    That is the difference from ``disease_loso``, whose validation set spans
    every disease in the pool and therefore scores a transfer the test does not
    measure. It is also the difference from validating on all of D at once,
    which would leave no D in training at all.

    One member per study of D, so every one of them is validated on exactly
    once and the members' predictions are worth averaging. ``--n-estimators``
    does not apply, as under ``loso``.

    When D has a single study left in the pool, that member holds it out and D
    is absent from its training part. The member still runs; its strategy says
    ``same_disease_only_study`` so those folds can be read apart afterwards.

    Parameters
    ----------
    job : dict
        Needs ``valid_studies`` -- the one study this member holds out, chosen
        in :func:`build_jobs` from the test disease's cohorts -- and
        ``study_col``.
    sids : numpy.ndarray
        Sample IDs of the fold's training table.
    md : pandas.DataFrame
        Metadata indexed by sample ID.

    Returns
    -------
    mask : numpy.ndarray of bool
        True where the sample belongs in validation.
    strategy : str
        ``'same_disease'``, or ``'same_disease_only_study'`` when holding this
        study out leaves none of the disease in training.
    """
    g = _column(md, sids, job["study_col"], "study")
    dis = _column(md, sids, job["disease_col"], "disease")
    valid_studies = [str(x) for x in job["valid_studies"]]
    in_study = np.isin(g, valid_studies)
    tag = f"{job['fold']}/{job['valid_unit']}"

    if job.get("cv_fold") is not None:
        # The disease has this one study, so holding all of it out would leave
        # the model with none of the disease it is about to be tested on.
        # Cut the study k ways instead and let this member validate on one
        # part: the other 80% stays in training, and across the k members
        # every sample is validated on exactly once.
        from sklearn.model_selection import StratifiedKFold
        idx = np.flatnonzero(in_study)
        k = int(job["cv_k"])
        skf = StratifiedKFold(n_splits=k, shuffle=True,
                              random_state=job["split_seed"])
        _, val_local = list(skf.split(idx, lab[idx]))[int(job["cv_fold"])]
        mask = np.zeros(len(sids), dtype=bool)
        mask[idx[val_local]] = True
        print(f"[sdis] {tag}: {', '.join(sorted(set(dis[in_study])))} has only "
              f"{', '.join(valid_studies)}, so it is cut {k} ways -- this "
              f"member validates on part {int(job['cv_fold']) + 1}/{k} "
              f"({int(mask.sum())} samples) and trains on the other "
              f"{int(in_study.sum()) - int(mask.sum())} plus every other "
              f"disease")
        return mask, "same_disease_cv"

    mask = in_study
    if not mask.any():
        raise ValueError(
            f"{tag}: study {valid_studies} holds no sample of the training "
            f"table, so validation would be empty")

    held_dis = sorted(set(dis[mask]))
    left = {d: int(((dis == d) & ~mask).sum()) for d in held_dis}
    stranded = [d for d, n in left.items() if n == 0]
    print(f"[sdis] {tag}: validating on {', '.join(valid_studies)} "
          f"({', '.join(held_dis)}), {int(mask.sum())}/{len(sids)} samples "
          f"({mask.sum() / len(sids):.1%}); "
          + ", ".join(f"{d} keeps {n} training samples"
                      for d, n in sorted(left.items())))
    if stranded:
        print(f"[sdis] {tag}: {', '.join(stranded)} had this one study only, "
              f"so it is absent from this member's training part")
        return mask, "same_disease_only_study"
    return mask, "same_disease"


def _disease_loso_mask(job, sids, md):
    """Hold out one randomly drawn study per disease, for this member.

    The unit held out is the ``(disease, study)`` cell, not the whole study. A
    study that spans several diseases -- PRJEB11419 carries OB, IBS, BD and
    T2DM -- therefore contributes only its samples of the disease it was drawn
    for, and keeps the rest in training. That keeps each disease represented in
    validation by exactly one cohort, at the cost of letting one study sit on
    both sides of the split, so a validation AUC from this mode is not a clean
    cross-study estimate for the diseases that share a cohort.

    Every disease is drawn from, including one with a single study in the
    training table -- holding that cohort out removes the disease from training
    altogether. The draw says so on its output line rather than skipping it.

    The draw uses ``job['bag_seed']``, which is ``--bag-seed`` plus the member
    index, so the members differ in which cohorts they hold out and their
    predictions are worth ensembling.
    """
    rng = np.random.default_rng(job["bag_seed"])
    dis = _column(md, sids, job["disease_col"], "disease")
    stu = _column(md, sids, job["study_col"], "study")

    mask = np.zeros(len(sids), dtype=bool)
    picks, lonely = [], []
    # sorted() so the sequence of draws depends on the seed alone, not on the
    # order the diseases happen to appear in the table.
    for d in sorted(set(dis)):
        in_d = dis == d
        studies = sorted(set(stu[in_d]))
        chosen = studies[int(rng.integers(len(studies)))]
        mask |= in_d & (stu == chosen)
        picks.append(f"{d}/{chosen}" + ("*" if len(studies) == 1 else ""))
        if len(studies) == 1:
            lonely.append(d)

    tag = f"{job['fold']}/{job['valid_unit']}"
    share = mask.sum() / len(sids)
    print(f"[dloso] {tag}: {len(picks)} diseases -> {' '.join(picks)} | "
          f"{int(mask.sum())}/{len(sids)} samples ({share:.0%}) held out for "
          f"validation")
    if share > 0.5:
        # One cohort per disease is a large slice when the diseases have two or
        # three cohorts each, which is what this data mostly has. Not an error
        # -- the members ensemble, so it trades the same way bagging does -- but
        # a member training on under half its fold is worth knowing about.
        print(f"[dloso] {tag}: WARNING: that leaves only {1 - share:.0%} of the "
              f"fold to train on; the diseases here have few cohorts each")
    if lonely:
        print(f"[dloso] {tag}: * = {', '.join(lonely)} had one study only, so "
              f"that disease is now absent from this member's training part")
    return mask, "disease_loso"


def survey_same_disease_studies(train_biom, meta_path, study_col, disease_col,
                                test_study):
    """Training-table studies carrying the disease(s) of the held-out study."""
    sids = biom.load_table(train_biom).ids(axis="sample")
    md = pd.read_csv(meta_path, sep="\t", index_col=SAMPLE_ID_COL,
                     dtype={SAMPLE_ID_COL: str}, low_memory=False)
    test_mask = md[study_col].astype(str) == str(test_study)
    if not test_mask.any():
        raise SystemExit(
            f"no metadata row has {study_col}=={test_study!r}; --inner-split "
            f"same_disease expects the fold name to be the held-out study id, "
            f"which is how loso_all names its folds")
    test_dis = set(md.loc[test_mask, disease_col].dropna().astype(str))
    sub = md.loc[sids]
    keep = sub[disease_col].astype(str).isin(test_dis)
    sub = sub.loc[keep]
    cohorts = sorted(set(sub[study_col].astype(str)))
    lab = pd.to_numeric(sub[LABELS_COL], errors="coerce")
    counts = {c: (int(((sub[study_col].astype(str) == c) & (lab == 1)).sum()),
                  int(((sub[study_col].astype(str) == c) & (lab == 0)).sum()))
              for c in cohorts}
    return cohorts, sorted(test_dis), counts


def split_paths(out_dir):
    return (os.path.join(out_dir, "train_part.biom"),
            os.path.join(out_dir, "valid_part.biom"))


def member_split(job):
    """Materialise this member's training and validation parts as BIOM tables.

    The library no longer splits the training table itself, so the partition is
    made here and handed over as files. Which partition depends on
    --inner-split; see `_valid_mask`.

    Returns ``(train_path, valid_path, strategy)``. Existing files are reused,
    so a retried job -- or, under 'joint', every member after the first -- does
    not pay for the filtering again.

    --feat-frac is applied to the training part only, and only to the taxa that
    are non-zero somewhere in it. The validation and test tables keep every
    taxon: the vocabulary is built from the training part, so a dropped taxon is
    absent from it and the held-out samples resolve it to '<unk>' and mask it
    out, which is what a tree that never saw a feature does with it. Training
    samples the subset empties are dropped, up to EMPTY_SAMPLE_LIMIT of them;
    beyond that the member is refused.
    """
    split_dir = job["split_dir"]
    train_part, valid_part = split_paths(split_dir)
    info_path = os.path.join(split_dir, "split.json")
    if all(os.path.exists(x) for x in (train_part, valid_part, info_path)):
        with open(info_path) as f:
            return train_part, valid_part, json.load(f)["strategy"]
    os.makedirs(split_dir, exist_ok=True)

    table = biom.load_table(job["train"])
    sids = np.asarray(table.ids(axis="sample"))
    md = pd.read_csv(job["meta"], sep="\t", index_col=SAMPLE_ID_COL,
                     dtype={SAMPLE_ID_COL: str}, low_memory=False)
    missing = [x for x in sids if x not in md.index]
    if missing:
        raise KeyError(f"{len(missing)} training samples have no metadata row; "
                       f"first offenders: {missing[:5]}")
    # Drop whole diseases from this fold before anything is split off. Removing
    # them here rather than from the training part alone is deliberate: a
    # disease that is not trained on has no business in the validation set
    # either, where it would be scored as an unseen disease and drag the epoch
    # selection around. The test table is a separate file and is untouched, so
    # a dropped disease is still evaluated on its own fold.
    if job["drop_disease"]:
        dis = _column(md, sids, job["disease_col"], "disease")
        keep = ~np.isin(dis, [str(d) for d in job["drop_disease"]])
        n_drop = int((~keep).sum())
        if n_drop:
            if not keep.any():
                raise ValueError(
                    f"--drop-disease {job['drop_disease']} removed every "
                    f"sample of {job['fold']}'s training table")
            gone = sorted(set(dis[~keep]))
            print(f"[drop] {job['fold']}/{job['valid_unit']}: removed "
                  f"{n_drop}/{len(sids)} samples of {', '.join(gone)} from "
                  f"this fold's training and validation pool")
            table = table.filter(sids[keep], axis="sample", inplace=False)
            sids = np.asarray(table.ids(axis="sample"))

    lab = pd.to_numeric(md.loc[sids, LABELS_COL], errors="coerce")
    if lab.isna().any():
        bad = lab.index[lab.isna()][:5].tolist()
        raise ValueError(f"non-numeric labels in {LABELS_COL!r}: {bad}")
    lab = lab.to_numpy().astype(int)

    is_v, strategy = _valid_mask(job, sids, lab, md)

    # Nothing in a plain random or leave-one-study-out draw guarantees both
    # classes land on both sides. A part that lost a class cannot be trained on
    # or early-stopped against, so the member fails here rather than training
    # on a meaningless validation curve.
    for name, m in (("training", ~is_v), ("validation", is_v)):
        if m.sum() == 0:
            hint = ""
            if job["inner_mode"] == "disease_loso" and name == "training":
                hint = (f"; every disease was held out whole, which is what "
                        f"happens when --inner-group ({job['study_col']!r}) and "
                        f"--disease-col ({job['disease_col']!r}) name the same "
                        f"column -- each disease is then its own only cohort")
            raise ValueError(f"the {name} part of {job['valid_unit']} is empty "
                             f"under --inner-split {job['inner_mode']}{hint}")
        if len(np.unique(lab[m])) < 2:
            raise ValueError(
                f"the {name} part of {job['valid_unit']} has only one class "
                f"under --inner-split {job['inner_mode']}")

    train_tbl = table.filter(sids[~is_v], axis="sample", inplace=False)
    if job["feat_frac"] < 1.0:
        obs = np.asarray(train_tbl.ids(axis="observation"))
        # The tables are aligned on feature IDs, so the training part carries a
        # row for every taxon in the study -- including taxa that are zero
        # across all of its samples and only appear in the held-out cohort.
        # Those rows are not features to subsample: they teach the member
        # nothing, and dropping one costs the vocabulary an SNE vector that the
        # test cohort does have a use for. Draw the subset from the taxa the
        # member can actually learn from, and keep the all-zero rows whole.
        present = np.asarray((train_tbl.matrix_data > 0).sum(axis=1)).ravel() > 0
        pool, always = obs[present], obs[~present]
        k = max(1, int(round(len(pool) * job["feat_frac"])))
        drawn = np.random.default_rng(job["bag_seed"] + 1).choice(
            pool, size=k, replace=False)
        print(f"[feat] {job['valid_unit']}: {len(pool)} taxa present in this "
              f"member's training part -> kept {k} ({job['feat_frac']:.0%}); "
              f"{len(always)} taxa absent from it kept whole")
        keep = np.concatenate([drawn, always])
        train_tbl = train_tbl.filter(keep, axis="observation", inplace=False)
        # A sample whose every retained taxon is zero has nothing left: it
        # encodes to all padding, is masked out end to end, and contributes a
        # row of noise to the loss. Drop those samples rather than the member.
        # Losing the member would shrink K on that fold alone and leave the
        # ensemble uneven across folds, which is worse than losing a few rows
        # -- and on a sparse table at a low --feat-frac it is near-certain.
        kept_sids = np.asarray(train_tbl.ids(axis="sample"))
        nz = np.asarray((train_tbl.matrix_data > 0).sum(axis=0)).ravel()
        empty = int((nz == 0).sum())
        if empty:
            if empty > EMPTY_SAMPLE_LIMIT * len(nz):
                raise ValueError(
                    f"--feat-frac {job['feat_frac']} left {empty} of "
                    f"{len(nz)} training samples with no taxa at all "
                    f"({empty / len(nz):.0%} > "
                    f"{EMPTY_SAMPLE_LIMIT:.0%}); raise --feat-frac")
            survivors = kept_sids[nz > 0]
            lab_kept = lab[~is_v][nz > 0]
            if len(np.unique(lab_kept)) < 2:
                raise ValueError(
                    f"--feat-frac {job['feat_frac']} emptied every sample of "
                    f"one class in {job['valid_unit']}; raise --feat-frac")
            print(f"[feat] {job['valid_unit']}: --feat-frac "
                  f"{job['feat_frac']} emptied {empty}/{len(nz)} training "
                  f"samples; dropped them")
            train_tbl = train_tbl.filter(survivors, axis="sample",
                                         inplace=False)

    _write_table(train_tbl, train_part)
    _write_table(table.filter(sids[is_v], axis="sample", inplace=False),
                 valid_part)
    with open(info_path, "w") as f:
        json.dump(dict(strategy=strategy, n_train=int((~is_v).sum()),
                       n_valid=int(is_v.sum()),
                       dropped=list(job["drop_disease"])), f)
    return train_part, valid_part, strategy

# =====================================================================
# 2. Worker: bind a GPU, run one job
# =====================================================================
_GPU = None


def _init_worker(gpu_queue):
    """Each worker takes one GPU from the queue and keeps it for its lifetime."""
    global _GPU
    _GPU = gpu_queue.get()
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    os.environ["CUDA_VISIBLE_DEVICES"] = str(_GPU)
    os.environ["MPLBACKEND"] = "Agg"      # Animator calls savefig; no GUI backend
    os.environ.setdefault("OMP_NUM_THREADS", "4")


def _build_kwargs(job, train_biom, valid_biom):
    kw = dict(
        metadata=job["meta"], train_biom=train_biom, valid_biom=valid_biom,
        test_biom=job["test"],
        embedding_birnn=os.path.join(job["out_dir"], "model.pth"),
        plotfile_loss=os.path.join(job["out_dir"], "loss.png"),
        plotfile_auc=os.path.join(job["out_dir"], "auc.png"),
        labels_col=LABELS_COL, sample_id_col=SAMPLE_ID_COL,
        pred_out=os.path.join(job["out_dir"], "pred"),
        numb=0,                            # only one GPU visible in this process
        linear_branch=job["linear_branch"],
        lin_weight_decay=job["lin_weight_decay"],
        glove_embedding=job["glove"],
        select_by=job["select_by"],
        abund_mode=job["abund_mode"],
        train_sampler=job["train_sampler"],
        valid_auc=job["valid_auc"],
        min_group_n=job["min_group_n"],
        logit_adjust_tau=job["logit_adjust_tau"],
        disease_col=job["disease_col"],
        max_ratio=job["max_ratio"],
        sample_seed=job["sample_seed"],
        **CFG,
    )
    if job["head_hidden"] is not None:
        kw["head_hidden"] = job["head_hidden"]
    kw.update(job["override"])
    return kw


def _is_job_dir(path):
    """Whether `path` sits under a ``members/`` directory.

    A cheap sanity check before anything recursive is deleted: every job
    directory this script builds has that component, so a mistyped artifact
    root or an empty fold name cannot end up handing rmtree something like
    ``Data/`` itself.
    """
    return "members" in os.path.normpath(path).split(os.sep)


def _prepare_out_dir(job):
    """Create the job directory, wiping any previous run of the same job first."""
    out_dir = job["out_dir"]
    if job["overwrite"] and os.path.isdir(out_dir):
        if not _is_job_dir(out_dir):
            raise RuntimeError(
                f"refusing to wipe {out_dir!r}: it does not look like a job "
                f"directory (no 'members' path component)")
        shutil.rmtree(out_dir)
    os.makedirs(out_dir, exist_ok=True)


def run_job(job):
    """Run one Attention_biom call. Returns (job_id, result | None, error)."""
    from membed.otu_attention import Attention_biom

    res_path = os.path.join(job["out_dir"], "result.json")
    try:
        # Wipes the directory first unless --resume. Inside the try: the guard
        # in _prepare_out_dir raises, and an exception escaping run_job would
        # take down imap_unordered and with it the entire run.
        _prepare_out_dir(job)

        if not job["overwrite"] and os.path.exists(res_path):
            with open(res_path) as f:
                return job["id"], json.load(f), None

        train_part, valid_part, strategy = member_split(job)
        rec, tm = Attention_biom(**_build_kwargs(job, train_part, valid_part))
        ep = rec["epoch"]
        if not isinstance(ep, int):
            raise RuntimeError(f"unexpected epoch value {ep!r}")
        out = dict(task=job["task"], fold=job["fold"],
                   member=job["member"], valid_unit=str(job["valid_unit"]),
                   inner_mode=job["inner_mode"], strategy=strategy,
                   epoch=int(ep),
                   valid_auc=float(rec["valid_auc"]),
                   n_valid=int(len(rec["valid_label"])),
                   test_auc=float(tm["auc"]),
                   bag_frac=job["bag_frac"], feat_frac=job["feat_frac"],
                   abund_mode=job["abund_mode"],
                   out_dir=job["out_dir"], gpu=_GPU)

        with open(res_path, "w") as f:
            json.dump(out, f)
        if not job["keep_ckpt"]:
            ck = os.path.join(job["out_dir"], "model.pth")
            if os.path.exists(ck):
                os.remove(ck)
        if not job["keep_tables"] and job["split_dir"] == job["out_dir"]:
            # Only after result.json exists, so a crash before this point
            # leaves the split in place for the retry to reuse. A shared split
            # -- 'joint', where every member reads the same pair -- is left
            # alone here: the first member to finish would be deleting a table
            # its siblings still need. The parent drops it after the pool.
            for ft in split_paths(job["out_dir"]):
                if os.path.exists(ft):
                    os.remove(ft)
        return job["id"], out, None

    except Exception:
        # Do not let one bad job kill the whole run. No result.json is written,
        # so re-running the script retries this job automatically.
        err = traceback.format_exc()
        try:
            os.makedirs(job["out_dir"], exist_ok=True)
            with open(os.path.join(job["out_dir"], "failed.txt"), "w") as f:
                f.write(err)
        except Exception:
            pass
        return job["id"], None, err

    finally:
        # Animator opens two figures per train_cls call and never closes them.
        # Workers live for the whole run, so without this the leak grows
        # unbounded and matplotlib starts warning after ten jobs.
        try:
            import matplotlib.pyplot as plt
            plt.close("all")
        except Exception:
            pass


# =====================================================================
# 3. Job construction (parent process, never touches CUDA)
# =====================================================================
def build_jobs(task, args):
    """One job per (fold, member)."""
    t = TASKS[task]
    if args.inner_split == "disease_loso" and task not in ("lodo", "loso_all"):
        raise SystemExit(
            f"--inner-split disease_loso needs a training table spanning "
            f"several diseases, and {task!r} is single-disease (its metadata is "
            f"one file per disease), where the mode would degrade to drawing "
            f"one random study. Use --tasks lodo and/or loso_all, or pick "
            f"--inner-split loso for {task!r}")
    run = pd.read_csv(t["run_tsv"], sep="\t")
    if args.drop_disease:
        # Check the names against the real metadata once, here, rather than
        # letting a typo quietly drop nothing across a whole run.
        probe = t["meta"].format(**run.iloc[0].to_dict())
        known = set(pd.read_csv(probe, sep="\t", low_memory=False)
                    [args.disease_col].dropna().astype(str).unique())
        bad = [d for d in args.drop_disease if str(d) not in known]
        if bad:
            raise SystemExit(
                f"--drop-disease {bad} not found in column "
                f"{args.disease_col!r} of {probe}; known values: "
                f"{sorted(known)}")
    jobs, folds = [], []

    # min_delta 0 under 'auc': any strictly higher AUC counts as a new best, so
    # the checkpoint lands on the true peak and the patience countdown starts
    # there. A positive threshold would let a genuinely better epoch be scored
    # as "no improvement".
    mdelta = 0.0 if args.select_by in ("auc", "last") else 0.001
    cfg_over = parse_set(args.set_cfg)

    for _, r in run.iterrows():
        row = r.to_dict()
        fold = t["fold"].format(**row)
        train = t["train"].format(**row)
        test = t["test"].format(**row)
        meta = t["meta"].format(**row)

        # loso only: unit -> studies it holds out, and unit -> a shared split
        # directory when --loso-k groups several repeats onto the same split.
        group_studies, group_split_dir, cv_folds = {}, {}, {}

        if args.inner_split == "loso":
            cohorts = sorted(str(x) for x in
                             survey_training_table(train, meta, args.inner_group))
            if args.loso_k is None:
                # One member per training cohort. The member count is the
                # fold's -- a fold whose training table spans eight cohorts
                # trains eight members and one spanning four trains four --
                # so no flag caps it.
                units = cohorts
                for u in units:
                    group_studies[u] = [u]
                print(f"[info] {task}/{fold}: {len(units)} cohort(s) -> "
                      f"leave-one-study-out, {len(units)} members")
            else:
                if args.loso_k >= len(cohorts):
                    raise SystemExit(
                        f"--loso-k {args.loso_k} leaves no training cohort "
                        f"for {task}/{fold}, which has only {len(cohorts)} "
                        f"cohort(s)")
                # One shuffle, cut into non-overlapping groups of --loso-k --
                # every cohort in the fold is validated on exactly once.
                # --n-estimators then repeats each group with a different
                # model_seed, sharing that group's split (see split_dir below).
                rng = np.random.default_rng(args.loso_k_seed)
                shuffled = list(cohorts)
                rng.shuffle(shuffled)
                chunks = [shuffled[i:i + args.loso_k]
                         for i in range(0, len(shuffled), args.loso_k)]
                units = [f"g{gi}_m{ri}" for gi in range(len(chunks))
                        for ri in range(args.n_estimators)]
                for gi, chunk in enumerate(chunks):
                    split_dir = os.path.join(t["artifact"], fold, "members",
                                             f"_split_g{gi}")
                    for ri in range(args.n_estimators):
                        u = f"g{gi}_m{ri}"
                        group_studies[u] = chunk
                        group_split_dir[u] = split_dir
                print(f"[info] {task}/{fold}: {len(cohorts)} cohort(s) -> "
                      f"{len(chunks)} group(s) of up to {args.loso_k}, "
                      f"x{args.n_estimators} repeat(s) = {len(units)} members")
        elif args.inner_split == "same_disease":
            cohorts, test_dis, ccount = survey_same_disease_studies(
                train, meta, args.inner_group, args.disease_col, fold)
            if not cohorts:
                raise SystemExit(
                    f"{task}/{fold}: the held-out study's disease(s) "
                    f"{test_dis} have no other study in the training table, "
                    f"so there is nothing to validate on. Drop this fold or "
                    f"use --inner-split disease_loso")
            if len(cohorts) > 1:
                units = cohorts
                for u in units:
                    group_studies[u] = [u]
                print(f"[info] {task}/{fold}: test carries "
                      f"{', '.join(test_dis)} -> {len(units)} member(s), one "
                      f"per study of it")
            else:
                # One study left: holding all of it out would strip the
                # disease from training entirely. Cut it k ways instead, so
                # 1-1/k of it stays in training and each member validates on a
                # different part. k is capped by the smaller class, since a
                # stratified fold cannot be made without one of each.
                only = cohorts[0]
                n_pos, n_neg = ccount[only]
                k = min(args.cv_folds, n_pos, n_neg)
                if k < 2:
                    units = cohorts
                    group_studies[only] = [only]
                    print(f"[info] {task}/{fold}: {', '.join(test_dis)} has "
                          f"only {only} and it holds {n_pos} case/{n_neg} "
                          f"control, too few to cut up; holding the study out "
                          f"whole, so this member trains without the disease")
                else:
                    units = [f"{only}_cv{i}" for i in range(k)]
                    for i, u in enumerate(units):
                        group_studies[u] = [only]
                        cv_folds[u] = (i, k)
                    print(f"[info] {task}/{fold}: {', '.join(test_dis)} has "
                          f"only {only} ({n_pos} case/{n_neg} control) -> "
                          f"{k}-fold split of it, {k} member(s), each "
                          f"validating on {100 / k:.0f}% and training on the "
                          f"rest")
        else:
            units = [f"m{i}" for i in range(args.n_estimators)]
        folds.append(dict(task=task, fold=fold, artifact=t["artifact"],
                          n_members=len(units), inner_mode=args.inner_split))

        # Under 'joint' every member reads the same split, so it is written
        # once beside their job directories rather than copied into each.
        # Under grouped 'loso' each group similarly shares one split across
        # its --n-estimators repeats, which differ only in model_seed.
        shared = os.path.join(t["artifact"], fold, "members", "_split")

        for i, u in enumerate(units):
            out_dir = os.path.join(t["artifact"], fold, "members", str(u))
            if args.inner_split == "joint":
                split_dir = shared
            elif u in group_split_dir:
                split_dir = group_split_dir[u]
            else:
                split_dir = out_dir
            jobs.append(dict(
                id=f"{task}|{fold}|{u}", task=task, fold=fold, member=i,
                inner_mode=args.inner_split, valid_unit=u,
                split_dir=split_dir,
                valid_studies=group_studies.get(u),
                cv_fold=cv_folds.get(u, (None, None))[0],
                cv_k=cv_folds.get(u, (None, None))[1],
                study_col=args.inner_group,
                val_ratio=args.val_ratio, split_seed=args.split_seed,
                out_dir=out_dir, train=train, test=test, meta=meta,
                head_hidden=args.head_hidden,
                linear_branch=args.linear_branch,
                lin_weight_decay=args.lin_weight_decay,
                glove=args.glove_embedding, select_by=args.select_by,
                abund_mode=args.abund_mode,
                train_sampler=args.train_sampler,
                valid_auc=args.valid_auc,
                drop_disease=list(args.drop_disease),
                min_group_n=args.min_group_n,
                logit_adjust_tau=args.logit_adjust_tau,
                disease_col=args.disease_col,
                max_ratio=args.max_ratio,
                sample_seed=args.sample_seed,
                bag_frac=args.bag_frac, feat_frac=args.feat_frac,
                bag_seed=args.bag_seed + i,
                keep_ckpt=args.keep_ckpt, keep_tables=args.keep_tables,
                overwrite=args.overwrite,
                override=dict(model_seed=args.ensemble_seed0 + i,
                              patience=args.patience, min_delta=mdelta,
                              **cfg_over)))
    return jobs, folds


def cleanup_shared_splits(jobs):
    """Delete the 'joint' shared tables once every member has had its turn.

    A per-member split is removed by the worker that owns it, but a shared one
    has no single owner: whichever member finished first would be deleting a
    table its siblings still have to read. ``split.json`` stays -- it is small,
    records which strategy ran, and its absence is what tells a later --resume
    to rebuild.
    """
    n = 0
    for d in sorted({j["split_dir"] for j in jobs
                     if j["split_dir"] != j["out_dir"]}):
        for p in split_paths(d):
            if os.path.exists(p):
                os.remove(p)
                n += 1
    if n:
        print(f"[split] removed {n} shared split table(s); "
              f"pass --keep-tables to keep them")


# =====================================================================
# 4. Parallel execution
# =====================================================================
def run_pool(jobs, gpus, workers_per_gpu, label):
    if not jobs:
        print(f"[{label}] nothing to run")
        return []
    slots = [g for g in gpus for _ in range(workers_per_gpu)]
    ctx = mp.get_context("spawn")     # CUDA requires spawn; fork corrupts context
    q = ctx.Queue()
    for s in slots:
        q.put(s)

    print(f"[{label}] {len(jobs)} jobs across {len(slots)} slots "
          f"(GPUs {gpus} x {workers_per_gpu})")
    results, failed, t0 = [], [], time.time()
    with ctx.Pool(processes=len(slots), initializer=_init_worker,
                  initargs=(q,)) as pool:
        for i, (jid, out, err) in enumerate(
                pool.imap_unordered(run_job, jobs), 1):
            el = time.time() - t0
            eta = el / i * (len(jobs) - i)
            if err is None:
                results.append(out)
                print(f"[{label} {i}/{len(jobs)}] {jid}  "
                      f"AUC={out['test_auc']:.4f} epoch={out['epoch']} "
                      f"gpu={out['gpu']}  | elapsed {el/60:.0f}m "
                      f"eta {eta/60:.0f}m", flush=True)
            else:
                failed.append(jid)
                print(f"[{label} {i}/{len(jobs)}] {jid}  FAILED", flush=True)
                print("   ", err.strip().splitlines()[-1], flush=True)
    if failed:
        print(f"\n[{label}] {len(failed)} job(s) failed:")
        for j in failed[:20]:
            print("   ", j)
        print("  see failed.txt in each job directory; re-running this script "
              "retries them automatically")
    return results


# =====================================================================
# 5. Combining the members
# =====================================================================
#: Combiners scored side by side. AUC reads only the ordering of the scores, so
#: what separates these is how each handles the members not being on a common
#: scale: member i's probabilities come from a threshold fitted on member i's
#: own out-of-bag set.
ENSEMBLE_COMBINERS = ("logit", "prob", "rank", "normal", "zlogit",
                      "trim", "median", "rankmed", "wauc")

#: Trimmed fraction cut from each end of the member axis by ``'trim'``.
TRIM_FRAC = 0.2


def ensemble_aucs(pred_paths, weights=None, eps=1e-6):
    """Test AUC of each combiner over a fold's member predictions.

    Parameters
    ----------
    pred_paths : list of str
        One ``pred_test.csv`` per member. All must score the same samples in
        the same order; a mismatch returns NaNs rather than a quietly
        misaligned average.
    weights : array-like, optional
        Per-member weight for ``'wauc'``, normally ``valid_auc - 0.5``.
        Negative entries are clipped to zero, so a member that did worse than
        chance on its own out-of-bag set does not vote.

    Returns
    -------
    dict
        ``{name: auc}`` over `ENSEMBLE_COMBINERS`, NaN where the combiner could
        not be formed.

    Notes
    -----
    ``'logit'`` is the geometric mean of the odds, which is how independent
    evidence adds, and the combiner most exposed to one confident member: the
    clip bounds a member's contribution at ``log(eps / (1 - eps))``, about
    13.8, which without it would let a single ``p = 0.9999999`` decide a
    sample.

    ``'rank'`` throws the probabilities away and averages each member's
    ordering, which is the only combiner wholly immune to the members
    disagreeing about scale. Since AUC scores an ordering anyway it gives up
    nothing the metric uses -- except the distinction between a member that
    separates the classes widely and one that barely does.

    ``'median'`` is the one that survives a degenerate member: a model that
    early-stopped at epoch 1, or one whose bag was unlucky, still votes under a
    mean and cannot move a median.

    The middle four exist for the case this task actually has -- members whose
    probability distributions differ because they were calibrated on different
    cohorts. They sit between ``'logit'``, which trusts the scales, and
    ``'rank'``, which discards them entirely:

    ``'zlogit'``
        Each member's logits are standardised over the test cohort before
        averaging, which removes that member's mean and spread but keeps the
        relative distances inside it. What ``'rank'`` throws away and
        ``'logit'`` is hurt by is exactly the difference between members here.
    ``'normal'``
        Ranks mapped through the normal quantile function -- van der Waerden
        scores. As shift- and scale-free as ``'rank'``, but the spacing grows
        towards the tails, so a sample every member puts first is separated
        from the second more than two mid-ranked samples are. Plain rank
        averaging treats those gaps as equal.
    ``'trim'``
        The mean of the logits after cutting `TRIM_FRAC` of the members from
        each end, per sample. Between ``'logit'`` and ``'median'``: outlier
        members cannot pull it, but the magnitudes of the rest still count.
    ``'rankmed'``
        The median of the ranks -- scale-free and outlier-resistant at once.

    ``'wauc'`` needs the members' validation scores to be comparable, which is
    true when they share a validation set or draw it the same way, and false
    under leave-one-study-out where each was scored on a different cohort. The
    caller decides by passing or withholding `weights`.
    """
    from sklearn.metrics import roc_auc_score
    from scipy.stats import rankdata, norm
    out = dict.fromkeys(ENSEMBLE_COMBINERS, np.nan)
    dfs = [pd.read_csv(p) for p in pred_paths]
    if len(dfs) < 2:
        return out
    ref = dfs[0]
    y = ref["true_label"].to_numpy()
    if len(np.unique(y)) < 2:
        return out
    P = []
    for d in dfs:
        if not np.array_equal(d["sample_id"].to_numpy(),
                              ref["sample_id"].to_numpy()):
            return out
        P.append(np.clip(d["prob"].to_numpy(), eps, 1 - eps))
    P = np.vstack(P)                                   # (members, samples)
    L = np.log(P / (1 - P))

    n_mem, n_samp = P.shape
    R = np.vstack([rankdata(row) for row in P])

    out["logit"] = float(roc_auc_score(y, L.mean(axis=0)))
    out["prob"] = float(roc_auc_score(y, P.mean(axis=0)))
    out["median"] = float(roc_auc_score(y, np.median(L, axis=0)))
    out["rank"] = float(roc_auc_score(y, R.mean(axis=0)))
    out["rankmed"] = float(roc_auc_score(y, np.median(R, axis=0)))
    # Standardise each member over the test cohort: subtracting its own mean
    # and dividing by its own spread is what removes the calibration it
    # inherited from whichever cohort it validated on.
    sd = L.std(axis=1, keepdims=True)
    Z = (L - L.mean(axis=1, keepdims=True)) / np.clip(sd, 1e-9, None)
    out["zlogit"] = float(roc_auc_score(y, Z.mean(axis=0)))
    # Ranks through the normal quantile, i.e. van der Waerden scores. Division
    # by n + 1 keeps the extremes off +/- infinity.
    out["normal"] = float(roc_auc_score(y, norm.ppf(R / (n_samp + 1)).mean(axis=0)))
    k = int(n_mem * TRIM_FRAC)
    if n_mem - 2 * k >= 1:
        out["trim"] = float(roc_auc_score(
            y, np.sort(L, axis=0)[k:n_mem - k].mean(axis=0)))
    if weights is not None:
        w = np.clip(np.asarray(weights, dtype=float), 0, None)
        if w.shape == (P.shape[0],) and w.sum() > 0:
            w = w / w.sum()
            out["wauc"] = float(roc_auc_score(y, (w[:, None] * L).sum(axis=0)))
    return out


def collect(task, folds, results, args):
    by_fold = {}
    for r in results:
        if r["task"] == task:
            by_fold.setdefault(r["fold"], []).append(r)
    n_members = {f["fold"]: f["n_members"] for f in folds if f["task"] == task}

    rows = []
    for fold, recs in sorted(by_fold.items()):
        recs = sorted(recs, key=lambda r: r["member"])
        paths = [os.path.join(r["out_dir"], "pred_test.csv") for r in recs]
        keep = [i for i, p in enumerate(paths) if os.path.exists(p)]
        preds = [paths[i] for i in keep]
        # Weighting by valid_auc needs those numbers to mean the same thing
        # across members. Under 'joint' they share one validation set and under
        # 'bag' they draw it the same way, so they do; under 'loso' and
        # 'disease_loso' each member was scored on a different set of cohorts
        # and a higher valid_auc may only mean an easier one, so the weighting
        # is withheld and ens_wauc stays NaN.
        wts = ([recs[i]["valid_auc"] - 0.5 for i in keep]
               if args.inner_split in ("bag", "joint") and keep else None)
        combo = ensemble_aucs(preds, weights=wts)

        aucs = [r["test_auc"] for r in recs]
        eps = sorted(r["epoch"] for r in recs)
        strat = sorted({r.get("strategy") or "" for r in recs} - {""})
        row = dict(fold=fold, abund_mode=args.abund_mode,
                   inner_mode=args.inner_split, strategy=",".join(strat),
                   K=len(recs), n_members=n_members.get(fold),
                   bag_frac=args.bag_frac, feat_frac=args.feat_frac,
                   epoch_min=eps[0], epoch_med=int(np.median(eps)),
                   epoch_max=eps[-1],
                   inner_mean=float(np.mean(aucs)),
                   inner_sd=(float(np.std(aucs, ddof=1)) if len(aucs) >= 2
                             else np.nan),
                   **{f"ens_{k}": combo[k] for k in ENSEMBLE_COMBINERS})
        # Every combiner is written to the CSV whatever this is set to; it
        # only decides which one the printed table and ens_gain follow.
        # NaN propagates, so a fold with fewer than two surviving members
        # reports no gain rather than a gain of zero.
        row["ens_gain"] = row[f"ens_{args.report}"] - row["inner_mean"]
        rows.append(row)

    if not rows:
        return None
    df = pd.DataFrame(rows)
    out = TASKS[task]["results"]
    if args.run_name:
        out = out.replace(".csv", f".{args.run_name}.csv")
    df.to_csv(out, index=False)

    print(f"\n=== {task} -> {out} ===")
    head = f"ens_{args.report}"
    cols = ["fold", "inner_mode", "strategy", "K",
            "epoch_min", "epoch_med", "epoch_max",
            "inner_mean", "inner_sd", head, "ens_gain"]
    print(df[cols].to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print(f"  mean: inner_mean {df.inner_mean.mean():.4f}  "
          f"{head} {df[head].mean():.4f}  "
          f"ens_gain {df.ens_gain.mean():+.4f}")

    alt = [f"ens_{k}" for k in ENSEMBLE_COMBINERS
           if df[f"ens_{k}"].notna().any()]
    if len(alt) > 1:
        print("\n  --- combiners, mean over folds ---")
        for c in sorted(alt, key=lambda c: -df[c].mean()):
            n_ok = int(df[c].notna().sum())
            print(f"   {c:<11} {df[c].mean():.4f}"
                  + (f"  ({n_ok}/{len(df)} folds)" if n_ok < len(df) else "")
                  + ("   <- --report-combiner" if c == head else ""))
        print("   All five are recomputed from the same pred_test.csv files, "
              "so the comparison costs nothing. Picking the best of five on "
              "one task is a choice; confirm it on another before reporting "
              "it.")

    short = df[df.K < df.n_members]
    if len(short):
        print(f"\n  warn: {len(short)} fold(s) lost members to failures, so "
              f"their ensemble is over fewer than asked for: "
              + ", ".join(f"{r.fold}(K={r.K}/{r.n_members})"
                          for r in short.itertuples()))
        print("   Look at failed.txt under the fold's members/ directories; "
              "the surviving members are not a random subset if the failures "
              "depended on the draw.")
    return df


# =====================================================================
# 6. Main
# =====================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tasks", nargs="+", default=["disease"],
                    help=f"{list(TASKS)} or 'all'")
    ap.add_argument("--gpus", type=int, nargs="+",
                    default=[0, 1, 2, 3, 4, 5, 6, 7])
    ap.add_argument("--workers-per-gpu", type=int, default=1)

    ap.add_argument("--inner-split", default="bag",
                    choices=["bag", "joint", "loso", "disease_loso",
                             "same_disease"],
                    help="how each member's training and validation parts are "
                         "cut from the fold's training table. 'bag' (default) "
                         "is the random-forest one: every member draws its own "
                         "random --bag-frac and validates on the rest. 'joint' "
                         "cuts one (study,label)-stratified split shared by "
                         "every member, so they differ only in model_seed and "
                         "the validation samples come from the training "
                         "cohorts -- in-distribution, which is not what the "
                         "test cohort measures. 'loso' gives one member per "
                         "training cohort, each validating on the cohort it "
                         "held out, so it selects the epoch that transfers "
                         "across studies; --n-estimators is ignored there "
                         "unless --loso-k is set. 'disease_loso' holds out one "
                         "randomly drawn study per disease (--disease-col), so "
                         "every disease is represented in validation by one "
                         "held-out cohort rather than the epoch being chosen "
                         "on whichever single disease 'loso' happened to hold "
                         "out; each member redraws, so run --n-estimators 8 "
                         "and ensemble. lodo and loso_all only. "
                         "'same_disease' is leave-one-study-out confined to "
                         "the disease under test: one member per study of that "
                         "disease, each validating on it and training on the "
                         "disease's other studies plus every other disease. "
                         "The epoch is then chosen on the same within-disease "
                         "cross-study step the test cohort asks for. "
                         "--n-estimators does not apply; loso_all only")
    ap.add_argument("--cv-folds", type=int, default=5,
                    help="under --inner-split same_disease, how many ways to "
                         "cut a disease that has a single study left in the "
                         "training pool (default 5, i.e. 20%% validation per "
                         "member). Holding that one study out whole would "
                         "leave the model with none of the disease it is "
                         "about to be tested on; cutting it keeps 1-1/k of it "
                         "in training and still validates on the same "
                         "disease. Capped by the smaller class, and a study "
                         "with fewer than two of either is held out whole")
    ap.add_argument("--val-ratio", type=float, default=0.2,
                    help="validation fraction for --inner-split joint")
    ap.add_argument("--split-seed", type=int, default=11,
                    help="seed for the joint split, shared by all members")
    ap.add_argument("--inner-group", default="study",
                    help="metadata column naming the cohort, used by "
                         "--inner-split loso and by the joint stratification")
    ap.add_argument("--loso-k", type=int, default=None,
                    help="under --inner-split loso, hold out this many "
                         "cohorts together as each member's validation set "
                         "instead of one. The fold's cohorts are shuffled "
                         "once (--loso-k-seed) and cut into non-overlapping "
                         "groups of this size (the last group may be "
                         "smaller); each group is then trained "
                         "--n-estimators times, differing only in "
                         "model_seed and sharing that group's split. Omit "
                         "(default) for plain leave-one-cohort-out, which "
                         "does not use --n-estimators")
    ap.add_argument("--loso-k-seed", type=int, default=11,
                    help="seed for the one-time random partition of a "
                         "fold's cohorts into --loso-k-sized groups")
    ap.add_argument("--n-estimators", type=int, default=1,
                    help="members per fold under --inner-split bag, joint or "
                         "disease_loso. Each sees only --bag-frac of the "
                         "samples and --feat-frac of the taxa, so it is weaker "
                         "than a model trained on everything; this is what buys "
                         "that back, and why a forest grows hundreds of trees. "
                         "Under disease_loso each member instead redraws which "
                         "study each disease holds out, and 8 is the count that "
                         "mode is meant to be run at. It does not apply under "
                         "'loso', where the member count is however many "
                         "cohorts the fold's training table holds")
    ap.add_argument("--bag-frac", type=float, default=BAG_FRAC_DEFAULT,
                    help=f"fraction of the training table each member trains "
                         f"on, the rest being its out-of-bag validation set. "
                         f"Default {BAG_FRAC_DEFAULT:.4f} = 1 - 1/e, the "
                         f"fraction of distinct samples a bootstrap of the "
                         f"same size contains. Lower means weaker but more "
                         f"diverse members")
    ap.add_argument("--feat-frac", type=float, default=1.0,
                    help="fraction of the taxa each member sees, counted over "
                         "the taxa non-zero in its own training part; taxa "
                         "zero throughout it are kept whole. 1.0 (default) is "
                         "off and costs no disk; below that a filtered "
                         "training table is written per member")
    ap.add_argument("--bag-seed", type=int, default=101,
                    help="member i draws its bag and its taxa with this plus i")
    ap.add_argument("--ensemble-seed0", type=int, default=11,
                    help="member i initialises its weights with this plus i")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                    dest="set_cfg",
                    help="override one CFG entry, repeatable: --set p_drop=0.3 "
                         f"--set d_ff=400. Keys must already exist in CFG "
                         f"({', '.join(sorted(CFG))}) and the value is read "
                         "with the type the default has, so a typo in either "
                         "is refused rather than silently ignored. Pair it "
                         "with --run-name: a sweep without one writes every "
                         "setting into the same directories")

    ap.add_argument("--abund-mode", default="multiply",
                    choices=["multiply", "none"],
                    help="how abundance is fused with the SNEs vectors. "
                         "'none' ignores it, leaving presence/absence; "
                         "'multiply' scales each OTU embedding by its "
                         "abundance")
    ap.add_argument("--linear-branch", dest="linear_branch",
                    action="store_true", default=True)
    ap.add_argument("--no-linear-branch", dest="linear_branch",
                    action="store_false")
    ap.add_argument("--lin-weight-decay", type=float, default=1e-3)
    ap.add_argument("--glove-embedding", dest="glove_embedding",
                    default=GLOVE_DEFAULT,
                    help="pretrained OTU embedding file; width must equal d_model")
    ap.add_argument("--no-glove", dest="glove_embedding",
                    action="store_const", const=None,
                    help="use fixed random codes instead (still frozen)")
    ap.add_argument("--select-by", default="auc",
                    choices=["loss", "auc", "last"],
                    help="select the checkpoint by validation loss or AUC, or "
                         "'last' to keep the final epoch and train the full "
                         "budget with no early stopping. Use 'last' when the "
                         "validation metric has been shown not to predict test "
                         "performance: an early peak in a noisy metric then "
                         "ends training after a couple of epochs, which is a "
                         "worse outcome than not selecting at all. Pair it "
                         "with --set num_epochs=N, since N now decides how "
                         "long training actually runs")
    ap.add_argument("--logit-adjust-tau", type=float, default=1.0,
                    help="strength of the logit adjustment, used only with "
                         "--set loss=LogitAdjusted. 1.0 (default) removes each "
                         "disease's own case/control base rate exactly (Menon "
                         "et al., ICLR 2021); 0 disables it. The correction is "
                         "applied to the training logits only -- validation "
                         "and test are scored on the raw, prior-free score, "
                         "which is the point under a left-out-disease split "
                         "where the test cohort's base rate is unknown. Note "
                         "it corrects base rates, not the fact that one "
                         "disease has far more samples than another: the "
                         "disease size cancels out of the ratio, and that "
                         "imbalance is what --train-sampler addresses. "
                         "lodo and loso_all only")
    ap.add_argument("--drop-disease", nargs="*", default=[], metavar="NAME",
                    help="remove these diseases (--disease-col) from every "
                         "fold's training AND validation pool, before the "
                         "inner split. Their own test fold still runs and is "
                         "still reported -- only their contribution to other "
                         "folds' training is removed. Use it for a cohort "
                         "whose case/control contrast is not comparable to "
                         "the rest: one that scores well below 0.5 is being "
                         "learnt backwards, and while it sits in the training "
                         "pool it teaches every other fold the opposite "
                         "relationship. Names are checked against the "
                         "metadata at startup, so a typo is refused")
    ap.add_argument("--valid-auc", default="pooled",
                    choices=["pooled", "macro"],
                    help="how the validation AUC that selects the epoch is "
                         "formed. 'pooled' (default) is one AUC over every "
                         "validation sample, which lets the disease with the "
                         "most samples decide the epoch and, worse, counts "
                         "cross-disease case/control pairs -- with K diseases "
                         "only about 1/K of the pooled pairs are the "
                         "within-disease ones anything is tested on. 'macro' "
                         "scores each disease (--disease-col) on its own and "
                         "averages, so every disease weighs the same. Pair it "
                         "with --inner-split disease_loso, whose validation "
                         "set holds one cohort per disease. Only the AUC "
                         "changes; the loss and the decision threshold stay "
                         "pooled")
    ap.add_argument("--min-group-n", type=int, default=0,
                    help="under --valid-auc macro, leave a disease out of the "
                         "average when its validation set has fewer than this "
                         "many samples. Default 0 keeps every disease holding "
                         "both classes. Equal votes equalise the level but not "
                         "the variance: an AUC's noise goes as 1/sqrt(n), so a "
                         "disease with ~24 validation samples swings ~0.13 per "
                         "epoch on noise alone and can end up choosing the "
                         "epoch on its own. The disease is still trained on -- "
                         "this only removes it from the selection signal. Read "
                         "the '[macro] scored:' line of a run to pick a value")
    ap.add_argument("--train-sampler", default="none",
                    choices=["none", "disease_undersample"],
                    help="how the training loader draws its batches. 'none' "
                         "(default) shuffles the whole training part, which is "
                         "what every result so far was produced with. "
                         "'disease_undersample' caps each disease at "
                         "--max-ratio times the smallest one's sample count "
                         "and redraws the subset each epoch, preserving each "
                         "disease's own case/control ratio. Only the training "
                         "loader is affected; validation and test are always "
                         "scored on their full cohorts. Consider raising "
                         "--patience to 20-25 alongside it: the changing "
                         "subset makes the validation curve noisier")
    ap.add_argument("--disease-col", default="disease_name_ab",
                    help="metadata column naming each sample's disease. Read "
                         "under --train-sampler disease_undersample, which "
                         "caps each disease's share of a training epoch, and "
                         "under --inner-split disease_loso, which holds out "
                         "one study per disease. On a single-disease task the "
                         "undersampling is a no-op")
    ap.add_argument("--max-ratio", type=float, default=10.0,
                    help="largest allowed ratio between the biggest and the "
                         "smallest disease after capping (default 10, one "
                         "order of magnitude)")
    ap.add_argument("--sample-seed", type=int, default=None,
                    help="seed for the per-epoch undersampling draws. Default "
                         "None uses each member's model_seed, so an ensemble "
                         "gets diverse subsets as well as diverse weights. Pin "
                         "it to hold every member's subsets identical")
    ap.add_argument("--patience", type=int, default=15,
                    help="stop this many epochs after the best validation "
                         f"metric. Capped by num_epochs={CFG['num_epochs']}")
    ap.add_argument("--head-hidden", type=int, default=None,
                    help=f"hidden width of the output head, one integer "
                         f"(default {CFG['head_hidden']}). Pass 0 for a bare "
                         f"linear head; omit for the library default of "
                         f"d_model // 2")

    ap.add_argument("--overwrite", dest="overwrite", action="store_true",
                    default=True,
                    help="wipe and rerun a job whose output directory already "
                         "exists (default)")
    ap.add_argument("--resume", dest="overwrite", action="store_false",
                    help="keep finished jobs and run only the missing ones")
    ap.add_argument("--report-combiner", dest="report", default="prob",
                    choices=list(ENSEMBLE_COMBINERS),
                    help="which combiner the printed table and ens_gain "
                         "follow. Every combiner is computed and written to "
                         "the CSV whatever this is set to, so it changes the "
                         "headline number and nothing else. 'zlogit' or "
                         "'normal' are the ones to reach for when the members "
                         "were calibrated on different cohorts")
    ap.add_argument("--run-name", default="")
    ap.add_argument("--keep-ckpt", action="store_true")
    ap.add_argument("--keep-tables", action="store_true",
                    help="keep each member's feature-subsampled training table "
                         "instead of deleting it on success")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    tasks = list(TASKS) if "all" in args.tasks else args.tasks
    for t in tasks:
        if t not in TASKS:
            raise SystemExit(f"unknown task {t}; choose from {list(TASKS)} or 'all'")
    if args.n_estimators < 1:
        raise SystemExit(f"--n-estimators must be at least 1, got "
                         f"{args.n_estimators}")
    if args.inner_split != "bag" and args.bag_frac != BAG_FRAC_DEFAULT:
        print(f"[warn] --bag-frac only applies to --inner-split bag; ignored "
              f"under {args.inner_split!r}")
    if args.train_sampler == "none" and args.max_ratio != 10.0:
        print("[warn] --max-ratio only applies under --train-sampler "
              "disease_undersample; ignored")
    if args.train_sampler != "none" and args.max_ratio <= 0:
        raise SystemExit(f"--max-ratio must be positive, got {args.max_ratio}")
    if args.loso_k is not None and args.loso_k < 1:
        raise SystemExit(f"--loso-k must be at least 1, got {args.loso_k}")
    if args.loso_k is not None and args.inner_split != "loso":
        print(f"[warn] --loso-k only applies to --inner-split loso; ignored "
              f"under {args.inner_split!r}")
    if args.valid_auc == "pooled" and args.min_group_n:
        print("[warn] --min-group-n only applies under --valid-auc macro; "
              "ignored")
    if args.min_group_n < 0:
        raise SystemExit(f"--min-group-n must be >= 0, got {args.min_group_n}")
    if args.valid_auc == "macro" and args.select_by == "loss":
        # The epoch is chosen by whichever quantity --select-by names, and
        # under 'loss' that is the pooled BCE -- a per-sample average, so a
        # disease with ten times the samples carries ten times the weight,
        # which is the imbalance --valid-auc macro exists to remove. The macro
        # AUC would still be computed and plotted, and still decide nothing.
        raise SystemExit(
            "--valid-auc macro with --select-by loss selects the epoch on the "
            "pooled validation loss, so the macro AUC is computed, plotted, "
            "and then ignored -- and the pooled loss weights each disease by "
            "its sample count, which is what macro is for. Use --select-by "
            "auc (the default) alongside --valid-auc macro, or drop the macro")
    if args.valid_auc == "macro":
        single = [t for t in tasks if t not in ("lodo", "loso_all")]
        if single:
            # One disease means one group, so the macro average is that single
            # disease's pooled AUC -- the same number, reached the long way.
            print(f"[warn] --valid-auc macro on {single}: those tasks are "
                  f"single-disease, so the average is over one group and "
                  f"equals the pooled AUC")
        if args.inner_split in ("bag", "joint"):
            print(f"[warn] --valid-auc macro under --inner-split "
                  f"{args.inner_split}: validation samples come from the same "
                  f"cohorts as training, so a per-disease AUC there measures "
                  f"within-cohort separation. --inner-split disease_loso is "
                  f"the split this pairs with")
    _la = parse_set(args.set_cfg).get("loss", CFG["loss"]) == "LogitAdjusted"
    if _la:
        bad = [t for t in tasks if t not in ("lodo", "loso_all")]
        if bad:
            raise SystemExit(
                f"loss=LogitAdjusted needs a training table spanning several "
                f"diseases, and {bad} is single-disease (its metadata is one "
                f"file per disease), where the per-disease base rate collapses "
                f"to one global number. Use --tasks lodo and/or loso_all")
        if args.logit_adjust_tau < 0:
            raise SystemExit(f"--logit-adjust-tau must be >= 0, got "
                             f"{args.logit_adjust_tau}")
        print(f"[cfg] loss=LogitAdjusted tau={args.logit_adjust_tau:g}: each "
              f"disease's case/control base rate is removed from the training "
              f"logits (column {args.disease_col!r}); validation and test see "
              f"raw logits")
    elif args.logit_adjust_tau != 1.0:
        print("[warn] --logit-adjust-tau only applies with "
              "--set loss=LogitAdjusted; ignored")
    if args.inner_split == "same_disease":
        bad = [t for t in tasks if t != "loso_all"]
        if bad:
            raise SystemExit(
                f"--inner-split same_disease validates on the other studies "
                f"of the test study's disease, and {bad} does not hold one "
                f"study out per fold. Use --tasks loso_all")
        if args.inner_group == args.disease_col:
            raise SystemExit(
                f"--inner-group and --disease-col are both "
                f"{args.inner_group!r}; under same_disease they must name the "
                f"cohort column and the disease column separately")
        print("[info] --n-estimators does not apply under --inner-split "
              "same_disease: a fold runs one member per study of the disease "
              "under test, so folds whose disease has more cohorts get more "
              "members, and their predictions are averaged")
    if args.inner_split == "disease_loso" and args.inner_group == args.disease_col:
        # Both columns pointing at the same thing makes every disease look like
        # it has exactly one cohort -- itself -- so the draw holds out the whole
        # disease, every disease, and the training part comes out empty. It
        # fails one job at a time deep inside the pool, so catch it here.
        raise SystemExit(
            f"--inner-group and --disease-col are both {args.inner_group!r}. "
            f"Under --inner-split disease_loso they name different things: "
            f"--inner-group is the cohort/study column (default 'study') and "
            f"--disease-col is the disease column (default 'disease_name_ab'). "
            f"Pointing both at one column makes each disease its own only "
            f"cohort, so every disease is held out whole and nothing is left "
            f"to train on. Drop --inner-group to take the default")
    if args.inner_split == "disease_loso" and args.n_estimators == 1:
        print("[warn] --inner-split disease_loso with --n-estimators 1 keeps "
              "one random study per disease and never redraws, so the epoch "
              "rides on a single draw. The mode is meant to be ensembled: "
              "pass --n-estimators 8")
    if args.inner_split == "loso" and args.loso_k is None:
        print("[info] --n-estimators does not apply under --inner-split loso: "
              "each fold runs one member per cohort subset its own training "
              "table supports, so folds with more cohorts get more members")
    elif args.inner_split == "loso":
        print(f"[info] --inner-split loso --loso-k {args.loso_k}: a fold's "
              f"cohorts are shuffled once (--loso-k-seed {args.loso_k_seed}) "
              f"and cut into non-overlapping groups of {args.loso_k}, each "
              f"trained --n-estimators={args.n_estimators} times, differing "
              f"only in model_seed")
    if not 0.05 <= args.bag_frac <= 0.95:
        raise SystemExit(
            f"--bag-frac {args.bag_frac} leaves too little on one side; the "
            f"complement is the member's validation set, so keep it within "
            f"[0.05, 0.95]")
    if not 0 < args.feat_frac <= 1:
        raise SystemExit(f"--feat-frac must be in (0, 1]; got {args.feat_frac}")
    if args.n_estimators == 1 and args.inner_split == "bag":
        print("[warn] --n-estimators 1 trains one model on a random part of "
              "the data, which is strictly worse than training it on all of "
              "it. The bagging only pays off across members")
    if args.run_name:
        for t in tasks:
            TASKS[t]["artifact"] = f"{TASKS[t]['artifact']}{args.run_name}"

    print(f"[cfg] tasks={tasks} gpus={args.gpus} x{args.workers_per_gpu}")
    if args.inner_split == "bag":
        print(f"[cfg] inner-split=bag: {args.n_estimators} members per fold, "
              f"each trained on a plain random {args.bag_frac:.1%} of the "
              f"training table and validated on the {1 - args.bag_frac:.1%} it "
              f"did not draw. Seeds bag_seed={args.bag_seed}.."
              f"{args.bag_seed + args.n_estimators - 1}")
    elif args.inner_split == "joint":
        print(f"[cfg] inner-split=joint: one (study,label)-stratified split at "
              f"{args.val_ratio:.0%} validation, shared by "
              f"{args.n_estimators} members that differ only in model_seed. "
              f"The validation samples come from the training cohorts, so this "
              f"selects for generalization within the training distribution")
    elif args.inner_split == "same_disease":
        print(f"[cfg] inner-split=same_disease: one member per study (column "
              f"{args.inner_group!r}) of the test study's disease (column "
              f"{args.disease_col!r}); each validates on that one study and "
              f"trains on the disease's remaining studies plus every other "
              f"disease. The member count is the fold's, not a flag's")
        print(f"[cfg] a disease down to one study is cut --cv-folds="
              f"{args.cv_folds} ways instead, so it stays in training")
    elif args.inner_split == "disease_loso":
        print(f"[cfg] inner-split=disease_loso: {args.n_estimators} members "
              f"per fold, each holding out one randomly drawn study (column "
              f"{args.inner_group!r}) per disease (column {args.disease_col!r}) "
              f"as its validation set. Seeds bag_seed={args.bag_seed}.."
              f"{args.bag_seed + args.n_estimators - 1}")
        print(f"[cfg] the unit held out is the (disease, study) cell, so a "
              f"study spanning several diseases keeps its other diseases in "
              f"training and does sit on both sides of that member's split")
    elif args.loso_k is None:
        print(f"[cfg] inner-split=loso: one member per training cohort "
              f"(column {args.inner_group!r}), each validating on the cohort "
              f"it held out. The member count is the fold's, not a flag's")
    else:
        print(f"[cfg] inner-split=loso: cohorts (column {args.inner_group!r}) "
              f"grouped --loso-k={args.loso_k} at a time (seed "
              f"{args.loso_k_seed}), each group trained "
              f"--n-estimators={args.n_estimators} times differing only in "
              f"model_seed")
    if args.inner_split != "loso":
        print(f"[cfg] model_seed={args.ensemble_seed0}.."
              f"{args.ensemble_seed0 + args.n_estimators - 1}")
    if args.feat_frac < 1.0:
        print(f"[cfg] feat_frac={args.feat_frac}: each member sees that "
              f"fraction of the taxa non-zero in its own training part "
              f"(taxa zero throughout it are kept whole, so the vocabulary "
              f"still covers what only the test cohort carries), and writes "
              f"its own training table")
    print(f"[cfg] abund_mode={args.abund_mode} "
          f"linear_branch={args.linear_branch} select_by={args.select_by} "
          f"patience={args.patience} train_sampler={args.train_sampler} "
          f"valid_auc={args.valid_auc}"
          + (f" min_group_n={args.min_group_n}"
             if args.valid_auc == "macro" else ""))
    if args.train_sampler != "none":
        seed_note = (args.sample_seed if args.sample_seed is not None
                     else "model_seed")
        print(f"[cfg] disease undersampling: column={args.disease_col!r} "
              f"max_ratio={args.max_ratio} sample_seed={seed_note}")
    # Parsed here as well as in build_jobs so a typo costs a second, not a
    # GPU-hour, and so the deviation from CFG is on the record above the run.
    cfg_over = parse_set(args.set_cfg)
    if cfg_over:
        print("[cfg] --set overrides: " + ", ".join(
            f"{k}={v!r} (default {CFG[k]!r})" for k, v in cfg_over.items()))
        if not args.run_name:
            print("[cfg] WARNING: --set without --run-name writes this "
                  "setting into the same directories as every other setting")
    if args.drop_disease:
        print(f"[cfg] dropping {', '.join(args.drop_disease)} from every "
              f"fold's training and validation pool; their own test folds "
              f"still run")
    print(f"[cfg] glove={args.glove_embedding or 'None (random codes)'}")
    if not args.run_name:
        print("[cfg] no --run-name: job directories and results are shared "
              "across settings, so a sweep overwrites itself")
    print(f"[cfg] existing job directories: "
          f"{'wiped and rerun' if args.overwrite else 'kept, jobs skipped'}")

    all_jobs, all_folds = [], []
    for t in tasks:
        j, f = build_jobs(t, args)
        all_jobs += j
        all_folds += f
    todo = all_jobs if args.overwrite else [
        j for j in all_jobs
        if not os.path.exists(os.path.join(j["out_dir"], "result.json"))]
    per_fold = [f["n_members"] for f in all_folds]
    if per_fold and min(per_fold) == max(per_fold):
        shape = f"{len(all_folds)} folds x {per_fold[0]} members"
    else:
        shape = (f"{len(all_folds)} folds, {min(per_fold)}..{max(per_fold)} "
                 f"members each")
    print(f"\n{len(all_jobs)} jobs ({shape}), {len(todo)} to run")
    if args.inner_split == "loso" and per_fold:
        print(f"[cfg] model_seed={args.ensemble_seed0}.."
              f"{args.ensemble_seed0 + max(per_fold) - 1}, member i of every "
              f"fold using the same one")

    if args.dry_run:
        for t in tasks:
            n = sum(1 for j in all_jobs if j["task"] == t)
            nf = sum(1 for f in all_folds if f["task"] == t)
            print(f"  {t:12s}: {nf} folds, {n} jobs")
        slots = len(args.gpus) * args.workers_per_gpu
        print(f"  at 5 min per job, this takes about "
              f"{len(todo) / slots * 5 / 60:.1f} h on {slots} slots")
        return

    run_pool(todo, args.gpus, args.workers_per_gpu, "members")

    if not args.keep_tables:
        cleanup_shared_splits(all_jobs)

    results = []
    for j in all_jobs:
        p = os.path.join(j["out_dir"], "result.json")
        if os.path.exists(p):
            with open(p) as f:
                results.append(json.load(f))
    print(f"\nDone: {len(results)}/{len(all_jobs)} members")

    for t in tasks:
        collect(t, all_folds, results, args)


if __name__ == "__main__":
    main()
