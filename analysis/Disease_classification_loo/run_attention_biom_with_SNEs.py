#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Ensemble of Attention_biom across all LOSO tasks, on multiple GPUs.

One protocol, cohort-held-out validation
----------------------------------------
Each fold trains several members, and every member's validation part is a whole
held-out cohort rather than a random slice of the training table: see
--inner-split for the modes. Selecting the epoch on a cohort the member never
trained on is the same step the test table asks for, which is why the random
in-distribution splits this script used to offer are gone. What separates a
member from its siblings is therefore which cohorts it holds out, and nothing
else: every member sees every taxon of whatever it does train on.

The split is made here rather than in the library, which takes the training and
validation cohorts as two separate BIOM tables and no longer partitions
anything itself. Each member therefore writes its own pair of tables into its
job directory, and they are deleted once it succeeds unless --keep-tables.

The fold's result is the ensemble of its members' test predictions, not any
single member. A member that gave a whole cohort up to validation is weaker
than one trained on everything; the member count is what buys that back.

Combining the members
---------------------
AUC reads only the ordering of the scores, and the members are not on a common
scale -- each one's probabilities come from a threshold fitted on the cohort it
held out. Several combiners are therefore scored side by side, all recomputed
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
  ens_wauc    logits weighted by valid_auc - 0.5. Always NaN now: every
              remaining --inner-split scores each member on a different set of
              cohorts, so a higher valid_auc may only mean an easier one and
              the weights would not be comparable.

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
  python run_attention_biom_with_SNEs.py --tasks disease \
      --gpus 0 1 2 3 4 5 6 7 --run-name rf_loso
  python run_attention_biom_with_SNEs.py --tasks lodo \
      --inner-split disease_loso --n-estimators 8 --run-name rf_dloso
  python run_attention_biom_with_SNEs.py --tasks lodo loso_all \
      --inner-split per_disease --run-name rf_pdis
  python run_attention_biom_with_SNEs.py --tasks lodo loso_all \
      --set loss=GroupBalanced --group-balance-beta 0.5 \
      --valid-auc macro --run-name rf_gb05
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import re
import shutil
import time
import traceback

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
    "shuffled_table": dict(
        run_tsv="run_shuffled_table.tsv",
        train="Data/shuffle_table_IBD_CRC/{disease}/{study}/train_loo.biom",
        test="Data/shuffle_table_IBD_CRC/{disease}/{study}/test_loo.biom",
        meta="Data/shuffle_table_IBD_CRC/{disease}/metadata.tsv",
        fold="{disease}_{study}",
        artifact="Data/shuffle_table_IBD_CRC/",
        results="disease.csv")
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


# =====================================================================
# 1. Materialising a member's split as two BIOM tables
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


def drop_study_mask(sids, md, study_col, drop_studies):
    """Boolean over `sids`; True keeps the sample, False drops its cohort.

    ``--drop-studies`` exists for a cohort whose labels cannot be trusted --
    self-reported diagnoses, say. Such a cohort has to leave every table at
    once: keeping it in training teaches the model those labels, and keeping it
    in test scores the model against them. Filtering in one place, applied at
    each of the three, is what keeps the two halves from drifting apart.

    Blanks are kept rather than refused, unlike :func:`_column`: a sample whose
    study is unrecorded is not a member of a named dropped cohort, and this is
    not the function that should be deciding the metadata is malformed.
    """
    if not drop_studies:
        return np.ones(len(sids), dtype=bool)
    g = md.loc[sids, study_col].astype(str).to_numpy()
    return ~np.isin(g, [str(s) for s in drop_studies])


def survey_training_table(train_biom, meta_path, study_col, drop_studies=()):
    """The cohorts present in a fold's training table."""
    sids = biom.load_table(train_biom).ids(axis="sample")
    md = pd.read_csv(meta_path, sep="\t", index_col=SAMPLE_ID_COL,
                     dtype={SAMPLE_ID_COL: str}, low_memory=False)
    sids = np.asarray(sids)[drop_study_mask(np.asarray(sids), md, study_col,
                                            drop_studies)]
    return md.loc[sids, study_col].unique()


def _valid_mask(job, sids, lab, md):
    """Boolean mask over `sids` (True = validation), plus a strategy label."""
    mode = job["inner_mode"]

    if mode == "same_disease":
        return _same_disease_mask(job, sids, lab, md)

    if mode == "per_disease":
        # One member per disease, validating on one drawn study of that disease
        # alone. Where 'disease_loso' gives one member a validation set spanning
        # every disease, this spreads the same draw across members: the fold
        # ends up covering every disease too, but through the ensemble rather
        # than inside any one member's split.
        return _per_disease_mask(job, sids, md)

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

    The draw uses ``job['member_seed']``, which is ``--member-seed`` plus the
    member index, so the members differ in which cohorts they hold out and
    their predictions are worth ensembling.
    """
    rng = np.random.default_rng(job["member_seed"])
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
        # -- the members ensemble, so the fold gets back across members what
        # each one gave up -- but a member training on under half its fold is
        # worth knowing about.
        print(f"[dloso] {tag}: WARNING: that leaves only {1 - share:.0%} of the "
              f"fold to train on; the diseases here have few cohorts each")
    if lonely:
        print(f"[dloso] {tag}: * = {', '.join(lonely)} had one study only, so "
              f"that disease is now absent from this member's training part")
    return mask, "disease_loso"


def _per_disease_mask(job, sids, md):
    """Hold out one study of *this member's* disease, and nothing else.

    The member is the disease: a fold whose training table spans twelve
    diseases trains twelve members, member *d* validating on the single study
    drawn for *d* in :func:`survey_per_disease_units`. Every other disease, and
    every other study of *d*, stays in training.

    Compared with ``disease_loso``, which hands one member a validation set
    holding one cohort of every disease at once, this cuts the same draw up
    across members. Each member therefore gives up far less of its fold -- one
    cohort of one disease rather than one cohort of each -- and the epoch it
    selects is the one that transfers across studies *for its own disease*.
    What makes the fold's answer speak for every disease is the ensemble over
    members, not any single split.

    As in ``disease_loso`` the unit held out is the ``(disease, study)`` cell,
    not the whole study, so a study spanning several diseases -- PRJEB11419
    carries OB, IBS, BD and T2DM -- contributes only its samples of this
    member's disease and keeps the rest in training. That leaves one study on
    both sides of the split, so this member's valid_auc is not a clean
    cross-study estimate for a disease that shares its cohort.

    The draw itself happened in the parent process, so this only applies it:
    ``job['valid_disease']`` names the disease and ``job['valid_studies']``
    holds the one study drawn for it.
    """
    dis = _column(md, sids, job["disease_col"], "disease")
    stu = _column(md, sids, job["study_col"], "study")
    d = str(job["valid_disease"])
    chosen = str(job["valid_studies"][0])

    in_d = dis == d
    mask = in_d & (stu == chosen)
    tag = f"{job['fold']}/{job['valid_unit']}"
    if not mask.any():
        raise ValueError(
            f"{tag}: no training sample has {job['disease_col']}=={d!r} and "
            f"{job['study_col']}=={chosen!r}, so validation would be empty")

    left = int((in_d & ~mask).sum())
    print(f"[pdis] {tag}: validating on {d}/{chosen}, {int(mask.sum())}/"
          f"{len(sids)} samples ({mask.sum() / len(sids):.1%}); {d} keeps "
          f"{left} training sample(s), the other "
          f"{len(set(dis)) - 1} disease(s) all of theirs")
    if left == 0:
        print(f"[pdis] {tag}: {d} had this one study, so it is absent from "
              f"this member's training part")
        return mask, "per_disease_only_study"
    return mask, "per_disease"


def survey_per_disease_units(train_biom, meta_path, study_col, disease_col,
                             seed, min_class_n=1, drop_studies=()):
    """Draw one study per disease from a fold's training table.

    Done in the parent process rather than in the worker, because the member
    list *is* the disease list: `build_jobs` has to know which diseases survive
    the draw before it can queue anything. Doing it here also puts the whole
    plan on the record above the run.

    A disease's candidates are the studies whose ``(disease, study)`` cell holds
    at least `min_class_n` cases **and** at least `min_class_n` controls. A cell
    short on either cannot carry the epoch decision: with none of a class the
    AUC is undefined outright, and with a handful the AUC is mostly sampling
    noise, so the member would early-stop on a coin flip. Candidates are
    filtered before the draw rather than after, which is the difference between
    never picking such a cell and failing the member that picked it. A disease
    left with no candidate at all gets no member, and says so.

    The threshold applies to each class separately on purpose. A cohort of 200
    controls and 4 cases has plenty of samples and still cannot say whether an
    epoch separates the classes, so a rule on the cell's total would let it
    through.

    Each disease is drawn for independently, in sorted order so the sequence
    depends on `seed` alone and not on the order the diseases happen to appear
    in the table.

    Returns
    -------
    picks : dict
        ``{disease: {'study', 'n_eligible', 'n_studies', 'n_valid'}}``, in
        sorted disease order. One entry per member the fold will run.
    skipped : dict
        ``{disease: reason}`` for the diseases that had no two-class study.
    """
    sids = np.asarray(biom.load_table(train_biom).ids(axis="sample"))
    md = pd.read_csv(meta_path, sep="\t", index_col=SAMPLE_ID_COL,
                     dtype={SAMPLE_ID_COL: str}, low_memory=False)
    # Before anything is counted: a dropped cohort must not be a candidate, and
    # must not count towards a disease's cohort tally either.
    keep = drop_study_mask(sids, md, study_col, drop_studies)
    # A disease living only in a dropped cohort leaves the table altogether, so
    # the loop below never reaches it and it would vanish without a word. Note
    # it here, while both sides of the filter are still in hand.
    vanished = (set(md.loc[sids, disease_col].astype(str))
                - set(md.loc[sids[keep], disease_col].astype(str)))
    sids = sids[keep]
    dis = _column(md, sids, disease_col, "disease")
    stu = _column(md, sids, study_col, "study")
    lab = pd.to_numeric(md.loc[sids, LABELS_COL], errors="coerce").to_numpy()

    rng = np.random.default_rng(seed)
    picks, skipped = {}, {}
    for d in sorted(vanished):
        skipped[d] = ("--drop-studies removed every cohort it had, so the "
                      "disease is no longer in this fold's training table at "
                      "all")
    for d in sorted(set(dis)):
        in_d = dis == d
        studies = sorted(set(stu[in_d]))
        eligible, sizes = [], []
        for s in studies:
            cell = in_d & (stu == s)
            n_pos = int((lab[cell] == 1).sum())
            n_neg = int((lab[cell] == 0).sum())
            sizes.append(f"{s} {n_pos}/{n_neg}")
            if n_pos >= min_class_n and n_neg >= min_class_n:
                eligible.append(s)
        if not eligible:
            skipped[d] = (
                f"none of its {len(studies)} study/studies reaches "
                f"{min_class_n} case(s) and {min_class_n} control(s) "
                f"(case/control by study: {', '.join(sizes)})")
            continue
        chosen = eligible[int(rng.integers(len(eligible)))]
        cell = in_d & (stu == chosen)
        picks[d] = dict(study=chosen, n_eligible=len(eligible),
                        n_studies=len(studies),
                        n_valid=int(cell.sum()),
                        n_case=int((lab[cell] == 1).sum()),
                        n_ctrl=int((lab[cell] == 0).sum()))
    return picks, skipped


def survey_same_disease_studies(train_biom, meta_path, study_col, disease_col,
                                test_study, drop_studies=()):
    """Training-table studies carrying the disease(s) of the held-out study."""
    sids = np.asarray(biom.load_table(train_biom).ids(axis="sample"))
    md = pd.read_csv(meta_path, sep="\t", index_col=SAMPLE_ID_COL,
                     dtype={SAMPLE_ID_COL: str}, low_memory=False)
    sids = sids[drop_study_mask(sids, md, study_col, drop_studies)]
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


def test_path(out_dir):
    return os.path.join(out_dir, "test_part.biom")


def _filtered_test_table(job, md):
    """The fold's test table with `--drop-studies` removed, as a path.

    Returns ``job['test']`` untouched when nothing is dropped, so a run without
    the flag reads the original file and writes nothing extra. Otherwise a
    filtered copy is written beside the member's split; `build_jobs` has
    already refused any fold this would empty.
    """
    if not job.get("drop_studies"):
        return job["test"]
    out = test_path(job["split_dir"])
    if os.path.exists(out):
        return out
    table = biom.load_table(job["test"])
    sids = np.asarray(table.ids(axis="sample"))
    keep = drop_study_mask(sids, md, job["study_col"], job["drop_studies"])
    if not keep.any():
        raise ValueError(
            f"{job['fold']}: every test sample belongs to --drop-studies "
            f"{list(job['drop_studies'])}, so there is nothing to score")
    _write_table(table.filter(sids[keep], axis="sample", inplace=False), out)
    print(f"[drop ] {job['fold']}/{job['valid_unit']}: test table "
          f"{len(sids)} -> {int(keep.sum())} samples")
    return out


def member_split(job):
    """Materialise this member's training and validation parts as BIOM tables.

    The library no longer splits the training table itself, so the partition is
    made here and handed over as files. Which partition depends on
    --inner-split; see `_valid_mask`.

    ``--drop-studies`` is applied first, to all three tables. It has to reach
    the test table as well as the training one: a cohort excluded because its
    labels are doubtful would otherwise still be scored against.

    Returns ``(train_path, valid_path, test_path, strategy)``. Existing files
    are reused, so a retried job -- or, under grouped 'loso', every repeat of a
    group after the first -- does not pay for the filtering again.
    """
    split_dir = job["split_dir"]
    train_part, valid_part = split_paths(split_dir)
    info_path = os.path.join(split_dir, "split.json")
    md = pd.read_csv(job["meta"], sep="\t", index_col=SAMPLE_ID_COL,
                     dtype={SAMPLE_ID_COL: str}, low_memory=False)
    if all(os.path.exists(x) for x in (train_part, valid_part, info_path)):
        with open(info_path) as f:
            return (train_part, valid_part, _filtered_test_table(job, md),
                    json.load(f)["strategy"])
    os.makedirs(split_dir, exist_ok=True)

    table = biom.load_table(job["train"])
    sids = np.asarray(table.ids(axis="sample"))
    missing = [x for x in sids if x not in md.index]
    if missing:
        raise KeyError(f"{len(missing)} training samples have no metadata row; "
                       f"first offenders: {missing[:5]}")

    # Before the split, so a dropped cohort can reach neither part. The surveys
    # in build_jobs filtered the same way, so the member plan already assumed
    # this table.
    keep = drop_study_mask(sids, md, job["study_col"], job.get("drop_studies"))
    if not keep.all():
        n_before = len(sids)
        sids = sids[keep]
        table = table.filter(sids, axis="sample", inplace=False)
        if len(sids) == 0:
            raise ValueError(
                f"{job['fold']}/{job['valid_unit']}: every training sample "
                f"belongs to --drop-studies {list(job['drop_studies'])}")
        print(f"[drop ] {job['fold']}/{job['valid_unit']}: training table "
              f"{n_before} -> {len(sids)} samples")

    lab = pd.to_numeric(md.loc[sids, LABELS_COL], errors="coerce")
    if lab.isna().any():
        bad = lab.index[lab.isna()][:5].tolist()
        raise ValueError(f"non-numeric labels in {LABELS_COL!r}: {bad}")
    lab = lab.to_numpy().astype(int)

    is_v, strategy = _valid_mask(job, sids, lab, md)

    # Nothing in a leave-one-cohort-out draw guarantees both classes land on
    # both sides. A part that lost a class cannot be trained on or
    # early-stopped against, so the member fails here rather than training on a
    # meaningless validation curve.
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
    _write_table(train_tbl, train_part)
    _write_table(table.filter(sids[is_v], axis="sample", inplace=False),
                 valid_part)
    with open(info_path, "w") as f:
        json.dump(dict(strategy=strategy, n_train=int((~is_v).sum()),
                       n_valid=int(is_v.sum()),
                       dropped_studies=list(job.get("drop_studies") or [])), f)
    return train_part, valid_part, _filtered_test_table(job, md), strategy

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


def _build_kwargs(job, train_biom, valid_biom, test_biom=None):
    kw = dict(
        metadata=job["meta"], train_biom=train_biom, valid_biom=valid_biom,
        test_biom=job["test"] if test_biom is None else test_biom,
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
        valid_auc=job["valid_auc"],
        min_group_n=job["min_group_n"],
        logit_adjust_tau=job["logit_adjust_tau"],
        group_balance_beta=job["group_balance_beta"],
        group_balance_max_ratio=job["group_balance_max_ratio"],
        disease_col=job["disease_col"],
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

        train_part, valid_part, test_part, strategy = member_split(job)
        rec, tm = Attention_biom(
            **_build_kwargs(job, train_part, valid_part, test_part))
        ep = rec["epoch"]
        if not isinstance(ep, int):
            raise RuntimeError(f"unexpected epoch value {ep!r}")
        out = dict(task=job["task"], fold=job["fold"],
                   member=job["member"], valid_unit=str(job["valid_unit"]),
                   inner_mode=job["inner_mode"], strategy=strategy,
                   epoch=int(ep),
                   valid_auc=float(rec["valid_auc"]),
                   # The train/valid pair at the selected epoch. train_auc is
                   # the one quantity that separates a member which fit its
                   # training part too well from one that never fit it at all,
                   # and until now it reached stdout and nowhere else -- which
                   # is unusable, since eight workers interleave their output.
                   train_auc=float(rec["train_auc"]),
                   train_loss=float(rec["train_loss"]),
                   valid_loss=float(rec["valid_loss"]),
                   n_valid=int(len(rec["valid_label"])),
                   test_auc=float(tm["auc"]),
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
            # -- grouped 'loso', where a group's repeats read the same pair --
            # is left alone here: the first member to finish would be deleting
            # a table its siblings still need. The parent drops it after the
            # pool.
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
def _test_survives_drop(test_biom, meta_path, study_col, drop, task, fold):
    """Whether a fold still has a test table once --drop-studies is applied.

    Under loso_all the fold *is* a study, so dropping that study empties its
    test table outright: there is no cohort left to score and the fold cannot
    mean anything. Caught here, in the parent, so the fold never reaches a GPU
    -- and reported, since a silently shorter results table is worse than a
    slightly noisier log.
    """
    sids = np.asarray(biom.load_table(test_biom).ids(axis="sample"))
    md = pd.read_csv(meta_path, sep="\t", index_col=SAMPLE_ID_COL,
                     dtype={SAMPLE_ID_COL: str}, low_memory=False)
    keep = drop_study_mask(sids, md, study_col, drop)
    if keep.all():
        return True
    if not keep.any():
        print(f"[warn] {task}/{fold}: the test table is empty after "
              f"--drop-studies {' '.join(drop)}, so this fold is skipped "
              f"entirely -- it will not appear in the results")
        return False
    print(f"[drop ] {task}/{fold}: test table {len(sids)} -> "
          f"{int(keep.sum())} samples")
    return True


def _unit_name(disease):
    """A disease name made safe to use as a job directory name.

    Under per_disease the member *is* a disease, so the disease is what names
    its directory. Disease labels come from the metadata and are free text, and
    a '/' in one would silently write the member a directory deeper.
    """
    return re.sub(r"[^0-9A-Za-z._-]+", "_", str(disease)).strip("_") or "unnamed"


def build_jobs(task, args):
    """One job per (fold, member)."""
    t = TASKS[task]
    if (args.inner_split in ("disease_loso", "per_disease")
            and task not in ("lodo", "loso_all")):
        raise SystemExit(
            f"--inner-split {args.inner_split} needs a training table spanning "
            f"several diseases, and {task!r} is single-disease (its metadata is "
            f"one file per disease), where the mode would degrade to drawing "
            f"one random study. Use --tasks lodo and/or loso_all, or pick "
            f"--inner-split loso for {task!r}")
    run = pd.read_csv(t["run_tsv"], sep="\t")
    jobs, folds = [], []
    drop = [str(x) for x in (args.drop_studies or [])]

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

        if drop and not _test_survives_drop(test, meta, args.inner_group,
                                            drop, task, fold):
            continue

        # loso only: unit -> studies it holds out, and unit -> a shared split
        # directory when --loso-k groups several repeats onto the same split.
        group_studies, group_split_dir, cv_folds = {}, {}, {}
        # per_disease only: unit -> the disease that unit is the member for.
        unit_disease = {}

        if args.inner_split == "loso":
            cohorts = sorted(str(x) for x in
                             survey_training_table(train, meta, args.inner_group,
                                                   drop_studies=drop))
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
        elif args.inner_split == "per_disease":
            picks, skipped = survey_per_disease_units(
                train, meta, args.inner_group, args.disease_col,
                args.member_seed, args.per_disease_min_class,
                drop_studies=drop)
            for d, why in sorted(skipped.items()):
                print(f"[warn] {task}/{fold}: no member for {d}: {why}, so "
                      f"there is nothing to select an epoch on")
            if not picks:
                raise SystemExit(
                    f"{task}/{fold}: no disease in the training table has a "
                    f"study with {args.per_disease_min_class} of each class, "
                    f"so --inner-split per_disease can build no member at all. "
                    f"Lower --per-disease-min-class (currently "
                    f"{args.per_disease_min_class}) or use --inner-split loso")
            units = []
            for d, info in picks.items():
                u = _unit_name(d)
                units.append(u)
                group_studies[u] = [info["study"]]
                unit_disease[u] = d
            if len(set(units)) != len(units):
                raise SystemExit(
                    f"{task}/{fold}: two diseases in {args.disease_col!r} "
                    f"reduce to the same job-directory name; rename them or "
                    f"they would share one directory. Diseases: "
                    f"{', '.join(sorted(picks))}")
            # * marks a disease with one study full stop, whose member therefore
            # trains without it. A disease with one *eligible* study out of
            # several is not stranded -- its single-class cohorts stay in
            # training, they just cannot be validated on.
            drawn = " ".join(
                f"{d}/{i['study']}" + ("*" if i["n_studies"] == 1 else "")
                for d, i in picks.items())
            lonely = [d for d, i in picks.items() if i["n_studies"] == 1]
            print(f"[info] {task}/{fold}: {len(units)} disease(s) -> "
                  f"{len(units)} member(s), one study drawn for each: {drawn}"
                  + (f"  (* = {', '.join(lonely)} had one study only, so that "
                     f"member trains without the disease)" if lonely else ""))
            print(f"[info] {task}/{fold}: validation sizes (case/control): "
                  + "  ".join(f"{d} {i['n_case']}/{i['n_ctrl']}"
                              for d, i in picks.items()))
            if len(units) < 2:
                # ensemble_aucs returns NaN below two members, so the fold's
                # ens_* columns would be empty and its answer would rest on
                # the single surviving disease's model.
                print(f"[warn] {task}/{fold}: only {len(units)} member "
                      f"survived --per-disease-min-class "
                      f"{args.per_disease_min_class}, so this fold cannot be "
                      f"ensembled and its ens_* columns will be NaN. Lower the "
                      f"threshold if that is not what you want")
        elif args.inner_split == "same_disease":
            cohorts, test_dis, ccount = survey_same_disease_studies(
                train, meta, args.inner_group, args.disease_col, fold,
                drop_studies=drop)
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

        # Under grouped 'loso' each group shares one split across its
        # --n-estimators repeats, which differ only in model_seed, so it is
        # written once beside their job directories rather than copied into
        # each.
        for i, u in enumerate(units):
            out_dir = os.path.join(t["artifact"], fold, "members", str(u))
            if u in group_split_dir:
                split_dir = group_split_dir[u]
            else:
                split_dir = out_dir
            jobs.append(dict(
                id=f"{task}|{fold}|{u}", task=task, fold=fold, member=i,
                inner_mode=args.inner_split, valid_unit=u,
                split_dir=split_dir,
                valid_studies=group_studies.get(u),
                valid_disease=unit_disease.get(u),
                cv_fold=cv_folds.get(u, (None, None))[0],
                cv_k=cv_folds.get(u, (None, None))[1],
                study_col=args.inner_group,
                split_seed=args.split_seed,
                drop_studies=drop,
                out_dir=out_dir, train=train, test=test, meta=meta,
                head_hidden=args.head_hidden,
                linear_branch=args.linear_branch,
                lin_weight_decay=args.lin_weight_decay,
                glove=args.glove_embedding, select_by=args.select_by,
                abund_mode=args.abund_mode,
                valid_auc=args.valid_auc,
                min_group_n=args.min_group_n,
                logit_adjust_tau=args.logit_adjust_tau,
                group_balance_beta=args.group_balance_beta,
                group_balance_max_ratio=args.group_balance_max_ratio,
                disease_col=args.disease_col,
                member_seed=args.member_seed + i,
                keep_ckpt=args.keep_ckpt, keep_tables=args.keep_tables,
                overwrite=args.overwrite,
                override=dict(model_seed=args.ensemble_seed0 + i,
                              patience=args.patience, min_delta=mdelta,
                              **cfg_over)))
    return jobs, folds


def cleanup_shared_splits(jobs):
    """Delete the shared split tables once every member has had its turn.

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
#: scale: member i's probabilities come from a threshold fitted on the cohort
#: member i held out.
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
        chance on the cohort it held out does not vote. `collect` passes None
        under every remaining --inner-split, since the members are scored on
        different cohorts and their AUCs are not comparable.

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
    early-stopped at epoch 1, or one whose held-out cohort was unlucky, still
    votes under a mean and cannot move a median.

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
        # across members. Every remaining --inner-split scores each member on
        # a different set of cohorts, so a higher valid_auc may only mean an
        # easier one; the weighting is withheld and ens_wauc stays NaN. The
        # modes it was defined for -- one shared validation set, or the same
        # random draw for every member -- were both in-distribution and are
        # gone.
        combo = ensemble_aucs(preds, weights=None)

        aucs = [r["test_auc"] for r in recs]
        eps = sorted(r["epoch"] for r in recs)
        strat = sorted({r.get("strategy") or "" for r in recs} - {""})
        row = dict(fold=fold, abund_mode=args.abund_mode,
                   inner_mode=args.inner_split, strategy=",".join(strat),
                   K=len(recs), n_members=n_members.get(fold),
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

    ap.add_argument("--inner-split", default="loso",
                    choices=["loso", "disease_loso", "per_disease",
                             "same_disease"],
                    help="how each member's training and validation parts are "
                         "cut from the fold's training table. Every mode holds "
                         "out whole cohorts, so the epoch is selected on the "
                         "same cross-study step the test table asks for. "
                         "'loso' (default) gives one member per "
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
                         "'per_disease' cuts that same draw across members "
                         "instead of putting it all in one: one member per "
                         "disease in the fold's training table, each "
                         "validating on a single study drawn for its own "
                         "disease and training on everything else. A member "
                         "therefore gives up one cohort of one disease rather "
                         "than one of each, and it is the ensemble over "
                         "members -- as many as the fold has diseases -- that "
                         "covers every disease. --n-estimators does not apply; "
                         "lodo and loso_all only. "
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
    ap.add_argument("--split-seed", type=int, default=11,
                    help="seed for the --cv-folds cut a single-study disease "
                         "gets under --inner-split same_disease")
    ap.add_argument("--inner-group", default="study",
                    help="metadata column naming the cohort, used by every "
                         "--inner-split to decide what is held out")
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
                    help="members per fold under --inner-split disease_loso, "
                         "where each member redraws which study every disease "
                         "holds out; 8 is the count that mode is meant to be "
                         "run at. A member gives a whole cohort up to "
                         "validation, so it is weaker than a model trained on "
                         "everything, and the ensemble is what buys that back. "
                         "It does not "
                         "apply under 'loso' (unless --loso-k is set), "
                         "'per_disease' or "
                         "'same_disease', where the member count is however "
                         "many cohorts or diseases the fold's training table "
                         "holds")
    ap.add_argument("--drop-studies", nargs="+", default=[], metavar="STUDY",
                    help="remove these cohorts (column --inner-group) from "
                         "every table before anything else happens: the "
                         "training part, the validation part drawn from it, "
                         "and the test table. For a cohort whose labels cannot "
                         "be trusted -- self-reported diagnoses, say -- "
                         "dropping it from training alone is not enough, since "
                         "the fold would still be scored against those labels. "
                         "Applies to every task and every --inner-split, and "
                         "the member plan is drawn up from the filtered table, "
                         "so a disease left with no usable cohort simply gets "
                         "no member. A fold whose test table this empties -- "
                         "under loso_all the fold is a study, so dropping that "
                         "study empties it -- is skipped and reported. Off by "
                         "default; pair it with --run-name, since the results "
                         "are not comparable with a run that kept the cohort")
    ap.add_argument("--per-disease-min-class", type=int, default=20,
                    help="under --inner-split per_disease, a study is only a "
                         "candidate for a disease if its (disease, study) cell "
                         "holds at least this many cases AND at least this "
                         "many controls (default 20; the test is >=, so pass "
                         "21 for a strict 'more than 20'). The epoch is "
                         "selected on that one cell's AUC, and an AUC over a "
                         "handful of either class is mostly sampling noise, so "
                         "the member would early-stop on a coin flip. The "
                         "threshold is per class, not on the cell's total: a "
                         "cohort of 200 controls and 4 cases is large and "
                         "still cannot say whether an epoch separates them. A "
                         "disease with no cohort reaching it gets no member "
                         "and is left out of the fold's ensemble, which is "
                         "reported on a [warn] line. Set 1 for the old "
                         "behaviour of accepting any cell holding both classes")
    ap.add_argument("--member-seed", type=int, default=101,
                    help="seed for the draws the split is made of. Under "
                         "--inner-split disease_loso, member i draws the "
                         "cohorts it holds out with this plus i, so the "
                         "members differ. Under per_disease the draw is one "
                         "per disease and made once per fold, so this seed "
                         "alone decides which study each disease validates on")
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
                         "disease size cancels out of the ratio and is left "
                         "untouched -- that is what --group-balance-beta is "
                         "for. lodo and loso_all only")
    ap.add_argument("--group-balance-beta", type=float, default=0.5,
                    help="how far each disease's contribution to the training "
                         "loss is evened out, used only with --set "
                         "loss=GroupBalanced or loss=GroupBalanced+"
                         "LogitAdjusted. A sample of a disease holding n "
                         "training samples is weighted n**-beta, rescaled so "
                         "the weights average 1 and the loss keeps the scale "
                         "--patience and --min-delta were tuned on. 0 weights "
                         "every sample alike and is exactly BCEWithLogits; 1 "
                         "gives every disease the same total whatever its "
                         "size. Default 0.5, deliberately not the 1.0 that "
                         "--logit-adjust-tau defaults to: an offset does not "
                         "change how hard a sample pulls, a weight multiplies "
                         "it, and this is the only thing in the model that can "
                         "make one sample's gradient twenty times another's. "
                         "On a pool running 40 to 820 samples per disease, "
                         "beta 0.5 lifts the smallest disease from 1.3%% to "
                         "3.7%% of the loss at a 12%% rise in gradient noise, "
                         "where beta 1 reaches 8.3%% but costs 62%%. Read the "
                         "[gbal] lines of a run. lodo and loso_all only")
    ap.add_argument("--group-balance-max-ratio", type=float, default=None,
                    help="cap the per-sample weight ratio between the rarest "
                         "disease and the most common one, under the same two "
                         "losses. Default None keeps the natural "
                         "(n_max/n_min)**beta. The counterpart of "
                         "--min-group-n, for the one case where reweighting "
                         "bites: a disease down to a handful of training "
                         "samples earns a weight big enough that its "
                         "contribution is mostly noise. 20 is a reasonable "
                         "safety net. The disease is still trained on -- this "
                         "only bounds how loud it gets. Prefer lowering "
                         "--group-balance-beta first")
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
    ap.add_argument("--disease-col", default="disease_name_ab",
                    help="metadata column naming each sample's disease. Read "
                         "under --valid-auc macro, which scores each disease "
                         "separately, under --set loss=LogitAdjusted, which "
                         "removes each disease's base rate, and under "
                         "--inner-split disease_loso, which holds out one "
                         "study per disease")
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
    ap.add_argument("--run-tsv", default=None,
                    help="read the fold list from this TSV instead of the "
                         "task's own. The columns are the task's, so this is "
                         "for running a subset of the folds -- a cut-down copy "
                         "of the task's TSV -- not a different task. Safe to "
                         "cut: a job's model_seed is --ensemble-seed0 plus its "
                         "index among its own fold's members, and every member "
                         "of every fold is enumerated from that fold's own "
                         "training table, so dropping rows leaves the "
                         "surviving jobs bit-for-bit the ones the full run "
                         "would have produced. One --tasks only, since which "
                         "task the file belongs to would otherwise be a guess")
    ap.add_argument("--run-name", default="")
    ap.add_argument("--keep-ckpt", action="store_true")
    ap.add_argument("--keep-tables", action="store_true",
                    help="keep each member's training and validation tables "
                         "instead of deleting them on success")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    tasks = list(TASKS) if "all" in args.tasks else args.tasks
    for t in tasks:
        if t not in TASKS:
            raise SystemExit(f"unknown task {t}; choose from {list(TASKS)} or 'all'")
    if args.n_estimators < 1:
        raise SystemExit(f"--n-estimators must be at least 1, got "
                         f"{args.n_estimators}")
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
    if args.valid_auc == "macro" and args.inner_split == "per_disease":
        # Same reason as the single-disease tasks below, reached from the other
        # end: the task's table spans many diseases but this member's
        # validation part holds exactly one, so the average is over one group.
        print("[warn] --valid-auc macro under --inner-split per_disease: a "
              "member validates on one disease and one study of it, so the "
              "macro average is over that one group and equals the pooled "
              "AUC. per_disease already stops one disease deciding the epoch "
              "-- it gives each its own member and ensembles them -- so the "
              "flag adds nothing here. It is what --inner-split disease_loso "
              "needs, whose validation set spans every disease at once")
    if args.valid_auc == "macro":
        single = [t for t in tasks if t not in ("lodo", "loso_all")]
        if single:
            # One disease means one group, so the macro average is that single
            # disease's pooled AUC -- the same number, reached the long way.
            print(f"[warn] --valid-auc macro on {single}: those tasks are "
                  f"single-disease, so the average is over one group and "
                  f"equals the pooled AUC")
    # Membership rather than equality: the combined name carries both halves,
    # so an exact == would let it slip past both guards below and then print a
    # warning saying tau was ignored when it was not.
    _loss_name = parse_set(args.set_cfg).get("loss", CFG["loss"])
    _la = _loss_name in ("LogitAdjusted", "GroupBalanced+LogitAdjusted")
    _gb = _loss_name in ("GroupBalanced", "GroupBalanced+LogitAdjusted")
    if _la:
        bad = [t for t in tasks if t not in ("lodo", "loso_all")]
        if bad:
            raise SystemExit(
                f"loss={_loss_name} needs a training table spanning several "
                f"diseases, and {bad} is single-disease (its metadata is one "
                f"file per disease), where the per-disease base rate collapses "
                f"to one global number. Use --tasks lodo and/or loso_all")
        if args.logit_adjust_tau < 0:
            raise SystemExit(f"--logit-adjust-tau must be >= 0, got "
                             f"{args.logit_adjust_tau}")
        print(f"[cfg] loss={_loss_name} tau={args.logit_adjust_tau:g}: each "
              f"disease's case/control base rate is removed from the training "
              f"logits (column {args.disease_col!r}); validation and test see "
              f"raw logits")
    elif args.logit_adjust_tau != 1.0:
        print("[warn] --logit-adjust-tau only applies with "
              "--set loss=LogitAdjusted or "
              "--set loss=GroupBalanced+LogitAdjusted; ignored")
    if _gb:
        bad = [t for t in tasks if t not in ("lodo", "loso_all")]
        if bad:
            # One disease means one group, so every sample's weight is the
            # same number and the reweighting is an expensive no-op.
            raise SystemExit(
                f"loss={_loss_name} evens out how much each disease "
                f"contributes to the loss, and {bad} is single-disease (its "
                f"metadata is one file per disease), where every sample falls "
                f"in one group and every weight comes out 1. Use --tasks lodo "
                f"and/or loso_all")
        if not 0 <= args.group_balance_beta <= 1:
            raise SystemExit(
                f"--group-balance-beta must lie in [0, 1], got "
                f"{args.group_balance_beta}; below zero it would amplify the "
                f"size imbalance instead of evening it out, and above one a "
                f"rare disease would carry more of the loss in total than a "
                f"common one, which is past the point of equal footing")
        if (args.group_balance_max_ratio is not None
                and args.group_balance_max_ratio < 1):
            raise SystemExit(
                f"--group-balance-max-ratio must be >= 1, got "
                f"{args.group_balance_max_ratio}; below one it would invert "
                f"the cap and hold the rarest disease under the most common")
        if args.select_by == "loss":
            # Not an error: unlike --valid-auc macro with --select-by loss,
            # where the macro number is computed and then ignored, here the
            # training objective genuinely changes and only the selection
            # disagrees with it.
            print("[warn] --select-by loss with a group-balanced training "
                  "loss: the validation loss is pooled and unweighted, since "
                  "the weights come from the training disease counts, so the "
                  "epoch would be chosen on a quantity training is not "
                  "optimising. Prefer --select-by auc (the default)")
        print(f"[cfg] loss={_loss_name} beta={args.group_balance_beta:g}: each "
              f"disease's share of the training loss is rescaled by n^-beta "
              f"(column {args.disease_col!r}), normalised to mean 1"
              + (f", capped at {args.group_balance_max_ratio:g}x"
                 if args.group_balance_max_ratio is not None else ""))
        # Only worth suggesting where validation actually spans several
        # diseases. Under per_disease a member validates on one disease, so a
        # macro average would be over a single group and equal the pooled AUC.
        if args.inner_split == "per_disease":
            print("[cfg] validation loss stays pooled and unweighted, but "
                  "--inner-split per_disease already keeps one disease from "
                  "deciding the epoch: each member is selected on its own "
                  "disease and the fold is their ensemble. --valid-auc macro "
                  "would do nothing on top of that")
        else:
            print("[cfg] validation loss stays pooled and unweighted; pair "
                  "this with --valid-auc macro so the epoch is chosen per "
                  "disease too")
    elif args.group_balance_beta != 0.5:
        print("[warn] --group-balance-beta only applies with "
              "--set loss=GroupBalanced or "
              "--set loss=GroupBalanced+LogitAdjusted; ignored")
    elif args.group_balance_max_ratio is not None:
        print("[warn] --group-balance-max-ratio only applies with "
              "--set loss=GroupBalanced or "
              "--set loss=GroupBalanced+LogitAdjusted; ignored")
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
    if args.inner_split == "per_disease":
        bad = [t for t in tasks if t not in ("lodo", "loso_all")]
        if bad:
            raise SystemExit(
                f"--inner-split per_disease runs one member per disease in the "
                f"training table, and {bad} is single-disease (its metadata is "
                f"one file per disease), so every fold would come down to one "
                f"member holding out one random study. Use --tasks lodo and/or "
                f"loso_all, or --inner-split loso for {bad}")
        if args.inner_group == args.disease_col:
            # Every disease would then be its own only cohort, so the drawn
            # "study" is the whole disease and the member validates on a
            # disease it never trains on -- leave-one-disease-out nested inside
            # the fold, not the cross-study step this mode is for.
            raise SystemExit(
                f"--inner-group and --disease-col are both "
                f"{args.inner_group!r}. Under --inner-split per_disease they "
                f"name different things: --inner-group is the cohort/study "
                f"column (default 'study') and --disease-col is the disease "
                f"column (default 'disease_name_ab'). Pointing both at one "
                f"column makes each disease its own only cohort, so the "
                f"member holds its whole disease out instead of one of its "
                f"studies. Drop --inner-group to take the default")
        if args.per_disease_min_class < 1:
            raise SystemExit(
                f"--per-disease-min-class must be at least 1, got "
                f"{args.per_disease_min_class}; a cell needs one of each class "
                f"before an AUC exists at all")
        if args.n_estimators != 1:
            print(f"[info] --n-estimators {args.n_estimators} does not apply "
                  f"under --inner-split per_disease: a fold runs one member "
                  f"per disease its training table holds, so the member count "
                  f"is the fold's and not a flag's")
    elif args.per_disease_min_class != 20:
        print("[warn] --per-disease-min-class only applies under "
              "--inner-split per_disease; ignored")
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
    if args.drop_studies:
        print(f"[cfg] --drop-studies {' '.join(args.drop_studies)}: removed "
              f"from the training, validation and test tables alike (column "
              f"{args.inner_group!r}). Member plans are drawn from the "
              f"filtered table, and a fold whose test table this empties is "
              f"skipped")
        if not args.run_name:
            print("[cfg] WARNING: --drop-studies without --run-name writes "
                  "into the same directories as a run that kept those "
                  "cohorts, and the two are not comparable")
    if args.run_tsv is not None:
        if len(tasks) != 1:
            raise SystemExit(
                f"--run-tsv replaces one task's fold list, and --tasks names "
                f"{len(tasks)} ({', '.join(tasks)}), so which one it belongs "
                f"to would be a guess. Run the tasks one at a time")
        if not os.path.exists(args.run_tsv):
            raise SystemExit(f"--run-tsv {args.run_tsv!r} does not exist")
        TASKS[tasks[0]]["run_tsv"] = args.run_tsv

    if args.run_name:
        for t in tasks:
            TASKS[t]["artifact"] = f"{TASKS[t]['artifact']}{args.run_name}"

    print(f"[cfg] tasks={tasks} gpus={args.gpus} x{args.workers_per_gpu}")
    for t in tasks:
        # On the record above every run, so a set of results carries the fold
        # list it was made from -- a cut-down TSV is otherwise invisible
        # afterwards, and the folds it dropped look like folds that failed.
        print(f"[cfg] {t}: folds from {TASKS[t]['run_tsv']}"
              + ("  <- --run-tsv" if args.run_tsv is not None else ""))
    if args.inner_split == "same_disease":
        print(f"[cfg] inner-split=same_disease: one member per study (column "
              f"{args.inner_group!r}) of the test study's disease (column "
              f"{args.disease_col!r}); each validates on that one study and "
              f"trains on the disease's remaining studies plus every other "
              f"disease. The member count is the fold's, not a flag's")
        print(f"[cfg] a disease down to one study is cut --cv-folds="
              f"{args.cv_folds} ways instead, so it stays in training")
    elif args.inner_split == "per_disease":
        print(f"[cfg] inner-split=per_disease: one member per disease (column "
              f"{args.disease_col!r}) in the fold's training table; member d "
              f"validates on one study (column {args.inner_group!r}) drawn for "
              f"d alone and trains on everything else. The member count is the "
              f"fold's, not a flag's, and the fold's answer is the ensemble "
              f"over them")
        print(f"[cfg] the draw uses --member-seed {args.member_seed} and is "
              f"made once per fold; a study is a candidate only if its "
              f"(disease, study) cell holds >= {args.per_disease_min_class} "
              f"case(s) and >= {args.per_disease_min_class} control(s) "
              f"(--per-disease-min-class), and a disease with no such study "
              f"gets no member and drops out of that fold's ensemble")
        print("[cfg] the unit held out is the (disease, study) cell, so a "
              "study spanning several diseases keeps its other diseases in "
              "training and does sit on both sides of that member's split")
    elif args.inner_split == "disease_loso":
        print(f"[cfg] inner-split=disease_loso: {args.n_estimators} members "
              f"per fold, each holding out one randomly drawn study (column "
              f"{args.inner_group!r}) per disease (column {args.disease_col!r}) "
              f"as its validation set. Seeds member_seed={args.member_seed}.."
              f"{args.member_seed + args.n_estimators - 1}")
        print("[cfg] the unit held out is the (disease, study) cell, so a "
              "study spanning several diseases keeps its other diseases in "
              "training and does sit on both sides of that member's split")
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
    if args.inner_split not in ("loso", "per_disease"):
        # Both of those size their member count from the fold rather than from
        # --n-estimators, so the seed range is only known once the folds are
        # surveyed; it is printed below instead.
        print(f"[cfg] model_seed={args.ensemble_seed0}.."
              f"{args.ensemble_seed0 + args.n_estimators - 1}")
    print(f"[cfg] abund_mode={args.abund_mode} "
          f"linear_branch={args.linear_branch} select_by={args.select_by} "
          f"patience={args.patience} valid_auc={args.valid_auc}"
          + (f" min_group_n={args.min_group_n}"
             if args.valid_auc == "macro" else "")
          + (f" group_balance_beta={args.group_balance_beta:g}"
             + (f" max_ratio={args.group_balance_max_ratio:g}"
                if args.group_balance_max_ratio is not None else "")
             if _gb else ""))
    # Parsed here as well as in build_jobs so a typo costs a second, not a
    # GPU-hour, and so the deviation from CFG is on the record above the run.
    cfg_over = parse_set(args.set_cfg)
    if cfg_over:
        print("[cfg] --set overrides: " + ", ".join(
            f"{k}={v!r} (default {CFG[k]!r})" for k, v in cfg_over.items()))
        if not args.run_name:
            print("[cfg] WARNING: --set without --run-name writes this "
                  "setting into the same directories as every other setting")
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
    if args.inner_split in ("loso", "per_disease") and per_fold:
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
