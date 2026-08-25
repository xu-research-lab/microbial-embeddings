#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Support Vector Machine protocol across all LOSO tasks, sequential execution.

This script mirrors the structure of run_rf.py:
  Step 1  Inner leave-one-study-out: K SVM models, each holding out one
          training cohort (or one study-level CV fold) as validation.
  Step 2  Retrain one SVM model on ALL training cohorts per seed.
  Bonus   Probability-average ensemble of the K inner models.

Fallback when the training set has too few cohorts
--------------------------------------------------
Same logic as the Attention/RF scripts: when a fold's training table carries
fewer than --cv-fallback-below usable cohorts, the inner split falls back
to stratified K-fold cross-validation over samples.

When the training table has more than --study-cv-above cohorts, the inner
split switches from LOSO to 5-fold CV at the study level, training exactly
five inner models for the ensemble.

Usage
-----
  python run_svm.py --dry-run
  python run_svm.py --tasks disease
  python run_svm.py --tasks all --seeds 11 22 33
  python run_svm.py --tasks all --skip-retrain
"""
from __future__ import annotations

import argparse
import json
import os
import time
import traceback

import biom
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder
from sklearn.svm import SVC

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
        results="svm_disease.csv",
    ),
    "ibd_subtype": dict(
        run_tsv="run_leave_one_study_out_IBD_subtype_list.tsv",
        train="Data/IBD_subtype_data/{disease}/{study}/train_loo.biom",
        test="Data/IBD_subtype_data/{disease}/{study}/test_loo.biom",
        meta="Data/IBD_subtype_data/{disease}/metadata.tsv",
        fold="{disease}_{study}",
        artifact="Data/IBD_subtype_data/",
        results="svm_ibd_subtype.csv",
    ),
    "lodo": dict(
        run_tsv="run_leave_one_disease_out.tsv",
        train="Data/loo_all_diseases/data/{left_out_disease}/train_loo.biom",
        test="Data/loo_all_diseases/data/{left_out_disease}/test_loo.biom",
        meta="Data/loo_all_diseases/data/metadata.tsv",
        fold="{left_out_disease}",
        artifact="Data/loo_all_diseases/",
        results="svm_lodo.csv",
    ),
    "loso_all": dict(
        run_tsv="run_leave_one_study_out.tsv",
        train="Data/loo_all_studies/data/{left_out_study}/train_loo.biom",
        test="Data/loo_all_studies/data/{left_out_study}/test_loo.biom",
        meta="Data/loo_all_studies/data/metadata.tsv",
        fold="{left_out_study}",
        artifact="Data/loo_all_studies/",
        results="svm_loso_all.csv",
    ),
}

SAMPLE_ID_COL = "sample"
LABELS_COL = "group"
STUDY_COL_CANDIDATES = ["study", "Study", "dataset", "Dataset",
                        "project", "bioproject", "cohort", "batch"]


# =====================================================================
# 1. Data helpers
# =====================================================================
def _load_normalized(biom_path):
    """Load a BIOM table and apply rankdata normalisation.

    Returns the BIOM Table object; the caller extracts sample ids and then
    converts to a dense sample-by-feature matrix.
    """
    table = biom.load_table(biom_path)
    return table.rankdata(axis="sample", inplace=False)


def _encode_labels(y_train, y_test=None):
    """Encode binary labels to {0, 1} based on the training set."""
    le = LabelEncoder()
    y_train_enc = le.fit_transform(y_train)
    if y_test is not None:
        y_test_enc = le.transform(y_test)
        return y_train_enc, y_test_enc, le
    return y_train_enc, le


def _auc(y_true, y_prob):
    """Return ROC-AUC if both classes are present, else NaN."""
    y_true = np.asarray(y_true)
    if len(np.unique(y_true)) < 2:
        return float(np.nan)
    return float(roc_auc_score(y_true, y_prob))


def _save_pred(path, sample_ids, true_labels, probs):
    """Write predictions in the same format as the Attention script."""
    pd.DataFrame({
        "sample_id": sample_ids,
        "true_label": true_labels,
        "prob": probs,
    }).to_csv(path, index=False)


# =====================================================================
# 2. Split helpers
# =====================================================================
def _check_split(lab, ti, vi):
    if len(vi) == 0 or len(ti) == 0:
        raise ValueError("split leaves one side empty")
    for nm, idx in (("train", ti), ("valid", vi)):
        if len(np.unique(lab[idx])) < 2:
            raise ValueError(f"{nm} split has only one class")
    return ti, vi


def split_by_study(labels, groups, valid_study):
    """Validation set = the held-out cohort(s)."""
    g = np.asarray(groups).astype(str)
    lab = np.asarray(labels)
    if isinstance(valid_study, (list, tuple, np.ndarray)):
        held = set(str(x) for x in valid_study)
        is_v = np.isin(g, list(held))
    else:
        is_v = g == str(valid_study)
    return _check_split(lab, np.where(~is_v)[0], np.where(is_v)[0])


def split_by_cv(labels, cv_splits, cv_index, cv_seed):
    """Validation set = one stratified CV fold over samples."""
    lab = np.asarray(labels)
    skf = StratifiedKFold(n_splits=cv_splits, shuffle=True,
                          random_state=cv_seed)
    ti, vi = list(skf.split(np.zeros(len(lab)), lab))[cv_index]
    return _check_split(lab, ti, vi)


# =====================================================================
# 3. Worker: run one SVM job
# =====================================================================
def _build_svm(job):
    return SVC(
        C=job["C"],
        kernel=job["kernel"],
        gamma=job["gamma"],
        probability=True,
        random_state=job["seed"])


def run_job(job):
    """Train one SVM model and write result.json."""
    os.makedirs(job["out_dir"], exist_ok=True)
    res_path = os.path.join(job["out_dir"], "result.json")
    if os.path.exists(res_path):
        with open(res_path) as f:
            return job["id"], json.load(f), None

    try:
        train_table = _load_normalized(job["train"])
        test_table = _load_normalized(job["test"])
        metadata = pd.read_csv(job["meta"], sep="\t", index_col=SAMPLE_ID_COL,
                               dtype={SAMPLE_ID_COL: str}, low_memory=False)

        train_idx = train_table.ids(axis="sample")
        test_idx = test_table.ids(axis="sample")

        X_train_all = train_table.matrix_data.multiply(
            1 / train_table.max(axis="sample")).toarray().T
        X_test = test_table.matrix_data.multiply(
            1 / test_table.max(axis="sample")).toarray().T

        y_train_raw = metadata.loc[train_idx, LABELS_COL]
        y_test_raw = metadata.loc[test_idx, LABELS_COL]
        y_train_all, y_test, _ = _encode_labels(y_train_raw, y_test_raw)

        svm = _build_svm(job)

        if job["kind"] == "inner":
            if job["inner_mode"] == "cv":
                ti, vi = split_by_cv(y_train_all, job["cv_splits"],
                                     job["cv_index"], job["cv_seed"])
            else:
                groups = metadata.loc[train_idx, job["study_col"]]
                ti, vi = split_by_study(y_train_all, groups,
                                        job["valid_study"])

            X_tr, y_tr = X_train_all[ti], y_train_all[ti]
            X_va, y_va = X_train_all[vi], y_train_all[vi]

            svm.fit(X_tr, y_tr)

            y_prob_test = svm.predict_proba(X_test)[:, 1]
            y_prob_valid = svm.predict_proba(X_va)[:, 1]

            _save_pred(os.path.join(job["out_dir"], "pred_test.csv"),
                       test_idx, y_test, y_prob_test)
            _save_pred(os.path.join(job["out_dir"], "pred_valid.csv"),
                       train_idx[vi], y_va, y_prob_valid)

            out = dict(kind="inner", task=job["task"], fold=job["fold"],
                       inner_mode=job["inner_mode"],
                       valid_unit=str(job["valid_unit"]),
                       valid_auc=_auc(y_va, y_prob_valid),
                       n_valid=int(len(vi)),
                       test_auc=_auc(y_test, y_prob_test),
                       seed=job["seed"])
        else:
            svm.fit(X_train_all, y_train_all)
            y_prob_test = svm.predict_proba(X_test)[:, 1]

            _save_pred(os.path.join(job["out_dir"], "pred_test.csv"),
                       test_idx, y_test, y_prob_test)

            out = dict(kind="retrain", task=job["task"], fold=job["fold"],
                       seed=job["seed"],
                       test_auc=_auc(y_test, y_prob_test))

        with open(res_path, "w") as f:
            json.dump(out, f)
        return job["id"], out, None
    except Exception:
        return job["id"], None, traceback.format_exc()


# =====================================================================
# 4. Job construction (parent process)
# =====================================================================
def survey_training_table(train_biom, meta_path, study_col):
    """Return the cohort names present in the training table."""
    sids = biom.load_table(train_biom).ids(axis="sample")
    md = pd.read_csv(meta_path, sep="\t", index_col=SAMPLE_ID_COL,
                     dtype={SAMPLE_ID_COL: str}, low_memory=False)
    return md.loc[sids, study_col].unique()


def build_inner_jobs(task, args):
    t = TASKS[task]
    run = pd.read_csv(t["run_tsv"], sep="\t")
    jobs, meta_by_fold = [], {}

    for _, r in run.iterrows():
        row = r.to_dict()
        fold = t["fold"].format(**row)
        train = t["train"].format(**row)
        test = t["test"].format(**row)
        meta = t["meta"].format(**row)

        missing = [p for p in (train, test, meta) if not os.path.exists(p)]
        if missing:
            print(f"[warn] {task}/{fold}: missing {missing}, skipping fold")
            continue

        sc = args.inner_group
        studies = survey_training_table(train, meta, sc)

        if len(studies) > args.study_cv_above:
            mode = "study_cv"
            n_study_folds = 5
            studies_sorted = sorted(studies)
            study_folds = [studies_sorted[i::n_study_folds]
                           for i in range(n_study_folds)]
            units = [f"study_cv{i}" for i in range(n_study_folds)]
            print(f"[info] {task}/{fold}: {len(studies)} usable cohorts -> "
                  f"{n_study_folds}-fold CV at study level")
        elif len(studies) >= args.cv_fallback_below:
            mode, units, study_folds = "loso", studies, None
        else:
            if args.no_cv_fallback:
                print(f"[warn] {task}/{fold}: only {len(studies)} usable "
                      f"cohort(s) and CV fallback disabled, skipping fold")
                continue
            mode, units = "cv", [f"cv{i}" for i in range(args.cv_splits)]
            study_folds = None
            print(f"[info] {task}/{fold}: only {len(studies)} usable "
                  f"cohort(s) -> {args.cv_splits}-fold stratified CV inner split "
                  f"(validation is in-distribution; read with care)")

        meta_by_fold[(task, fold)] = dict(
            train=train, test=test, meta=meta, study_col=sc,
            units=units, inner_mode=mode, artifact=t["artifact"],
            study_folds=study_folds)

        for i, u in enumerate(units):
            out_dir = os.path.join(t["artifact"], fold, "inner", str(u))
            jobs.append(dict(
                id=f"{task}|{fold}|inner|{u}", kind="inner", task=task,
                fold=fold, inner_mode=mode, valid_unit=u, out_dir=out_dir,
                train=train, test=test, meta=meta, study_col=sc,
                C=args.C, kernel=args.kernel, gamma=args.gamma,
                seed=args.seed,
                valid_study=(study_folds[i] if mode == "study_cv"
                             else (u if mode == "loso" else None)),
                cv_index=i if mode == "cv" else None,
                cv_splits=len(units) if mode == "cv" else None,
                cv_seed=args.cv_seed))

    return jobs, meta_by_fold


def build_retrain_jobs(inner_results, meta_by_fold, args):
    """Queue retrain jobs: one SVM per seed trained on all cohorts."""
    by_fold = {}
    for r in inner_results:
        by_fold.setdefault((r["task"], r["fold"]), []).append(r)

    jobs = []
    for key, recs in sorted(by_fold.items()):
        if key not in meta_by_fold:
            continue
        task, fold = key
        m = meta_by_fold[key]
        for seed in args.seeds:
            out_dir = os.path.join(m["artifact"], fold, "retrain",
                                   f"seed{seed}")
            jobs.append(dict(
                id=f"{task}|{fold}|retrain|s{seed}", kind="retrain",
                task=task, fold=fold, seed=seed, out_dir=out_dir,
                train=m["train"], test=m["test"], meta=m["meta"],
                study_col=m["study_col"],
                C=args.C, kernel=args.kernel, gamma=args.gamma))
    return jobs


# =====================================================================
# 5. Sequential execution
# =====================================================================
def run_jobs(jobs, label):
    if not jobs:
        print(f"[{label}] nothing to run")
        return []

    print(f"[{label}] {len(jobs)} jobs (sequential)")
    results, failed, t0 = [], [], time.time()
    for i, job in enumerate(jobs, 1):
        jid, out, err = run_job(job)
        el = time.time() - t0
        eta = el / i * (len(jobs) - i)
        if err is None:
            results.append(out)
            extra = (f"valid_auc={out['valid_auc']:.4f}"
                     if out["kind"] == "inner"
                     else "retrain")
            print(f"[{label} {i}/{len(jobs)}] {jid}  "
                  f"test_AUC={out['test_auc']:.4f} {extra}  "
                  f"| elapsed {el/60:.0f}m eta {eta/60:.0f}m", flush=True)
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
# 6. Collection
# =====================================================================
def prob_avg_auc(pred_paths):
    """Average predicted probabilities across models and compute AUC."""
    dfs = [pd.read_csv(p) for p in pred_paths]
    ref = dfs[0]
    y = ref["true_label"].to_numpy()
    if len(np.unique(y)) < 2:
        return float(np.nan)
    probs = []
    for d in dfs:
        if not np.array_equal(d["sample_id"].to_numpy(),
                              ref["sample_id"].to_numpy()):
            return float(np.nan)
        probs.append(d["prob"].to_numpy())
    return float(roc_auc_score(y, np.vstack(probs).mean(axis=0)))


def collect(task, meta_by_fold, inner_results, retrain_results, args):
    by_fold = {}
    for r in inner_results:
        if r["task"] == task:
            by_fold.setdefault(r["fold"], []).append(r)
    re_by = {}
    for r in retrain_results:
        if r["task"] == task:
            re_by.setdefault(r["fold"], {})[r["seed"]] = r["test_auc"]

    rows = []
    for fold, recs in sorted(by_fold.items()):
        key = (task, fold)
        if key not in meta_by_fold:
            continue
        m = meta_by_fold[key]
        preds = [os.path.join(m["artifact"], fold, "inner",
                              r["valid_unit"], "pred_test.csv") for r in recs]
        preds = [p for p in preds if os.path.exists(p)]
        row = dict(fold=fold, inner_mode=m["inner_mode"], K=len(recs),
                   inner_mean=float(np.mean([r["test_auc"] for r in recs])),
                   ens_prob=prob_avg_auc(preds) if len(preds) >= 2 else float(np.nan))
        for s in args.seeds:
            row[f"retrain_seed{s}"] = re_by.get(fold, {}).get(s, float(np.nan))
        rows.append(row)

    if not rows:
        return None
    df = pd.DataFrame(rows)
    sc = [c for c in df.columns if c.startswith("retrain_seed")]
    df["retrain"] = df[sc].mean(axis=1)
    df["retrain_sd"] = df[sc].std(axis=1) if len(sc) > 1 else float(np.nan)

    out = TASKS[task]["results"]
    if args.run_name:
        out = out.replace(".csv", f".{args.run_name}.csv")
    df.to_csv(out, index=False)

    print(f"\n=== {task} -> {out} ===")
    cols = ["fold", "inner_mode", "K", "ens_prob",
            "inner_mean", "retrain", "retrain_sd"]
    print(df[cols].to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print(f"  mean: ens_prob {df.ens_prob.mean():.4f}  "
          f"inner_mean {df.inner_mean.mean():.4f}  "
          f"retrain {df.retrain.mean():.4f}")
    if (df.inner_mode == "cv").any():
        sub = df[df.inner_mode == "cv"]
        print(f"  note: {len(sub)}/{len(df)} fold(s) used the CV fallback "
              f"(in-distribution validation, lower ensemble diversity): "
              f"{', '.join(sub.fold)}")
    return df


# =====================================================================
# 7. Main
# =====================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tasks", nargs="+", default=["disease"],
                    help=f"{list(TASKS)} or 'all'")
    ap.add_argument("--seeds", type=int, nargs="+", default=[11],
                    help="seeds used for the retrain phase")
    ap.add_argument("--inner-group", default="study",
                    help="cohort column for the inner split")
    ap.add_argument("--study-cv-above", type=int, default=15,
                    help="when the training table has more cohorts than this, "
                         "use 5-fold CV at the study level instead of LOSO")
    ap.add_argument("--cv-fallback-below", type=int, default=2,
                    help="use stratified CV when the training table has fewer "
                         "usable cohorts than this")
    ap.add_argument("--cv-splits", type=int, default=5,
                    help="number of CV folds used by the fallback")
    ap.add_argument("--cv-seed", type=int, default=11)
    ap.add_argument("--no-cv-fallback", action="store_true",
                    help="skip such folds instead of falling back to CV")
    ap.add_argument("--run-name", default="",
                    help="suffix added to result directories and CSV files")
    ap.add_argument("--skip-retrain", action="store_true",
                    help="inner models and ensemble only")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--C", type=float, default=1.0,
                    help="regularisation parameter passed to SVC")
    ap.add_argument("--kernel", default="rbf",
                    choices=["linear", "poly", "rbf", "sigmoid"],
                    help="kernel passed to SVC")
    ap.add_argument("--gamma", default="scale",
                    help="kernel coefficient passed to SVC")
    ap.add_argument("--seed", type=int, default=11,
                    help="random seed used for the inner SVM models")
    args = ap.parse_args()

    tasks = list(TASKS) if "all" in args.tasks else args.tasks
    for t in tasks:
        if t not in TASKS:
            raise SystemExit(f"unknown task {t}; choose from {list(TASKS)} or 'all'")
    if args.run_name:
        for t in tasks:
            TASKS[t]["artifact"] = f"{TASKS[t]['artifact']}_{args.run_name}"

    print(f"[cfg] tasks={tasks} seeds={args.seeds}")
    print(f"[cfg] C={args.C} kernel={args.kernel} gamma={args.gamma} "
          f"seed={args.seed} (sequential execution)")
    print(f"[cfg] study_cv_above={args.study_cv_above} "
          f"cv_fallback_below={args.cv_fallback_below} "
          f"cv_splits={args.cv_splits} cv_seed={args.cv_seed}")

    # ---- Phase A ----
    all_inner, meta_by_fold = [], {}
    for t in tasks:
        j, m = build_inner_jobs(t, args)
        all_inner += j
        meta_by_fold.update(m)
    todo_inner = [j for j in all_inner
                  if not os.path.exists(os.path.join(j["out_dir"],
                                                     "result.json"))]
    print(f"\nPhase A (inner): {len(all_inner)} jobs, "
          f"{len(all_inner) - len(todo_inner)} done, {len(todo_inner)} to run")

    if args.dry_run:
        for t in tasks:
            n = sum(1 for j in all_inner if j["task"] == t)
            folds = [k for k in meta_by_fold if k[0] == t]
            n_cv = sum(1 for k in folds if meta_by_fold[k]["inner_mode"] == "cv")
            n_study_cv = sum(1 for k in folds if meta_by_fold[k]["inner_mode"] == "study_cv")
            print(f"  {t:12s}: {len(folds)} folds "
                  f"({len(folds) - n_cv - n_study_cv} LOSO, "
                  f"{n_study_cv} study-CV, {n_cv} CV fallback), "
                  f"{n} inner jobs")
        print(f"  at 1 min per job, phase A takes about "
              f"{len(todo_inner) / 60:.1f} h")
        return

    run_jobs(todo_inner, "A/inner")

    inner_results = []
    for j in all_inner:
        p = os.path.join(j["out_dir"], "result.json")
        if os.path.exists(p):
            with open(p) as f:
                inner_results.append(json.load(f))
    print(f"\nPhase A done: {len(inner_results)}/{len(all_inner)} inner results")

    # ---- Phase B ----
    retrain_results = []
    re_jobs = build_retrain_jobs(inner_results, meta_by_fold, args)
    if not args.skip_retrain:
        todo_re = [j for j in re_jobs
                   if not os.path.exists(os.path.join(j["out_dir"],
                                                      "result.json"))]
        print(f"\nPhase B (retrain): {len(re_jobs)} jobs, "
              f"{len(todo_re)} to run")
        run_jobs(todo_re, "B/retrain")
        for j in re_jobs:
            p = os.path.join(j["out_dir"], "result.json")
            if os.path.exists(p):
                with open(p) as f:
                    retrain_results.append(json.load(f))

    # ---- Collect ----
    for t in tasks:
        collect(t, meta_by_fold, inner_results, retrain_results, args)


if __name__ == "__main__":
    main()
