#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Random Forest across all LOSO tasks — simplest possible protocol.

For every fold:
  1. normalise each table: rankdata + norm per sample, over the full
     feature set,
  2. select: keep only the TOP-K largest entries of each sample (the rest
     are set to zero), then drop features that are empty everywhere,
  3. merge train and test onto one shared feature table,
  4. train one RF on the whole training table, evaluate on the test table.

No inner splits, no CV fallback, no ensembles, no retrain phase.

Usage
-----
  python run_rf.py --dry-run
  python run_rf.py --tasks disease
  python run_rf.py --tasks all --n-jobs 16 --top-k 600
  python run_rf.py --tasks all --top-k 0          # disable the filter
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
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import LabelEncoder

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
        results="rf_disease.csv",
    ),
    "ibd_subtype": dict(
        run_tsv="run_leave_one_study_out_IBD_subtype_list.tsv",
        train="Data/IBD_subtype_data/{disease}/{study}/train_loo.biom",
        test="Data/IBD_subtype_data/{disease}/{study}/test_loo.biom",
        meta="Data/IBD_subtype_data/{disease}/metadata.tsv",
        fold="{disease}_{study}",
        artifact="Data/IBD_subtype_data/",
        results="rf_ibd_subtype.csv",
    ),
    "lodo": dict(
        run_tsv="run_leave_one_disease_out.tsv",
        train="Data/loo_all_diseases/data/{left_out_disease}/train_loo.biom",
        test="Data/loo_all_diseases/data/{left_out_disease}/test_loo.biom",
        meta="Data/loo_all_diseases/data/metadata.tsv",
        fold="{left_out_disease}",
        artifact="Data/loo_all_diseases/",
        results="rf_lodo.csv",
    ),
    "loso_all": dict(
        run_tsv="run_leave_one_study_out.tsv",
        train="Data/loo_all_studies/data/{left_out_study}/train_loo.biom",
        test="Data/loo_all_studies/data/{left_out_study}/test_loo.biom",
        meta="Data/loo_all_studies/data/metadata.tsv",
        fold="{left_out_study}",
        artifact="Data/loo_all_studies/",
        results="rf_loso_all.csv",
    ),
}

SAMPLE_ID_COL = "sample"
LABELS_COL = "group"


# =====================================================================
# 1. Data helpers
# =====================================================================
def _keep_top_k_per_sample(M, k):
    """Zero every entry outside the k largest values of each sample.

    M is an observation-by-sample dense matrix (BIOM's native orientation),
    so "each sample" means "each column".  Ties at the cut-off are broken
    arbitrarily but deterministically by argpartition.
    """
    n_feat = M.shape[0]
    if not k or k >= n_feat:
        return M
    # rows *outside* the top-k of every column
    drop_rows = np.argpartition(-M, k - 1, axis=0)[k:]
    np.put_along_axis(M, drop_rows, 0.0, axis=0)
    return M


def _load_one(biom_path, top_k, renorm):
    """Step 1 + 2 for a single table.

    1. normalise: rankdata + L1 norm per sample, computed on the table's
       FULL feature set (identical to the original
       `table.rankdata(axis='sample').norm(axis='sample')`);
    2. select: keep only the top-k entries of each sample, zero the rest,
       then drop features that are empty in every sample.

    Returns a sample-by-feature DataFrame.
    """
    table = biom.load_table(biom_path)
    # table = table.rankdata(axis="sample", inplace=False)
    # table = table.norm(axis="sample", inplace=False)

    obs_ids = np.asarray(table.ids(axis="observation"), dtype=object)
    sample_ids = np.asarray(table.ids(axis="sample"), dtype=object)
    
    ranked = table.rankdata(axis="sample", inplace=False)
    max_values = np.asarray(ranked.max(axis="sample"), dtype=float)   
    M = ranked.matrix_data.multiply(1.0 / max_values).toarray()   # obs x sample

    M = _keep_top_k_per_sample(M, top_k)
    # if renorm:                                # restore per-sample sum == 1
    #     col = M.sum(axis=0, keepdims=True)
    #     col[col == 0] = 1.0
    #     M = M / col

    keep = M.sum(axis=1) > 0                                  # empty features
    return pd.DataFrame(M[keep].T, index=sample_ids, columns=obs_ids[keep])


def _build_fold_matrix(train_path, test_path, top_k, feature_space, renorm):
    """Step 3: merge the two already-normalised, already-filtered tables
    onto one shared feature space."""
    tr = _load_one(train_path, top_k, renorm)
    te = _load_one(test_path, top_k, renorm)

    if feature_space == "union":
        feats = tr.columns.union(te.columns, sort=False)
    elif feature_space == "intersection":
        feats = tr.columns.intersection(te.columns, sort=False)
    else:                                     # "train" (default)
        feats = tr.columns

    tr = tr.reindex(columns=feats, fill_value=0.0)
    te = te.reindex(columns=feats, fill_value=0.0)

    return (tr.index.to_numpy(), tr.to_numpy(dtype=float),
            te.index.to_numpy(), te.to_numpy(dtype=float),
            len(feats))


def _encode_labels(y_train, y_test):
    """Encode binary labels to {0, 1} based on the training set."""
    le = LabelEncoder()
    return le.fit_transform(y_train), le.transform(y_test), le


def _auc(y_true, y_prob):
    """Return ROC-AUC if both classes are present, else NaN."""
    y_true = np.asarray(y_true)
    if len(np.unique(y_true)) < 2:
        return float(np.nan)
    return float(roc_auc_score(y_true, y_prob))


def _save_pred(path, sample_ids, true_labels, probs):
    pd.DataFrame({
        "sample_id": sample_ids,
        "true_label": true_labels,
        "prob": probs,
    }).to_csv(path, index=False)


# =====================================================================
# 2. Worker: train on train, evaluate on test
# =====================================================================
def run_job(job):
    """Train one RF model on the full training table and write result.json."""
    os.makedirs(job["out_dir"], exist_ok=True)
    res_path = os.path.join(job["out_dir"], "result.json")
    if os.path.exists(res_path) and not job.get("overwrite"):
        with open(res_path) as f:
            return job["id"], json.load(f), None

    try:
        train_idx, X_train, test_idx, X_test, n_feat = _build_fold_matrix(
            job["train"], job["test"], job["top_k"], job["feature_space"],
            job["renorm_after_topk"])

        metadata = pd.read_csv(job["meta"], sep="\t", index_col=SAMPLE_ID_COL,
                               dtype={SAMPLE_ID_COL: str}, low_memory=False)
        y_train, y_test, _ = _encode_labels(
            metadata.loc[train_idx, LABELS_COL],
            metadata.loc[test_idx, LABELS_COL])

        rf = RandomForestClassifier(
            n_estimators=job["n_estimators"],
            max_depth=job["max_depth"],
            random_state=job["seed"],
            n_jobs=job["n_jobs"])
        rf.fit(X_train, y_train)
        y_prob_test = rf.predict_proba(X_test)[:, 1]

        _save_pred(os.path.join(job["out_dir"], "pred_test.csv"),
                   test_idx, y_test, y_prob_test)

        out = dict(task=job["task"], fold=job["fold"], seed=job["seed"],
                   top_k=job["top_k"], n_features=int(n_feat),
                   n_train=int(len(train_idx)), n_test=int(len(test_idx)),
                   test_auc=_auc(y_test, y_prob_test))
        with open(res_path, "w") as f:
            json.dump(out, f)
        return job["id"], out, None
    except Exception:
        return job["id"], None, traceback.format_exc()


# =====================================================================
# 3. Job construction
# =====================================================================
def build_jobs(task, args):
    t = TASKS[task]
    run = pd.read_csv(t["run_tsv"], sep="\t")
    jobs = []
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

        jobs.append(dict(
            id=f"{task}|{fold}", task=task, fold=fold,
            out_dir=os.path.join(t["artifact"], fold, "rf"),
            train=train, test=test, meta=meta,
            top_k=args.top_k, feature_space=args.feature_space,
            renorm_after_topk=args.renorm_after_topk,
            n_estimators=args.n_estimators, max_depth=args.max_depth,
            n_jobs=args.n_jobs, seed=args.seed,
            overwrite=args.overwrite))
    return jobs


# =====================================================================
# 4. Sequential execution
# =====================================================================
def run_jobs(jobs, label="rf"):
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
            print(f"[{label} {i}/{len(jobs)}] {jid}  "
                  f"test_AUC={out['test_auc']:.4f}  "
                  f"F={out['n_features']}  "
                  f"| elapsed {el/60:.0f}m eta {eta/60:.0f}m", flush=True)
        else:
            failed.append(jid)
            print(f"[{label} {i}/{len(jobs)}] {jid}  FAILED", flush=True)
            print("   ", err.strip().splitlines()[-1], flush=True)
    if failed:
        print(f"\n[{label}] {len(failed)} job(s) failed:")
        for j in failed[:20]:
            print("   ", j)
    return results


# =====================================================================
# 5. Collection
# =====================================================================
def collect(task, results, args):
    rows = [dict(fold=r["fold"], n_train=r["n_train"], n_test=r["n_test"],
                 n_features=r.get("n_features"), test_auc=r["test_auc"])
            for r in results if r["task"] == task]
    if not rows:
        return None
    df = pd.DataFrame(rows).sort_values("fold")
    out = TASKS[task]["results"]
    if args.run_name:
        out = out.replace(".csv", f".{args.run_name}.csv")
    df.to_csv(out, index=False)
    print(f"\n=== {task} -> {out} ===")
    print(df.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print(f"  mean test_auc: {df.test_auc.mean():.4f}")
    return df


# =====================================================================
# 6. Main
# =====================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tasks", nargs="+", default=["disease"],
                    help=f"{list(TASKS)} or 'all'")
    ap.add_argument("--top-k", type=int, default=600,
                    help="keep only the K most abundant features per sample "
                         "(0 disables the filter)")
    ap.add_argument("--feature-space", default="union",
                    choices=["train", "union", "intersection"],
                    help="feature set of the merged table used for modelling")
    ap.add_argument("--renorm-after-topk", action="store_true",
                    help="re-scale each sample to sum 1 after the top-k cut "
                         "(off by default: normalisation happens before "
                         "selection, so the kept values are the original "
                         "full-table ranks)")
    ap.add_argument("--run-name", default="",
                    help="suffix added to result directories and CSV files")
    ap.add_argument("--overwrite", action="store_true",
                    help="re-run jobs even if result.json already exists")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--n-estimators", type=int, default=200)
    ap.add_argument("--max-depth", type=int, default=None)
    ap.add_argument("--n-jobs", type=int, default=10,
                    help="n_jobs passed to RandomForestClassifier")
    ap.add_argument("--seed", type=int, default=11)
    args = ap.parse_args()

    tasks = list(TASKS) if "all" in args.tasks else args.tasks
    for t in tasks:
        if t not in TASKS:
            raise SystemExit(f"unknown task {t}; choose from {list(TASKS)} or 'all'")
    if args.run_name:
        for t in tasks:
            TASKS[t]["artifact"] = f"{TASKS[t]['artifact']}_{args.run_name}"

    print(f"[cfg] tasks={tasks} top_k={args.top_k} "
          f"feature_space={args.feature_space}")
    print(f"[cfg] n_estimators={args.n_estimators} max_depth={args.max_depth} "
          f"n_jobs={args.n_jobs} seed={args.seed}")

    all_jobs = []
    for t in tasks:
        all_jobs += build_jobs(t, args)

    todo = all_jobs if args.overwrite else [
        j for j in all_jobs
        if not os.path.exists(os.path.join(j["out_dir"], "result.json"))]
    print(f"\n{len(all_jobs)} jobs, {len(all_jobs) - len(todo)} done, "
          f"{len(todo)} to run")

    if args.dry_run:
        for t in tasks:
            n = sum(1 for j in all_jobs if j["task"] == t)
            print(f"  {t:12s}: {n} folds")
        return

    run_jobs(todo)

    results = []
    for j in all_jobs:
        p = os.path.join(j["out_dir"], "result.json")
        if os.path.exists(p):
            with open(p) as f:
                results.append(json.load(f))
    print(f"\nDone: {len(results)}/{len(all_jobs)} results")

    for t in tasks:
        collect(t, results, args)


if __name__ == "__main__":
    main()