#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Random Forest baseline for leave-one-continent-out CRC.

The comparison arm for ``run_attention_biom_CRC_continent.py``: same folds,
same tables, same metadata, a Random Forest instead of the attention model.

Everything below the fold list is ``run_rf.py`` imported unchanged -- the
rankdata/norm step, the top-k filter, the merge onto a shared feature space,
the RF itself, ``pred_test.csv``, ``result.json`` and the results CSV. Only the
task entry is new, so the two scripts' numbers are comparable by construction.

The settings are ``run_rf.py``'s own defaults, written out below rather than
parsed: top_k 600, feature_space union, 200 trees, no depth cap, seed 11. There
is nothing to pass on the command line.

Reads (all written by run_attention_biom_CRC_continent.py, which must have run):

  Data/disease_data/CRC_continent/folds.tsv
  Data/disease_data/CRC_continent/<continent>/{train,test}_loo.biom
  Data/disease_data/CRC/metadata.tsv

Writes:

  Data/disease_data/CRC_continent/<continent>/rf/{result.json,pred_test.csv}
  rf_crc_continent.csv

Usage
-----
  python run_rf_CRC_continent.py
"""
from __future__ import annotations

import argparse
import json
import os

import run_rf as RF

TASK = "crc_continent"
DATA_ROOT = "Data/disease_data/CRC_continent"
META = "Data/disease_data/CRC/metadata.tsv"

RF.TASKS[TASK] = dict(
    run_tsv=os.path.join(DATA_ROOT, "folds.tsv"),
    train=os.path.join(DATA_ROOT, "{continent}", "train_loo.biom"),
    test=os.path.join(DATA_ROOT, "{continent}", "test_loo.biom"),
    meta=META,
    fold="{continent}",
    artifact=DATA_ROOT.rstrip("/") + "/",
    results=f"rf_{TASK}.csv",
)

#: run_rf.py's defaults, one for one. build_jobs and collect read an argparse
#: Namespace, so that is what they get -- no parser, since nothing here varies.
ARGS = argparse.Namespace(
    tasks=[TASK], top_k=600, feature_space="union", renorm_after_topk=False,
    n_estimators=200, max_depth=None, n_jobs=10, seed=11,
    run_name="", overwrite=False, dry_run=False)


def main():
    t = RF.TASKS[TASK]
    if not os.path.exists(t["run_tsv"]):
        raise SystemExit(
            f"{t['run_tsv']} is not there. It is written by "
            f"run_attention_biom_CRC_continent.py, which builds the merged CRC "
            f"table and cuts it by continent; run that first")

    print(f"[cfg] task={TASK} top_k={ARGS.top_k} "
          f"feature_space={ARGS.feature_space}")
    print(f"[cfg] n_estimators={ARGS.n_estimators} max_depth={ARGS.max_depth} "
          f"n_jobs={ARGS.n_jobs} seed={ARGS.seed}")

    jobs = RF.build_jobs(TASK, ARGS)
    todo = [j for j in jobs
            if not os.path.exists(os.path.join(j["out_dir"], "result.json"))]
    print(f"\n{len(jobs)} jobs, {len(jobs) - len(todo)} done, "
          f"{len(todo)} to run")

    RF.run_jobs(todo)

    results = []
    for j in jobs:
        p = os.path.join(j["out_dir"], "result.json")
        if os.path.exists(p):
            with open(p) as f:
                results.append(json.load(f))
    print(f"\nDone: {len(results)}/{len(jobs)} results")
    RF.collect(TASK, results, ARGS)


if __name__ == "__main__":
    main()

