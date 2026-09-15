#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Sample efficiency of CRC: how much training data the test cohort needs.

The disease task's CRC folds, kept exactly as they are on the test side, with
the training side cut down. For every CRC row of
``run_leave_one_study_out_each_diease_list.tsv`` -- one held-out study S each --
the fold's own ``test_loo.biom`` is used untouched, while its training table is
reduced to 2, 3, 4 and 5 randomly drawn studies, 5 draws per size. The curve of
test AUC against the number of training cohorts is the answer.

Everything below the fold list is ``run_attention_biom_with_SNEs`` imported
unchanged -- members, splits, ensembling, GPU pool, results table -- so a point
on this curve is the same computation the full-data run does, only with fewer
training cohorts. Every flag of that script works here.

The draws
---------
A study is a candidate only if it holds at least --min-class of each class:
under --inner-split loso every drawn study is validated on exactly once, and a
single-class cohort has no AUC to select an epoch by. Draws of one size are
distinct; when a size has fewer than --reps distinct combinations, all of them
are run and the count is reported. The draw depends on --draw-seed alone.

With --inner-split loso a draw of k studies trains k members, each validating
on one of the k and training on the other k-1. At k=2 that is a member trained
on a single cohort, which is the point of the smallest size rather than a
problem with it.

Reads (the disease task's own data, nothing else):

  run_leave_one_study_out_each_diease_list.tsv     the CRC rows
  Data/disease_data/CRC/<study>/train_loo.biom     drawn from
  Data/disease_data/CRC/<study>/test_loo.biom      used as-is, never copied
  Data/disease_data/CRC/metadata.tsv

Writes:

  Data/disease_data/CRC_sample_efficiency/folds.tsv
  Data/disease_data/CRC_sample_efficiency/<study>_n<k>_r<i>/train_loo.biom
  crc_sample_efficiency.csv

The fold name carries the draw -- ``PRJEB6070_n3_r2`` is held-out study
PRJEB6070, three training cohorts, draw 2 -- and folds.tsv records which
studies each draw took, so the results CSV can be grouped by size afterwards.

Usage
-----
  python run_attention_biom_CRC_sample_efficiency.py --dry-run
  python run_attention_biom_CRC_sample_efficiency.py --gpus 0 1 2 3 4 5 6 7 \
      --run-name crc_se
  python run_attention_biom_CRC_sample_efficiency.py --self-check  # no Data/

Prep-only flags (everything else goes to run_attention_biom_with_SNEs):
  --data-root DIR   where the drawn training tables and folds.tsv are written
                    (default Data/disease_data/CRC_sample_efficiency)
  --sizes K [K ...] training-cohort counts to draw (default 2 3 4 5)
  --reps N          draws per size (default 5)
  --draw-seed N     seed for the draws (default 101)
  --min-class N     a study is a candidate only with >= N of each class
                    (default 1, the bare minimum for an AUC to exist)
  --held-out S [S ...]  restrict to these held-out studies (default: every CRC
                    row, which is len(sizes) x reps folds each)
  --rebuild         redraw even if folds.tsv is already there
  --self-check      run the draw asserts on synthetic data and exit
"""
from __future__ import annotations

import argparse
import itertools
import os
import sys

import biom
import numpy as np
import pandas as pd

import run_attention_biom_with_SNEs as R

DISEASE = "CRC"
TASK = "crc_sample_efficiency"


def crc_rows(run_tsv, held_out):
    """The disease task's CRC rows, as dicts, optionally restricted."""
    run = pd.read_csv(run_tsv, sep="\t")
    for col in ("disease", "study"):
        if col not in run.columns:
            raise SystemExit(f"{run_tsv} has no {col!r} column; found "
                             f"{list(run.columns)}")
    hit = run[run["disease"].astype(str) == DISEASE]
    if hit.empty:
        raise SystemExit(f"{run_tsv} has no {DISEASE} row; diseases present: "
                         f"{sorted(set(run['disease'].astype(str)))}")
    if held_out:
        want = [str(s) for s in held_out]
        hit = hit[hit["study"].astype(str).isin(want)]
        if hit.empty:
            raise SystemExit(
                f"--held-out {' '.join(want)} matches no {DISEASE} row; the "
                f"held-out studies are "
                f"{sorted(run[run.disease.astype(str) == DISEASE].study.astype(str))}")
    return [r.to_dict() for _, r in hit.iterrows()]


def candidate_studies(sids, md, study_col, min_class):
    """``(candidates, counts)`` -- the studies a draw may take, and why.

    A candidate holds at least `min_class` cases and `min_class` controls. Under
    --inner-split loso every drawn study is the validation set of one member, so
    a cohort short of a class cannot carry that member's epoch decision -- and
    at k=2 it would take its sibling down with it, since the other member would
    then train on it alone.
    """
    sids = np.asarray(sids)
    stu = R._column(md, sids, study_col, "study")
    lab = pd.to_numeric(md.loc[sids, R.LABELS_COL], errors="coerce")
    if lab.isna().any():
        raise SystemExit(f"non-numeric {R.LABELS_COL!r} for "
                         f"{int(lab.isna().sum())} sample(s); first offenders: "
                         f"{lab.index[lab.isna()][:5].tolist()}")
    lab = lab.to_numpy().astype(int)
    counts = {s: (int((lab[stu == s] == 1).sum()),
                  int((lab[stu == s] == 0).sum()))
              for s in sorted(set(stu))}
    cands = [s for s, (p, n) in counts.items()
             if p >= min_class and n >= min_class]
    return cands, counts


def draw_subsets(candidates, k, reps, rng):
    """Up to `reps` distinct k-subsets of `candidates`, as sorted tuples.

    Enumerating the combinations and drawing from them is what makes the draws
    distinct: with five candidates and k=4 there are five subsets in total, so
    five random draws with replacement would repeat two or three of them and
    spend GPU hours re-running the same training table. A size with fewer than
    `reps` combinations runs all of them instead, and the caller reports the
    shortfall.

    The candidate list of one disease is a handful of cohorts, so the
    enumeration is small; this would need rewriting for a pool of hundreds.
    """
    if k > len(candidates):
        return []
    combos = sorted(itertools.combinations(sorted(candidates), k))
    if len(combos) <= reps:
        return combos
    idx = rng.choice(len(combos), size=reps, replace=False)
    return [combos[i] for i in sorted(idx)]


def prepare(root, run_tsv, meta_tmpl, train_tmpl, study_col, sizes, reps,
            seed, min_class, held_out, rebuild):
    """Draw the reduced training tables. Returns ``(folds.tsv, metadata)``."""
    fold_tsv = os.path.join(root, "folds.tsv")
    rows = crc_rows(run_tsv, held_out)
    meta = meta_tmpl.format(**rows[0])
    metas = {meta_tmpl.format(**r) for r in rows}
    if len(metas) != 1:
        raise SystemExit(f"the {DISEASE} rows point at several metadata files "
                         f"{sorted(metas)}; this task assumes one")

    if os.path.exists(fold_tsv) and not rebuild:
        print(f"[prep ] {fold_tsv} is already there, reusing it and the tables "
              f"beside it (--rebuild to redraw)")
        return fold_tsv, meta

    md = pd.read_csv(meta, sep="\t", index_col=R.SAMPLE_ID_COL,
                     dtype={R.SAMPLE_ID_COL: str}, low_memory=False)
    os.makedirs(root, exist_ok=True)
    rng = np.random.default_rng(seed)
    out_rows = []

    for row in rows:
        S = str(row["study"])
        train_b = train_tmpl.format(**row)
        table = biom.load_table(train_b)
        sids = np.asarray(table.ids(axis="sample"))
        missing = [s for s in sids if s not in md.index]
        if missing:
            raise SystemExit(f"{len(missing)} samples of {train_b} have no row "
                             f"in {meta}; first offenders: {missing[:5]}")
        cands, counts = candidate_studies(sids, md, study_col, min_class)
        print(f"[prep ] held-out {S}: {len(counts)} training cohort(s), "
              f"{len(cands)} with >= {min_class} of each class "
              + "  ".join(f"{s} {p}/{n}" for s, (p, n) in counts.items()))
        if len(cands) < min(sizes):
            print(f"[warn ] held-out {S}: only {len(cands)} usable cohort(s), "
                  f"fewer than the smallest size {min(sizes)}; no fold")
            continue

        for k in sizes:
            picks = draw_subsets(cands, k, reps, rng)
            if not picks:
                print(f"[warn ] held-out {S}: k={k} needs {k} usable "
                      f"cohort(s), {len(cands)} available; skipped")
                continue
            if len(picks) < reps:
                print(f"[prep ] held-out {S}: k={k} has only {len(picks)} "
                      f"distinct combination(s) of {len(cands)} cohorts, so "
                      f"all of them are run instead of {reps} draws")
            for i, combo in enumerate(picks):
                fold = f"{S}_n{k}_r{i}"
                d = os.path.join(root, fold)
                os.makedirs(d, exist_ok=True)
                keep = sids[np.isin(R._column(md, sids, study_col, "study"),
                                    list(combo))]
                R._write_table(
                    table.filter(keep, axis="sample", inplace=False),
                    os.path.join(d, "train_loo.biom"))
                n_case = int((pd.to_numeric(md.loc[keep, R.LABELS_COL]) == 1).sum())
                out_rows.append(dict(fold=fold, study=S, n_studies=k, rep=i,
                                     studies=";".join(combo),
                                     n_train=len(keep), n_train_case=n_case,
                                     n_train_ctrl=len(keep) - n_case))

    if not out_rows:
        raise SystemExit("no draw was possible; see the warnings above")
    df = pd.DataFrame(out_rows)
    df.to_csv(fold_tsv, sep="\t", index=False)
    print(f"[prep ] {len(df)} fold(s) -> {fold_tsv}; "
          f"{df.n_studies.min()}..{df.n_studies.max()} training cohorts, "
          f"{int(df.n_train.min())}..{int(df.n_train.max())} training samples")
    print(f"[prep ] {df.n_studies.sum()} members in total under "
          f"--inner-split loso (one per drawn cohort)")
    return fold_tsv, meta


def self_check():
    """The draws: distinct, right size, exhaustive when they have to be."""
    md = pd.DataFrame(
        dict(study=["a", "a", "b", "b", "c", "c", "d", "d", "e", "e"],
             group=[1, 0, 1, 0, 1, 0, 1, 0, 1, 1]),
        index=[f"s{i}" for i in range(10)])
    cands, counts = candidate_studies(md.index, md, "study", 1)
    # 'e' is two cases and no control, so no member could be selected on it.
    assert cands == ["a", "b", "c", "d"], cands
    assert counts["e"] == (2, 0), counts
    assert candidate_studies(md.index, md, "study", 2)[0] == [], "min_class 2"

    rng = np.random.default_rng(0)
    d2 = draw_subsets(cands, 2, 5, rng)
    # C(4,2) = 6 > 5, so five distinct pairs are drawn.
    assert len(d2) == 5 and len(set(d2)) == 5, d2
    assert all(len(c) == 2 for c in d2), d2
    # C(4,4) = 1 < 5: run the one combination, do not repeat it five times.
    assert draw_subsets(cands, 4, 5, rng) == [("a", "b", "c", "d")]
    assert draw_subsets(cands, 5, 5, rng) == [], "k above the candidate count"
    # Same seed, same draws.
    assert draw_subsets(cands, 2, 5, np.random.default_rng(0)) == d2
    print("self-check ok")


def main():
    if {"-h", "--help"} & set(sys.argv[1:]):
        # Before prep, which would otherwise need Data/ just to print help.
        print(__doc__)
        sys.argv = [sys.argv[0], "--help"]
        R.main()
        return

    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--data-root",
                     default="Data/disease_data/CRC_sample_efficiency")
    pre.add_argument("--sizes", type=int, nargs="+", default=[2, 3, 4, 5])
    pre.add_argument("--reps", type=int, default=5)
    pre.add_argument("--draw-seed", type=int, default=101)
    pre.add_argument("--min-class", type=int, default=1)
    pre.add_argument("--held-out", nargs="+", default=[])
    pre.add_argument("--rebuild", action="store_true")
    pre.add_argument("--self-check", action="store_true")
    # Read here only to draw the same cohorts the run will hold out; passed on
    # to R.main() untouched.
    pre.add_argument("--inner-group", default="study")
    mine, rest = pre.parse_known_args()

    if mine.self_check:
        self_check()
        return
    if min(mine.sizes) < 2:
        raise SystemExit(f"--sizes must all be at least 2, got {mine.sizes}: "
                         f"--inner-split loso holds one drawn cohort out per "
                         f"member, so a draw of one leaves nothing to train on")
    if mine.reps < 1:
        raise SystemExit(f"--reps must be at least 1, got {mine.reps}")

    src = R.TASKS["disease"]
    fold_tsv, meta = prepare(mine.data_root, src["run_tsv"], src["meta"],
                             src["train"], mine.inner_group, mine.sizes,
                             mine.reps, mine.draw_seed, mine.min_class,
                             mine.held_out, mine.rebuild)

    R.TASKS[TASK] = dict(
        run_tsv=fold_tsv,
        train=os.path.join(mine.data_root, "{fold}", "train_loo.biom"),
        # The fold's own test table, untouched: only the training side is cut.
        test=src["test"].format(disease=DISEASE, study="{study}"),
        meta=meta,
        fold="{fold}",
        artifact=mine.data_root.rstrip("/") + "/",
        results=f"{TASK}.csv",
    )
    sys.argv = [sys.argv[0], "--tasks", TASK,
                "--inner-group", mine.inner_group] + rest
    R.main()


if __name__ == "__main__":
    main()

