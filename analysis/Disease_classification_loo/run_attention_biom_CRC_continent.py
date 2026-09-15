#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Leave-one-continent-out CRC, on top of run_attention_biom_with_SNEs.

The question is geographic transfer rather than cohort transfer: one continent
is the test set and every other continent trains. Everything below the fold
split -- members, splits, ensembling, GPU pool, results table -- is imported
from ``run_attention_biom_with_SNEs`` unchanged, so this file only has to build
the folds and hand them over. Every flag of that script therefore works here.

The folds
---------
CRC lives in nine countries (Argentina, Chile, UK, India, Vietnam, China,
Japan, France, USA), which fall into four of the seven continents:

    South America  Argentina, Chile
    Europe         UK, France
    Asia           India, Vietnam, China, Japan
    North America  USA

Controls may add more countries, so `CONTINENTS` below covers rather more than
those nine; an unmapped country stops the run with its name rather than being
guessed at or quietly dropped.

The table
---------
There is no single CRC BIOM on disk -- the disease task stores one train/test
pair per held-out study. The first CRC row of
``run_leave_one_study_out_each_diease_list.tsv`` is read and its two tables
concatenated, which by construction is every CRC sample exactly once. That is
the only file this script needs from you; everything below it is written here:

  Data/disease_data/CRC_continent/CRC_all.biom          the merged CRC table
  Data/disease_data/CRC_continent/folds.tsv             the continent list
  Data/disease_data/CRC_continent/<continent>/train_loo.biom
  Data/disease_data/CRC_continent/<continent>/test_loo.biom

The continents are worked out here from the metadata's country column -- there
is no fold list to supply. `folds.tsv` is an output, written so a run carries
the folds it was made from; --rebuild redoes all of it. The metadata stays the
disease task's own CRC metadata.tsv, which already covers every merged sample.

Validation
----------
The disease task's protocol, unchanged: --inner-split loso by default, one
member per training cohort (--inner-group, default 'study'), each validating on
the cohort it held out. The epoch is therefore still selected across studies,
while the test asks for a step across continents.

Usage
-----
  python run_attention_biom_CRC_continent.py --dry-run
  python run_attention_biom_CRC_continent.py --gpus 0 1 2 3 4 5 6 7 --run-name crc_geo
  python run_attention_biom_CRC_continent.py --self-check   # no Data/ needed

Prep-only flags (everything else goes to run_attention_biom_with_SNEs):
  --data-root DIR   where the merged table, the fold list and the
                    per-continent tables are written
                    (default Data/disease_data/CRC_continent)
  --rebuild         rebuild all of that even if it is already there
  --self-check      run the fold-planning asserts on synthetic data and exit
"""
from __future__ import annotations

import argparse
import os
import sys

import biom
import numpy as np
import pandas as pd

import run_attention_biom_with_SNEs as R

DISEASE = "CRC"
TASK = "crc_continent"

#: Continent -> the country spellings that land in it. Aliases are extra
#: entries rather than a second table, and the lookup is case-insensitive.
CONTINENTS = {
    "Europe": ["UK", "United Kingdom", "Great Britain", "England", "France",
               "Germany", "Austria", "Denmark", "Spain", "Italy", "Sweden",
               "Norway", "Finland", "Netherlands", "Belgium", "Switzerland",
               "Ireland", "Poland", "Hungary", "Estonia", "Latvia", "Iceland",
               "Luxembourg", "Portugal", "Czech Republic", "Slovenia",
               "Croatia", "Greece", "Russia"],
    "Asia": ["China", "Japan", "India", "Vietnam", "Viet Nam", "Singapore",
             "South Korea", "Korea", "Republic of Korea", "Israel",
             "Kazakhstan", "Mongolia", "Bangladesh", "Indonesia", "Malaysia",
             "Thailand", "Philippines", "Pakistan", "Nepal", "Sri Lanka",
             "Taiwan", "Hong Kong", "Turkey", "Iran", "Saudi Arabia"],
    "North_America": ["USA", "United States", "United States of America", "US",
                      "Canada", "Mexico", "Costa Rica", "Cuba"],
    "South_America": ["Argentina", "Chile", "Brazil", "Peru", "Colombia",
                      "Venezuela", "Ecuador", "Bolivia", "Uruguay"],
    "Africa": ["Tanzania", "Ghana", "Madagascar", "Egypt", "South Africa",
               "Nigeria", "Ethiopia", "Burkina Faso", "Cameroon", "Kenya",
               "Malawi", "Morocco", "Gambia", "Mali", "Botswana"],
    "Oceania": ["Australia", "New Zealand", "Fiji", "Papua New Guinea"],
    # Antarctica has no cohort and never will; it is here so the seven are all
    # named and an entry has somewhere obvious to go if one ever appears.
    "Antarctica": [],
}
_LOOKUP = {c.strip().lower(): cont
           for cont, countries in CONTINENTS.items() for c in countries}


def continent_of(country):
    """Continent of one country name. Raises KeyError for an unmapped one."""
    key = str(country).strip().lower()
    if key not in _LOOKUP:
        raise KeyError(country)
    return _LOOKUP[key]


def crc_row(run_tsv):
    """The first CRC row of the disease task's fold list, as a dict.

    Its train/test pair is what gets merged back into the whole CRC table. Any
    CRC row would do -- each is the same samples cut a different way -- so the
    first is taken and printed, since which one it was is otherwise invisible
    afterwards.
    """
    run = pd.read_csv(run_tsv, sep="\t")
    if "disease" not in run.columns:
        raise SystemExit(f"{run_tsv} has no 'disease' column; found "
                         f"{list(run.columns)}")
    hit = run[run["disease"].astype(str) == DISEASE]
    if hit.empty:
        raise SystemExit(f"{run_tsv} has no {DISEASE} row; diseases present: "
                         f"{sorted(set(run['disease'].astype(str)))}")
    return hit.iloc[0].to_dict()


def plan_folds(sids, md, country_col, study_col):
    """``(folds, skipped)`` -- which samples test and train under each continent.

    Parameters
    ----------
    sids : array of str
        Every sample of the merged CRC table.
    md : pandas.DataFrame
        Metadata indexed by sample id, carrying `country_col`, `study_col` and
        `R.LABELS_COL`.

    Returns
    -------
    folds : dict
        ``{continent: {'test': ids, 'train': ids, 'countries': [...],
        'n_train_studies': int}}``, in sorted continent order.
    skipped : dict
        ``{continent: reason}`` for the continents that cannot be a fold. A
        side with one class has no AUC, and a training side with a single
        cohort has nothing for --inner-split loso to validate on.
    """
    sids = np.asarray(sids)
    cty = R._column(md, sids, country_col, "country")
    unmapped = sorted({c for c in cty if c.strip().lower() not in _LOOKUP})
    if unmapped:
        raise SystemExit(
            f"no continent for {', '.join(repr(c) for c in unmapped)}; add "
            f"each to CONTINENTS in {os.path.basename(__file__)} -- one line, "
            f"and better than this script guessing")
    cont = np.array([continent_of(c) for c in cty])
    lab = pd.to_numeric(md.loc[sids, R.LABELS_COL], errors="coerce")
    if lab.isna().any():
        raise SystemExit(f"non-numeric {R.LABELS_COL!r} for "
                         f"{int(lab.isna().sum())} sample(s); first offenders: "
                         f"{lab.index[lab.isna()][:5].tolist()}")
    lab = lab.to_numpy().astype(int)
    stu = R._column(md, sids, study_col, "study")

    folds, skipped = {}, {}
    for c in sorted(set(cont)):
        is_test = cont == c
        n_cls = (len(np.unique(lab[is_test])), len(np.unique(lab[~is_test])))
        if n_cls[0] < 2:
            skipped[c] = (f"its {int(is_test.sum())} sample(s) are all one "
                          f"class, so a test AUC does not exist")
            continue
        if n_cls[1] < 2:
            skipped[c] = ("holding it out leaves the training side with one "
                          "class")
            continue
        n_stud = len(set(stu[~is_test]))
        if n_stud < 2:
            skipped[c] = (f"the training side is a single cohort ({n_stud}), "
                          f"so --inner-split loso has nothing to hold out")
            continue
        folds[c] = dict(test=sids[is_test], train=sids[~is_test],
                        countries=sorted(set(cty[is_test])),
                        n_train_studies=n_stud,
                        n_test=int(is_test.sum()),
                        n_test_case=int((lab[is_test] == 1).sum()),
                        n_test_ctrl=int((lab[is_test] == 0).sum()),
                        n_train=int((~is_test).sum()),
                        n_train_case=int((lab[~is_test] == 1).sum()),
                        n_train_ctrl=int((lab[~is_test] == 0).sum()))
    return folds, skipped


def prepare(root, run_tsv, meta_tmpl, train_tmpl, test_tmpl, country_col,
            study_col, rebuild):
    """Merge the CRC table and cut it by continent. Returns ``(tsv, metadata)``.

    Both are written under `root`, along with the per-continent tables. Nothing
    here is an input you have to have: the only file read is the disease task's
    own ``run_leave_one_study_out_each_diease_list.tsv``.
    """
    fold_tsv = os.path.join(root, "folds.tsv")
    merged = os.path.join(root, f"{DISEASE}_all.biom")
    row = crc_row(run_tsv)
    meta = meta_tmpl.format(**row)

    if os.path.exists(fold_tsv) and not rebuild:
        print(f"[prep ] {fold_tsv} is already there, reusing it and the "
              f"tables beside it (--rebuild to redo the split)")
        return fold_tsv, meta

    os.makedirs(root, exist_ok=True)
    if os.path.exists(merged) and not rebuild:
        print(f"[prep ] reusing the merged {DISEASE} table {merged}")
        table = biom.load_table(merged)
    else:
        train_b, test_b = train_tmpl.format(**row), test_tmpl.format(**row)
        print(f"[prep ] merging {DISEASE} from the first fold of {run_tsv} "
              f"(study {row.get('study')!r}): {train_b} + {test_b}")
        # Disjoint samples, union of taxa: the fold's two halves are every CRC
        # sample exactly once, so this is the whole disease back in one table.
        table = biom.load_table(train_b).concat([biom.load_table(test_b)],
                                                axis="sample")
        R._write_table(table, merged)
        print(f"[prep ] wrote {merged}")
    sids = np.asarray(table.ids(axis="sample"))
    md = pd.read_csv(meta, sep="\t", index_col=R.SAMPLE_ID_COL,
                     dtype={R.SAMPLE_ID_COL: str}, low_memory=False)
    missing = [s for s in sids if s not in md.index]
    if missing:
        raise SystemExit(f"{len(missing)} merged samples have no row in "
                         f"{meta}; first offenders: {missing[:5]}")
    print(f"[prep ] merged table: {len(sids)} samples x "
          f"{len(table.ids(axis='observation'))} taxa")

    folds, skipped = plan_folds(sids, md, country_col, study_col)
    for c, why in sorted(skipped.items()):
        print(f"[warn ] no fold for {c}: {why}")
    if not folds:
        raise SystemExit("no continent can be a fold; see the warnings above")

    for c, f in folds.items():
        d = os.path.join(root, c)
        os.makedirs(d, exist_ok=True)
        R._write_table(table.filter(f["train"], axis="sample", inplace=False),
                       os.path.join(d, "train_loo.biom"))
        R._write_table(table.filter(f["test"], axis="sample", inplace=False),
                       os.path.join(d, "test_loo.biom"))
        print(f"[prep ] {c:<14} test {f['n_test']:>4} "
              f"({f['n_test_case']}/{f['n_test_ctrl']} case/control) from "
              f"{', '.join(f['countries'])}  |  train {f['n_train']:>4} "
              f"({f['n_train_case']}/{f['n_train_ctrl']}) over "
              f"{f['n_train_studies']} cohort(s) -> {f['n_train_studies']} "
              f"members")
    pd.DataFrame({"continent": sorted(folds)}).to_csv(fold_tsv, sep="\t",
                                                      index=False)
    print(f"[prep ] {len(folds)} fold(s) -> {fold_tsv} (written by this "
          f"script; there is no continent fold list to supply)")
    return fold_tsv, meta


def self_check():
    """Fold planning on synthetic data: mapping, skipping, and who trains."""
    sids = [f"s{i}" for i in range(10)]
    md = pd.DataFrame(
        dict(country=["UK", "France", "China", "Japan", "USA", "USA",
                      "Chile", "Argentina", "Ghana", "Ghana"],
             study=["e1", "e2", "a1", "a2", "n1", "n1", "s1", "s2", "f1", "f1"],
             group=[1, 0, 1, 0, 1, 0, 1, 0, 1, 1]), index=sids)
    folds, skipped = plan_folds(sids, md, "country", "study")

    assert set(folds) == {"Europe", "Asia", "North_America", "South_America"}, folds
    # Africa is one class, so it is not a fold -- but its samples still train.
    assert "Africa" in skipped and "one class" in skipped["Africa"], skipped
    assert sorted(folds["Europe"]["test"]) == ["s0", "s1"]
    assert "s8" in folds["Europe"]["train"], "a skipped continent still trains"
    assert len(folds["Europe"]["train"]) == 8
    assert folds["Asia"]["countries"] == ["China", "Japan"]
    # USA is one cohort, so holding it out leaves the other 7 to train on.
    assert folds["North_America"]["n_train_studies"] == 7, folds["North_America"]
    assert continent_of("united states") == "North_America"
    try:
        continent_of("Atlantis")
    except KeyError:
        pass
    else:
        raise AssertionError("an unmapped country must not resolve")
    print("self-check ok")


def main():
    if {"-h", "--help"} & set(sys.argv[1:]):
        # Before prep, which would otherwise need Data/ just to print help.
        print(__doc__)
        sys.argv = [sys.argv[0], "--help"]
        R.main()
        return

    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--data-root", default="Data/disease_data/CRC_continent")
    pre.add_argument("--rebuild", action="store_true")
    pre.add_argument("--self-check", action="store_true")
    # Read here only to split the table the same way the run will; both are
    # passed on to R.main() untouched.
    pre.add_argument("--country-col", default="country")
    pre.add_argument("--inner-group", default="study")
    mine, rest = pre.parse_known_args()

    if mine.self_check:
        self_check()
        return

    src = R.TASKS["disease"]
    fold_tsv, meta = prepare(mine.data_root, src["run_tsv"], src["meta"],
                             src["train"], src["test"], mine.country_col,
                             mine.inner_group, mine.rebuild)

    R.TASKS[TASK] = dict(
        run_tsv=fold_tsv,
        train=os.path.join(mine.data_root, "{continent}", "train_loo.biom"),
        test=os.path.join(mine.data_root, "{continent}", "test_loo.biom"),
        meta=meta,
        fold="{continent}",
        artifact=mine.data_root.rstrip("/") + "/",
        results=f"{TASK}.csv",
    )
    # --inner-group is consumed by the pre-parser above, so hand it back.
    sys.argv = [sys.argv[0], "--tasks", TASK,
                "--inner-group", mine.inner_group] + rest
    R.main()


if __name__ == "__main__":
    main()
