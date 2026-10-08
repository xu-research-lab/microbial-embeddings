#!/usr/bin/env python
"""Fig. 2A input: how strongly each taxon points along the "oxygen direction".

For every (facultative anaerobe i, obligate anaerobe j) pair the difference
d_ij = v_i - v_j of their SNEs is one estimate of the direction that separates
the two lifestyles (the "king - man + woman" arithmetic of word embeddings).
Each taxon t is scored by its mean cosine to all of them,

    s_t = mean_ij cos(v_t, d_ij),

which equals cos-to-|v_t| of the mean *unit* difference vector, so the
169 x 655 = 110,695 differences never have to be scored one by one.

Oxygen groups come from Traitar's three oxygen phenotypes in
data/trait_predcit.csv (Anaerobe, Aerobe, Facultative), where 3 means both of
Traitar's classifiers call the phenotype and 0 that neither does. A taxon is
kept only when exactly one of the three is called:

    Anaerobe    == 3, the other two 0  -> anaerobic       (655 taxa)
    Aerobe      == 3, the other two 0  -> aerobic         (144 taxa)
    Facultative == 3, the other two 0  -> facultatively   (169 taxa)

so a taxon Traitar calls both anaerobic and facultative is left out.

Run from analysis/traits:  python oxygen_vector.py
Writes data/Traitar_Facultatively_Anaerobic_Anaerobic_Aerobic_all_co.csv, the
file traits_results.ipynb reads for Fig. 2A. `--check` compares against the
committed copy instead of overwriting it.
"""
import argparse
import sys

import numpy as np
import pandas as pd

EMB = "../../data/social_niche_embedding_100.txt"
TRAITS = "data/trait_predcit.csv"
OUT = "data/Traitar_Facultatively_Anaerobic_Anaerobic_Aerobic_all_co.csv"


def oxygen_scores(emb_path=EMB, traits_path=TRAITS):
    emb = pd.read_csv(emb_path, sep=" ", header=None, index_col=0)
    emb = emb.drop(index="<unk>", errors="ignore")
    tr = pd.read_csv(traits_path, index_col=0)
    tr = tr[tr.index.isin(emb.index)]

    an, ae, fa = tr["Anaerobe"], tr["Aerobe"], tr["Facultative"]
    group = pd.Series(np.select([(an == 3) & (ae == 0) & (fa == 0),
                                 (ae == 3) & (an == 0) & (fa == 0),
                                 (fa == 3) & (an == 0) & (ae == 0)],
                                ["anaerobic", "aerobic", "facultatively"], default=""),
                      index=tr.index)
    group = group[group != ""]

    fac = emb.loc[group.index[group == "facultatively"]].to_numpy()
    obl = emb.loc[group.index[group == "anaerobic"]].to_numpy()
    d = fac[:, None, :] - obl[None, :, :]
    u = (d / np.linalg.norm(d, axis=2, keepdims=True)).mean(axis=(0, 1))

    v = emb.loc[group.index].to_numpy()
    cosine = v @ u / np.linalg.norm(v, axis=1)
    return pd.DataFrame({"fid": group.index, "cosine": cosine, "Type": group.to_numpy()})


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true",
                    help="compare with the committed file instead of overwriting it")
    args = ap.parse_args()

    res = oxygen_scores()
    print("taxa per group:", res.Type.value_counts().to_dict())
    if not args.check:
        res.to_csv(OUT, index=False)
        print("wrote", OUT)
        return

    ref = pd.read_csv(OUT)
    same_rows = res[["fid", "Type"]].equals(ref[["fid", "Type"]])
    err = np.abs(res.cosine.to_numpy() - ref.cosine.to_numpy()).max() if same_rows else np.inf
    print("same taxa, order and groups: %s; max |cosine difference| %.1e" % (same_rows, err))
    sys.exit(0 if same_rows and err < 1e-12 else 1)


if __name__ == "__main__":
    main()
