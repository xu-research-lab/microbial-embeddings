#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Shuffle each column of the social niche embedding matrix independently
(row labels stay in place), and save the result as a space-separated txt file.

Usage:
    python shuffle_columns.py \
        --input ../../data/social_niche_embedding_removing_disease_samples_100.txt \
        --output ../../data/social_niche_embedding_removing_disease_samples_100_shuffled.txt \
        --seed 11
"""

import argparse

import numpy as np
import pandas as pd


def shuffle_columns(df: pd.DataFrame, seed) -> pd.DataFrame:
    """Shuffle every column independently (within-column permutation)."""
    rng = np.random.default_rng(seed)
    values = df.to_numpy(copy=True)
    for j in range(values.shape[1]):
        rng.shuffle(values[:, j])          # in-place shuffle of column j
    return pd.DataFrame(values, index=df.index, columns=df.columns)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", "-i",
                   default="../../data/social_niche_embedding_removing_disease_samples_100.txt")
    p.add_argument("--output", "-o",
                   default="../../data/social_niche_embedding_removing_disease_samples_100_shuffled.txt")
    p.add_argument("--seed", "-s", type=int, default=5,
                   help="random seed for reproducibility")
    args = p.parse_args()

    # First column holds the node names (index); no header row.
    snes_embed = pd.read_csv(args.input, sep=" ", index_col=0, header=None)
    print(f"Loaded matrix: {snes_embed.shape[0]} rows x {snes_embed.shape[1]} columns")

    shuffled = shuffle_columns(snes_embed, seed=args.seed)

    # Save: space-separated, keep row names, no header.
    shuffled.to_csv(args.output, sep=" ", header=False, index=True)
    print(f"Saved to: {args.output}")

    # Sanity check: each column must contain the same values as before, only reordered.
    same_set = all(
        np.allclose(np.sort(snes_embed.iloc[:, j].to_numpy()),
                    np.sort(shuffled.iloc[:, j].to_numpy()))
        for j in range(snes_embed.shape[1])
    )
    print(f"Check - per-column value sets match the original: {same_set}")


if __name__ == "__main__":
    main()