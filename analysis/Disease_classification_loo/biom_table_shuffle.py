#!/usr/bin/env python3
"""Create shuffled SNE or BIOM controls."""

import pandas as pd
import argparse
import random
from pathlib import Path

def shuffle_biom(input_path, output_path, seed):
    import biom
    import numpy as np

    table = biom.load_table(str(input_path))
    data = table.matrix_data.toarray()
    rng = np.random.default_rng(seed)

    shape = data.shape
    data = rng.permutation(data.ravel()).reshape(shape)

    shuffled = biom.Table(
        data, table.ids(axis="observation"), table.ids(axis="sample"),
        table.metadata(axis="observation"), table.metadata(axis="sample"),
        table_id=table.table_id, type=table.type)
    with biom.util.biom_open(str(output_path), "w") as handle:
        shuffled.to_hdf5(handle, f"global shuffle (all values); seed={seed}")


def main():
    run_list = pd.read_csv("run_shuffled_table.tsv", sep="\t")
    for i in range(run_list.shape[0]):
        disease = run_list.iloc[i, 0]
        study = run_list.iloc[i, 1]
        input_path = f"Data/disease_data/{disease}/{study}/train_loo.biom"
        output_path = f"Data/shuffle_table_IBD_CRC/{disease}/{study}/train_loo.biom"
        # if disease == "IBD":
        shuffle_biom(input_path, output_path, seed=42)

        input_path = f"Data/disease_data/{disease}/{study}/test_loo.biom"
        output_path = f"Data/shuffle_table_IBD_CRC/{disease}/{study}/test_loo.biom"
        # if disease == "IBD":
        shuffle_biom(input_path, output_path, seed=42)




if __name__ == "__main__":
    main()
