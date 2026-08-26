#!/usr/bin/env python3
"""Create shuffled SNE or BIOM controls."""

import argparse
import random
from pathlib import Path


def shuffle_sne(input_path, output_path, seed):
    rng = random.Random(seed)
    with input_path.open() as source:
        rows = [line.split() for line in source if line.split()]
    if not rows or len(rows[0]) < 2:
        raise ValueError("SNE file contains no embedding vectors")
    if any(len(row) != len(rows[0]) for row in rows):
        raise ValueError("SNE vectors have inconsistent dimensions")

    for dimension in range(1, len(rows[0])):
        values = [row[dimension] for row in rows]
        rng.shuffle(values)
        for row, value in zip(rows, values):
            row[dimension] = value

    with output_path.open("x") as output:
        for row in rows:
            output.write(" ".join(row) + "\n")


def shuffle_biom(input_path, output_path, seed):
    import biom
    import numpy as np

    table = biom.load_table(str(input_path))
    data = table.matrix_data.toarray()
    rng = np.random.default_rng(seed)
    for sample in range(data.shape[1]):
        rng.shuffle(data[:, sample])

    shuffled = biom.Table(
        data, table.ids(axis="observation"), table.ids(axis="sample"),
        table.metadata(axis="observation"), table.metadata(axis="sample"),
        table_id=table.table_id, type=table.type)
    with biom.util.biom_open(str(output_path), "w") as handle:
        shuffled.to_hdf5(handle, f"within-sample shuffle; seed={seed}")


def shuffle_all(sne_input, biom_input, output_dir, seed):
    output_dir.mkdir()
    shuffled_sne = output_dir / "shuffled_sne.txt"
    shuffled_biom = output_dir / "shuffled_table.biom"
    shuffle_sne(sne_input, shuffled_sne, seed)
    shuffle_biom(biom_input, shuffled_biom, seed)

    conditions = [
        ("sne_only", shuffled_sne, biom_input),
        ("biom_only", sne_input, shuffled_biom),
        ("both", shuffled_sne, shuffled_biom),
    ]
    with (output_dir / "conditions.tsv").open("x") as output:
        output.write("condition\tsne\tbiom\tseed\n")
        for condition, sne, table in conditions:
            output.write(
                f"{condition}\t{sne.resolve()}\t{table.resolve()}\t{seed}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=("sne", "biom", "all"))
    parser.add_argument("paths", nargs="+", type=Path)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    expected = 3 if args.kind == "all" else 2
    if len(args.paths) != expected:
        parser.error(f"{args.kind} requires {expected} paths")

    if args.kind == "all":
        sne_input, biom_input, output_dir = args.paths
        if not sne_input.is_file() or not biom_input.is_file():
            parser.error("SNE or BIOM input file not found")
        if output_dir.exists():
            parser.error(f"output directory already exists: {output_dir}")
        shuffle_all(sne_input, biom_input, output_dir, args.seed)
        return

    input_path, output_path = args.paths
    if input_path.resolve() == output_path.resolve():
        parser.error("input and output must be different files")
    if not input_path.is_file():
        parser.error(f"input file not found: {input_path}")
    if output_path.exists():
        parser.error(f"output already exists: {output_path}")

    {"sne": shuffle_sne, "biom": shuffle_biom}[args.kind](
        input_path, output_path, args.seed)


if __name__ == "__main__":
    main()
