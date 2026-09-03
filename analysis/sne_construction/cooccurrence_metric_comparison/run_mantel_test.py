import argparse
import re
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from pathlib import Path

import numpy as np
import pandas as pd
from skbio.stats.distance import mantel


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = SCRIPT_DIR / "results/similarity_matrices"
DEFAULT_OUTPUT = SCRIPT_DIR / "results/mantel_test_results.csv"
NAME_PATTERN = re.compile(r"similarity_(.+?)(?:_(\d+))?_100\.csv$")


def describe_matrix(path):
    match = NAME_PATTERN.fullmatch(path.name)
    if not match:
        return None
    metric, study = match.groups()
    study = int(study) if study is not None else 0
    return metric, study, f"{metric.replace('_', ' ')}-{study}", path


def compare_pair(pair, permutations):
    left, right = pair
    left_table = pd.read_csv(left[3], index_col=0)
    right_table = pd.read_csv(right[3], index_col=0)
    common = left_table.index.intersection(right_table.index)
    if len(common) < 3:
        return None

    left_distance = 1 - left_table.loc[common, common].to_numpy()
    right_distance = 1 - right_table.loc[common, common].to_numpy()
    np.fill_diagonal(left_distance, 0)
    np.fill_diagonal(right_distance, 0)
    r_value, p_value, _ = mantel(
        left_distance, right_distance, permutations=permutations
    )
    return left[2], right[2], r_value, p_value


def main():
    parser = argparse.ArgumentParser(
        description="Compare each study embedding with its full-cohort baseline."
    )
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()

    matrices = [
        item
        for path in sorted(args.input_dir.glob("similarity_*.csv"))
        if (item := describe_matrix(path)) is not None
    ]
    baselines = {item[0]: item for item in matrices if item[1] == 0}
    pairs = [
        (baselines[item[0]], item)
        for item in matrices
        if item[1] != 0 and item[0] in baselines
    ]
    if not pairs:
        parser.error(f"no baseline/study pairs found in {args.input_dir}")

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        results = executor.map(compare_pair, pairs, repeat(args.permutations))
        rows = [result for result in results if result is not None]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        rows, columns=["Matrix 1", "Matrix 2", "R_value", "p_value"]
    ).to_csv(args.output, index=False)
    print(args.output)


if __name__ == "__main__":
    main()
