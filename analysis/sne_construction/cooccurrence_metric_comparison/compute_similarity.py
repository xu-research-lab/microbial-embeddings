import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = (
    SCRIPT_DIR.parent
    / "data/cooccurrence_metric_comparison/difference_study_training/embeding_list"
)
DEFAULT_OUTPUT = SCRIPT_DIR / "results/similarity_matrices"


def compute_similarity(embedding_file, output_dir):
    embedding_file = Path(embedding_file)
    table = pd.read_csv(embedding_file, sep=r"\s+", header=None).dropna()
    feature_ids = table.iloc[:, 0]
    similarities = cosine_similarity(table.iloc[:, 1:])
    output_file = output_dir / f"similarity_{embedding_file.stem}.csv"
    pd.DataFrame(similarities, index=feature_ids, columns=feature_ids).to_csv(output_file)
    return output_file


def main():
    parser = argparse.ArgumentParser(
        description="Compute dense cosine-similarity matrices for embedding files."
    )
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()

    embedding_files = sorted(args.input_dir.glob("*.txt"))
    if not embedding_files:
        parser.error(f"no embedding files found in {args.input_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        outputs = executor.map(
            compute_similarity,
            embedding_files,
            [args.output_dir] * len(embedding_files),
        )
        for output in outputs:
            print(output)


if __name__ == "__main__":
    main()
