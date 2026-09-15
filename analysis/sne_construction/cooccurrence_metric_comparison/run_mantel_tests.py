import os
import re
import glob
import pandas as pd
import numpy as np
from skbio.stats.distance import mantel
from concurrent.futures import ProcessPoolExecutor

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MATRIX_FOLDER = os.path.join(SCRIPT_DIR, "results/similarity_matrix")
OUTPUT_CSV = os.path.join(SCRIPT_DIR, "results/mantel_test_results.csv")
NUM_CORES_TO_USE = None
MANTEL_PERMUTATIONS = 999
NAME_PATTERN = re.compile(r"similarity_([\w_]+?)(?:_(\d+))?_100\.csv$")


def parse_filename(filepath):
    """Return the metric name, study index, label, and path."""
    match = NAME_PATTERN.fullmatch(os.path.basename(filepath))
    if not match:
        return None

    metric_name, study_index = match.groups()
    study_index = int(study_index) if study_index is not None else 0
    metric_name = metric_name.replace("_", " ").strip()
    return metric_name, study_index, f"{metric_name}-{study_index}", filepath


def run_mantel_for_pair(pair_of_files):
    """Run one Mantel test for a pair of similarity matrices."""
    filepath1, filepath2 = pair_of_files
    try:
        df1 = pd.read_csv(filepath1, index_col=0)
        df2 = pd.read_csv(filepath2, index_col=0)

        common_labels = df1.index.intersection(df2.index)
        if len(common_labels) < 3:  # Mantel requires at least three observations.
            return None
        
        df1_aligned = df1.loc[common_labels, common_labels]
        df2_aligned = df2.loc[common_labels, common_labels]

        # Mantel expects distances rather than similarities.
        dist_matrix1 = 1 - df1_aligned.to_numpy()
        dist_matrix2 = 1 - df2_aligned.to_numpy()
        
        # Set both distance-matrix diagonals to zero.
        np.fill_diagonal(dist_matrix1, 0)
        np.fill_diagonal(dist_matrix2, 0)

        r_val, p_val, _ = mantel(dist_matrix1, dist_matrix2, permutations=MANTEL_PERMUTATIONS)
        
        _, _, id1, _ = parse_filename(filepath1)
        _, _, id2, _ = parse_filename(filepath2)
        
        return (id1, id2, r_val, p_val)
    except Exception as e:
        print(f"Error processing pair ({os.path.basename(filepath1)}, {os.path.basename(filepath2)}): {e}")
        return None


def main():
    """Run Mantel tests against each metric's baseline matrix."""
    print(f"Scanning for matrices in: {MATRIX_FOLDER}")
    all_files = sorted(glob.glob(os.path.join(MATRIX_FOLDER, 'similarity_*.csv')))

    grouped_files = {}
    baseline_files = {}

    for filepath in all_files:
        parsed = parse_filename(filepath)
        if parsed is None:
            continue

        metric, study_index, _, path = parsed
        if study_index == 0:
            baseline_files[metric] = path
        else:
            grouped_files.setdefault(metric, []).append(path)

    print(f"\nFound {len(grouped_files)} groups of matrices and {len(baseline_files)} baseline matrices.")

    all_pairs_to_test = []
    for metric, files in sorted(grouped_files.items()):
        if metric in baseline_files:
            baseline_file = baseline_files[metric]
            for file_to_test in files:
                all_pairs_to_test.append((baseline_file, file_to_test))
        else:
            print(f"Warning: No baseline file found for metric '{metric}'. Skipping tests for this metric.")

    if not all_pairs_to_test:
        print("Could not generate any pairs for testing. Please check file names and counts.")
        return

    print(f"Generated a total of {len(all_pairs_to_test)} pairs for Mantel testing (relative to baseline).")
    print(f"Using {MANTEL_PERMUTATIONS} permutations for each test.")
    if NUM_CORES_TO_USE:
        print(f"Running on {NUM_CORES_TO_USE} CPU cores...")
    else:
        print("Running on all available CPU cores...")
    
    with ProcessPoolExecutor(max_workers=NUM_CORES_TO_USE) as executor:
        results = list(executor.map(run_mantel_for_pair, all_pairs_to_test))

    valid_results = [res for res in results if res is not None]
    
    if not valid_results:
        print("\nNo valid results were generated from the Mantel tests. Please check error messages above.")
        return

    result_df = pd.DataFrame(valid_results, columns=['Matrix 1', 'Matrix 2', 'R_value', 'p_value'])
    result_df.to_csv(OUTPUT_CSV, index=False)
    
    print(f"\nProcessing complete! All results saved to: {OUTPUT_CSV}")

if __name__ == "__main__":
    main()
