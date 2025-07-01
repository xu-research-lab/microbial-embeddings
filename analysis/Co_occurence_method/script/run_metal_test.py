import os
import re
import glob
import pandas as pd
import numpy as np
from skbio.stats.distance import mantel
from itertools import combinations
from concurrent.futures import ProcessPoolExecutor
from functools import partial

# --- CONFIGURATION ---
MATRIX_FOLDER = '/home/cjj/projects/membed_local/membed/analysis_distance_true_datasize/similarity_matrix'
OUTPUT_CSV = '/home/cjj/projects/membed_local/membed/analysis_distance_true_datasize/mantel_test_results.csv'
NUM_CORES_TO_USE = None # Set to an integer for a specific number of cores, or None to use all available
MANTEL_PERMUTATIONS = 999


def parse_filename(filepath):
    """Parses the filepath to extract the metric name and a custom label."""
    filename = os.path.basename(filepath)
    # Match for files like similarity_abundance_percentile_X_100.csv
    match_numbered = re.search(r'similarity_([\w_]+?)_(\d+)_100\.csv', filename)
    if match_numbered:
        metric_name = match_numbered.group(1).replace('_', ' ').strip()
        label_num = match_numbered.group(2)
        identifier = f"{metric_name}-{label_num}"
        return metric_name, identifier, filepath
    # Match for files like similarity_abundance_percentile_100.csv (the baseline)
    match_baseline = re.search(r'similarity_([\w_]+?)_100\.csv', filename)
    if match_baseline:
        metric_name = match_baseline.group(1).replace('_', ' ').strip()
        identifier = f"{metric_name}-baseline" # Use a distinct identifier for the baseline
        return metric_name, identifier, filepath
    return None, None, filepath

def run_mantel_for_pair(pair_of_files):
    """
    This function is the core workload for each process.
    It takes a tuple of two file paths, performs the Mantel test, and returns the result.
    """
    filepath1, filepath2 = pair_of_files
    try:
        df1 = pd.read_csv(filepath1, index_col=0)
        df2 = pd.read_csv(filepath2, index_col=0)

        common_labels = df1.index.intersection(df2.index)
        if len(common_labels) < 3: # Mantel test requires at least 3 common observations
            return None
        
        df1_aligned = df1.loc[common_labels, common_labels]
        df2_aligned = df2.loc[common_labels, common_labels]

        # Convert similarity to distance: distance = 1 - similarity
        dist_matrix1 = 1 - df1_aligned.to_numpy()
        dist_matrix2 = 1 - df2_aligned.to_numpy()
        
        # Enforce hollow property: diagonal elements must be zero
        np.fill_diagonal(dist_matrix1, 0)
        np.fill_diagonal(dist_matrix2, 0)

        r_val, p_val, _ = mantel(dist_matrix1, dist_matrix2, permutations=MANTEL_PERMUTATIONS)
        
        _, id1, _ = parse_filename(filepath1)
        _, id2, _ = parse_filename(filepath2)
        
        return (id1, id2, r_val, p_val)
    except Exception as e:
        print(f"Error processing pair ({os.path.basename(filepath1)}, {os.path.basename(filepath2)}): {e}")
        return None


def main():
    """Main function to orchestrate the whole process."""
    print(f"Scanning for matrices in: {MATRIX_FOLDER}")
    all_files = glob.glob(os.path.join(MATRIX_FOLDER, 'similarity_*.csv'))

    grouped_files = {}
    baseline_files = {} # To store the 'percentile_100' files separately

    for f in all_files:
        metric, identifier, _ = parse_filename(f)
        if metric:
            # Check if it's the baseline file (e.g., 'similarity_abundance_percentile_100.csv')
            if 'baseline' in identifier:
                baseline_files[metric] = f
            else:
                if metric not in grouped_files:
                    grouped_files[metric] = []
                grouped_files[metric].append(f)

    print(f"\nFound {len(grouped_files)} groups of matrices and {len(baseline_files)} baseline matrices.")

    all_pairs_to_test = []
    for metric, files in grouped_files.items():
        if metric in baseline_files:
            baseline_file = baseline_files[metric]
            for file_to_test in files:
                # Add pairs where one is the baseline and the other is a numbered percentile file
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