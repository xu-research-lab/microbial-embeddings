import os
import glob
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
from concurrent.futures import ProcessPoolExecutor
from functools import partial

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ANALYSIS_DIR = os.path.dirname(SCRIPT_DIR)
INPUT_FOLDER = os.path.join(
    ANALYSIS_DIR,
    "data/cooccurrence_metric_comparison/difference_study_training/embeding_list",
)
OUTPUT_FOLDER = os.path.join(SCRIPT_DIR, "results/similarity_matrix")

# Number of worker processes.
NUM_CORES_TO_USE = 24


def calculate_and_save_similarity(filepath, output_dir):
    """Compute and save a cosine-similarity matrix for one embedding file."""
    try:
        basename = os.path.basename(filepath)
        
        output_filename = f"similarity_{os.path.splitext(basename)[0]}.csv"
        output_path = os.path.join(output_dir, output_filename)
        
        df = pd.read_csv(filepath, sep=' ', header=None)
        df.dropna(inplace=True)
        
        if df.empty:
            return f"Empty file skipped: {basename}"
            
        otu_ids = df.iloc[:, 0].values
        embeddings = df.iloc[:, 1:].values
        
        similarity_matrix = cosine_similarity(embeddings)
        
        result_df = pd.DataFrame(similarity_matrix, index=otu_ids, columns=otu_ids)
        
        result_df.to_csv(output_path)
        
        return f"Completed: {basename} -> {output_filename}"
        
    except Exception as e:
        return f"Failed: {os.path.basename(filepath)} - Error: {e}"


def main():
    """Scan input files and dispatch tasks."""
    print(f"Input directory: {INPUT_FOLDER}")
    
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    print(f"Results will be saved to: {OUTPUT_FOLDER}")
    
    embedding_files = sorted(glob.glob(os.path.join(INPUT_FOLDER, '*.txt')))
    
    if not embedding_files:
        print("\nError: No .txt files were found in the input directory. Check the path.")
        return
        
    print(f"\nFound {len(embedding_files)} embedding files.")
    if NUM_CORES_TO_USE is None:
        print("Using all available CPU cores...")
    else:
        print(f"Using {NUM_CORES_TO_USE} CPU cores...")
    
    task_function = partial(calculate_and_save_similarity, output_dir=OUTPUT_FOLDER)
    with ProcessPoolExecutor(max_workers=NUM_CORES_TO_USE) as executor:
        results = executor.map(task_function, embedding_files)
        
        for result in results:
            print(result)
            
    print("\nAll files processed.")


if __name__ == "__main__":
    main()