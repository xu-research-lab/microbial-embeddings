import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.metrics.pairwise import cosine_similarity
from skbio.stats.distance import mantel
import matplotlib.pyplot as plt
import os

# --- 1. Load the reference interaction matrix ---

df_randomE = pd.read_csv('data/E_global.csv', index_col=0)
species_randomE_dict = {}
for microbe_name, row_data in df_randomE.iterrows():
    resource_vector = row_data.tolist()
    species_randomE_dict[microbe_name] = resource_vector

# --- 2. Define helper functions ---
def load_pretrained_embeddings(vectors_file):
    """Load a GloVe embedding file."""
    if not os.path.exists(vectors_file):
        return None
    with open(vectors_file, 'r', encoding='utf-8') as f:
        vectors = {}
        for line in f:
            vals = line.rstrip().split(' ')
            # Require a token and at least one vector value; skip <unk>.
            if len(vals) > 1 and vals[0] != '<unk>':
                try:
                    vectors[vals[0]] = [float(x) for x in vals[1:]]
                except ValueError:
                    continue
    return vectors

def get_upper_triangle_values(df):
    """Return the matrix upper triangle, excluding the diagonal."""
    matrix = df.values
    # k=1 excludes the diagonal.
    return matrix[np.triu_indices_from(matrix, k=1)]

# --- 3. Define dataset sizes, replicates, and input paths ---
dataset_sizes = [1000, 2000, 5000, 10000, 20000, 40000, 80000, 160000, 200000]
repetitions = [1, 2, 3, 4, 5]

# embedding_file_template = os.path.join(EMBEDDING_DIR, 'processed_size_{size}', "subset_{rep}/result/embeddings_100.txt") 
embedding_file_template = os.path.join(
    "data/datasize_embeddings/embeddings_100_size{size}_{rep}.txt"
)


# --- 4. Calculate correlations for each dataset size and replicate ---
all_correlations_by_size = {size: [] for size in dataset_sizes} # Five replicates per dataset size.

print("\nCalculating correlations...")
for size_val in tqdm(dataset_sizes, desc="Dataset sizes"):
    correlations_for_current_size = []
    for rep_val in repetitions:
        current_embedding_file = embedding_file_template.format(size=size_val, rep=rep_val)
        
        embedding_dict_current = load_pretrained_embeddings(current_embedding_file)
        
        if embedding_dict_current is None or not embedding_dict_current:
            print(f"Warning: Could not load {current_embedding_file}, or the file is empty. Skipping.")
            correlations_for_current_size.append(np.nan) # Mark missing data as NaN.
            continue
            
        # Align the reference and embedding data by species name.
        common_keys = sorted(list(set(species_randomE_dict.keys()).intersection(embedding_dict_current.keys())))
        
        if len(common_keys) < 2: # Pearson correlation requires at least two data points.
            print(f"Warning: {current_embedding_file} shares fewer than two species with the reference ({len(common_keys)}); correlation cannot be calculated. Skipping.")
            correlations_for_current_size.append(np.nan)
            continue
            
        # Filter both vector sets in the same order.
        filtered_randomE_vectors_list = []
        for key in common_keys:
            filtered_randomE_vectors_list.append(species_randomE_dict[key])
        
        filtered_embedding_vectors_list = []
        for key in common_keys:
            filtered_embedding_vectors_list.append(embedding_dict_current[key])

        # Convert to NumPy arrays for cosine similarity.
        filtered_randomE_np = np.array(filtered_randomE_vectors_list)
        filtered_embedding_np = np.array(filtered_embedding_vectors_list)

        # Calculate cosine similarity matrices.
        df_sim_randomE_current = pd.DataFrame(cosine_similarity(filtered_randomE_np), index=common_keys, columns=common_keys)
        df_sim_embedding_current = pd.DataFrame(cosine_similarity(filtered_embedding_np), index=common_keys, columns=common_keys)
        
        # Calculate the Mantel correlation.
        try:
            array_sim_randomE = df_sim_randomE_current.to_numpy()
            np.fill_diagonal(array_sim_randomE, 0)
            array_sim_embedding = df_sim_embedding_current.to_numpy()
            np.fill_diagonal(array_sim_embedding, 0)
            corr, p_value, n = mantel(array_sim_randomE, array_sim_embedding, method='pearson', permutations=999)
            correlations_for_current_size.append(corr)
        except ValueError as e:
            print(f"Warning: Pearson correlation failed for {current_embedding_file}: {e}. The input array may have zero variance. Skipping.")
            correlations_for_current_size.append(np.nan)
            
    all_correlations_by_size[size_val] = correlations_for_current_size

# --- 5. Plot the dataset-size curve with error bars ---
print("\nPlotting the dataset-size curve...")

# Prepare plotting data.
plot_sizes = []
mean_correlations = []
std_dev_correlations = []

for size_val in dataset_sizes:
    corrs_for_size = all_correlations_by_size[size_val]
    # Exclude missing values from summary statistics.
    valid_corrs = [c for c in corrs_for_size if not np.isnan(c)]
    
    if valid_corrs: # Plot only sizes with valid correlations.
        plot_sizes.append(size_val)
        mean_correlations.append(np.mean(valid_corrs))
        std_dev_correlations.append(np.std(valid_corrs))
    else:
        print(f"Warning: Dataset size {size_val} has no valid correlations and will not be plotted.")

df_plot = pd.DataFrame({'size': plot_sizes, 'mean_correlation': mean_correlations, 'std_dev': std_dev_correlations})
df_plot.to_csv('results/df_plot_datasetsize.csv', index=False)
if not plot_sizes:
    print("Error: Not enough data to plot. Check the embedding files and paths.")
else:
    fig, ax = plt.subplots(figsize=(4, 3))

    # Plot means with standard-deviation error bars.
    ax.errorbar(plot_sizes, mean_correlations, yerr=std_dev_correlations, ecolor='black', elinewidth=2,
                fmt='-', color='black', alpha=0.5, capsize=3)


    ax.set_xlabel("Dataset size", fontsize=14)
    ax.set_ylabel("R", fontsize=14)
    # ax.legend(fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.6)

    ax.set_xscale('log')
    # ['1,000', '2,000', '5,000', '10,000', '20,000', '40,000', '80,000', '160,000', '200,000']
    custom_xtick_labels = ['1,000', '2,000', '5,000', '10,000', '20,000', '40,000', '80,000', '160,000', '200,000']
    ax.set_xticks(plot_sizes)
    ax.set_xticklabels(custom_xtick_labels, rotation=45, ha="right", fontsize=10)
    plt.xticks(fontsize=9, rotation=45, ha="right")

    plt.yticks(fontsize=12)
    
    # Add 10% padding to the observed y-axis range.
    if mean_correlations:
        min_val = min(m - s for m, s in zip(mean_correlations, std_dev_correlations) if not (np.isnan(m) or np.isnan(s)))
        max_val = max(m + s for m, s in zip(mean_correlations, std_dev_correlations) if not (np.isnan(m) or np.isnan(s)))
        if not (np.isnan(min_val) or np.isnan(max_val)):
            padding = (max_val - min_val) * 0.1
            ax.set_ylim(min_val - padding, max_val + padding)


    plt.tight_layout()
    plt.savefig('results/correlation_by_datasetsize.png', dpi=300)

print("\nScript completed.")
