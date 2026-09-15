from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from matplotlib.colors import LogNorm
import seaborn as sns
from collections import Counter

# Run this script from the repository root.
REPO_ROOT = Path.cwd().resolve()
if not (REPO_ROOT / "analysis/sne_construction").is_dir():
    raise RuntimeError("Run this script from the repository root.")

DATA_DIR = REPO_ROOT / "data"
ANALYSIS_DIR = REPO_ROOT / "analysis/sne_construction"
OVERVIEW_DIR = ANALYSIS_DIR / "embedding_overview"
OVERVIEW_DATA_DIR = ANALYSIS_DIR / "data/embedding_overview"
FIGURES_DIR = OVERVIEW_DIR / "results/figures"

EMBEDDING_FILE = DATA_DIR / "social_niche_embedding_100.txt"
PREVALENCE_FILE = OVERVIEW_DATA_DIR / "OTU_prevalence_abundance.csv"
OUTPUT_FILE = FIGURES_DIR / "tsne_prevalence.png"


def load_pretrained_embeddings(vectors_file):
    with open(vectors_file, 'r') as f:
        vectors = {}
        for line in f:
            vals = line.rstrip().split(' ')
            if vals[0] != '<unk>':
                vectors[vals[0]] = [float(x) for x in vals[1:]]
    return vectors


embedding_dict = load_pretrained_embeddings(EMBEDDING_FILE)
prevalence_table = np.genfromtxt(
    PREVALENCE_FILE, delimiter=",", names=True, dtype=None, encoding="utf-8"
)
# Match prevalence values to the embedding row order.
prevalence_dict = dict(zip(prevalence_table["OTU"], prevalence_table["prevalence"]))
embedding_array = np.array(list(embedding_dict.values()))
prev_list = [prevalence_dict[fid] for fid in embedding_dict.keys()]

# Project the 100-dimensional SNE into two dimensions.
tsne = TSNE(n_components=2, random_state=42)
embedding_2d = tsne.fit_transform(embedding_array)
# Create a figure with a color bar
plt.figure(figsize=(4, 3))
cmap = plt.get_cmap('viridis')  # Other colormaps include 'plasma' and 'magma'
# Map point colors to prevalence
scatter = plt.scatter(
    embedding_2d[:, 0], 
    embedding_2d[:, 1], 
    c=prev_list,
    cmap=cmap,
    norm=LogNorm(),
    alpha=0.6,
    s=10,
    edgecolor='w',
    linewidth=0.3
)
# Add the color bar
cbar = plt.colorbar(scatter, pad=0.03)
cbar.set_label('Prevalence', fontsize=12)
cbar.ax.tick_params(labelsize=10)
# Add labels
# plt.title('t-SNE Visualization Colored by Prevalence', fontsize=14, pad=20)
plt.xlabel('t-SNE 1', fontsize=12)
plt.ylabel('t-SNE 2', fontsize=12)
plt.grid(alpha=0.2, linestyle='--')
# Adjust the layout
plt.tight_layout()
OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUTPUT_FILE, dpi=300, bbox_inches="tight")
plt.show()






















