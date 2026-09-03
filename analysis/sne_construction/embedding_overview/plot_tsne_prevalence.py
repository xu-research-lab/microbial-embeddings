from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from sklearn.manifold import TSNE


matplotlib.use("Agg")
import matplotlib.pyplot as plt


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
EMBEDDING_FILE = REPO_ROOT / "data/social_niche_embedding_100.txt"
PREVALENCE_FILE = (
    SCRIPT_DIR.parent / "data/embedding_overview/OTU_prevalence_abundance.csv"
)
OUTPUT_FILE = SCRIPT_DIR / "results/figures/tsne_prevalence.png"


embedding = pd.read_csv(EMBEDDING_FILE, sep=r"\s+", header=None, index_col=0)
embedding.drop(index="<unk>", errors="ignore", inplace=True)
prevalence = pd.read_csv(PREVALENCE_FILE, index_col="OTU")["prevalence"]

feature_ids = embedding.index.intersection(prevalence.index)
coordinates = TSNE(n_components=2, random_state=42).fit_transform(
    embedding.loc[feature_ids].to_numpy()
)

fig, ax = plt.subplots(figsize=(4, 3))
points = ax.scatter(
    coordinates[:, 0],
    coordinates[:, 1],
    c=prevalence.loc[feature_ids],
    cmap="viridis",
    norm=LogNorm(),
    alpha=0.6,
    s=10,
    edgecolors="white",
    linewidth=0.3,
)
fig.colorbar(points, ax=ax, pad=0.03, label="Prevalence")
ax.set(xlabel="t-SNE 1", ylabel="t-SNE 2")
ax.grid(alpha=0.2, linestyle="--")
fig.tight_layout()

OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUTPUT_FILE, dpi=300)
print(OUTPUT_FILE)
