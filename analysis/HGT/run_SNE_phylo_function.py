import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics.pairwise import cosine_similarity
from scipy.spatial.distance import pdist, squareform

from scipy.stats import norm
from sklearn.metrics import pairwise_distances
import json
import random
from phylodm import PhyloDM
import dendropy
from scipy.stats import ttest_ind

from scipy import stats
from tqdm import tqdm
from scipy.stats import gaussian_kde
from skbio.stats.distance import mantel, DistanceMatrix

plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

df_embedding = pd.read_csv("../../data/social_niche_embedding_100.txt",
                          header=None, sep=" ", low_memory=False, index_col=0)
df_embedding = df_embedding.drop("<unk>")

type = sys.argv[1]
### KO, EC, CAZY
func_file = f"Data/picrust/bac_{type}_predicted.tsv"
func_table = pd.read_csv(func_file, index_col=0, sep="\t", low_memory=False)

id_ = np.intersect1d(df_embedding.index.values, func_table.index.values)
id_ = np.random.choice(id_, size=1000, replace=False)

## function
func_table = func_table.loc[id_]
# func_table[func_table > 0] = 1
func_table = func_table.loc[:, (func_table != 0).any(axis=0)]
non_zero_counts = (func_table != 0).sum()
func_table = func_table.loc[:, non_zero_counts >= 100]


func_dist = pairwise_distances(func_table.values, metric='braycurtis', n_jobs=1) # braycurtis, jaccard
func_dist = pd.DataFrame(data=func_dist, index=func_table.index, columns=func_table.index)

## SNE
df_embedding = df_embedding.loc[func_dist.index.values]
cos_embed_pick = 1 - cosine_similarity(df_embedding)
cos_embed_pick = np.triu(cos_embed_pick, k=1)
cos_embed_pick = pd.DataFrame(data=cos_embed_pick, index=df_embedding.index, columns=df_embedding.index)
cos_embed_pick['id_1'] = cos_embed_pick.index.tolist()
cos_embed_pick = pd.melt(cos_embed_pick, id_vars="id_1")
cos_embed_pick = cos_embed_pick.loc[cos_embed_pick.value != 0]
cos_embed_pick.columns = ["id_1", "id_2", "cosine_dis"]

## phy
tree = dendropy.Tree.get_from_path('../SNE_overview/Data/tree.tre', schema='newick')
pdm = PhyloDM.load_from_dendropy(tree)
dm = pdm.dm(norm=False)
labels = pdm.taxa()
dm = pd.DataFrame(data=dm, index=labels, columns=labels) 
dm = dm.loc[df_embedding.index.values, df_embedding.index.values]

cosine_dist = 1 - cosine_similarity(df_embedding)
np.fill_diagonal(cosine_dist, 0)

dm_mat = DistanceMatrix(dm.values, ids=dm.index)
func_mat = DistanceMatrix(func_dist.values, ids=func_dist.index)
cosine_mat = DistanceMatrix(cosine_dist, ids=df_embedding.index.values)

mantel_phylo_func = mantel(dm_mat, func_mat, method='pearson', permutations=100)
mantel_cosine_func = mantel(cosine_mat, func_mat, method='pearson', permutations=100)

sne_dist = cos_embed_pick["cosine_dis"].values
func_dist_list = []
phy_dist = []
for i in range(cos_embed_pick.shape[0]):
    func_dist_list.append(func_dist.loc[cos_embed_pick.id_1.values[i], cos_embed_pick.id_2.values[i]])
    phy_dist.append(dm.loc[cos_embed_pick.id_1.values[i], cos_embed_pick.id_2.values[i]])
plot_df = pd.DataFrame({"phylo_dist" : phy_dist, "func_dist":func_dist_list, "sne_dist":sne_dist})

x = plot_df['phylo_dist']
y = plot_df['func_dist']

g = sns.JointGrid(data=plot_df, x='phylo_dist', y='func_dist', space=0, height=2)
hb = g.ax_joint.hexbin(
    x, y,
    gridsize=100,
    cmap='viridis',
    bins='log',
    mincnt=0,
    edgecolors='none',
    alpha=0.85
)
g.ax_joint.set_xlim(0, 4)
g.ax_joint.set_ylim(0, 1)

cax = g.fig.add_axes([0.92, 0.1, 0.02, 0.8])
plt.colorbar(hb, cax=cax, label='Log10(Density)')
sns.kdeplot(x, ax=g.ax_marg_x, fill=True, color='blue', alpha=0.2, bw_method=0.5, gridsize=50)
sns.kdeplot(y, ax=g.ax_marg_y, fill=True, color='blue', alpha=0.2, bw_method=0.5, gridsize=50, vertical=True)

coef = np.polyfit(x, y, 1)
poly1d_fn = np.poly1d(coef)
sorted_x = np.sort(x)
g.ax_joint.plot(sorted_x, poly1d_fn(sorted_x), linestyle='--', color='red')
x_min, x_max = x.min(), x.max()
y_min, y_max = y.min(), y.max()
text_x = x_min + 0.05 * (x_max - x_min)
text_y = 0.85
g.ax_joint.text(
    text_x, text_y,
    f'r = {mantel_phylo_func[0]:.4f}\np = {mantel_phylo_func[1]:.4f}',
    color='red', fontsize=7
)

for ax in [g.ax_marg_x, g.ax_marg_y]:
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.tick_params(axis='both', which='major', labelsize=8)

g.ax_marg_x.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
g.ax_marg_y.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

g.set_axis_labels('Phylogenetic Distance', 'Functional Distance', fontsize=10)
g.fig.suptitle(f'{type}', fontsize=16, y=1.02)

g.ax_joint.set_facecolor('#f0f0f0')
g.ax_joint.grid(True, linestyle='--', alpha=0.3)

plt.subplots_adjust(right=0.9)
plt.savefig(f'Figures/Phylo_func_{type}.pdf', format='pdf', dpi=300, bbox_inches='tight', pad_inches=0.1)


x = plot_df['sne_dist']
y = plot_df['func_dist']

g = sns.JointGrid(data=plot_df, x='sne_dist', y='func_dist', space=0, height=2)
hb = g.ax_joint.hexbin(
    x, y,
    gridsize=100,
    cmap='viridis',
    bins='log',
    mincnt=0,
    edgecolors='none',
    alpha=0.85
)
g.ax_joint.set_ylim(0, 1)
g.ax_joint.set_xlim(0, 2)
cax = g.fig.add_axes([0.92, 0.1, 0.02, 0.8])
plt.colorbar(hb, cax=cax, label='Log10(Density)')
sns.kdeplot(x, ax=g.ax_marg_x, fill=True, color='blue', alpha=0.2, bw_method=0.5, gridsize=50)
sns.kdeplot(y, ax=g.ax_marg_y, fill=True, color='blue', alpha=0.2, bw_method=0.5, gridsize=50, vertical=True)

coef = np.polyfit(x, y, 1)
poly1d_fn = np.poly1d(coef)
sorted_x = np.sort(x)
g.ax_joint.plot(sorted_x, poly1d_fn(sorted_x), linestyle='--', color='red')
x_min, x_max = x.min(), x.max()
y_min, y_max = y.min(), y.max()
text_x = x_min + 0.05 * (x_max - x_min)
text_y = 0.85
g.ax_joint.text(
    text_x, text_y,
    f'r = {mantel_cosine_func[0]:.4f}\np = {mantel_cosine_func[1]:.4f}',
    color='red', fontsize=7
)

for ax in [g.ax_marg_x, g.ax_marg_y]:
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.tick_params(axis='both', which='major', labelsize=8)

g.ax_marg_x.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
g.ax_marg_y.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

g.set_axis_labels('1 - SNE cosine', 'Functional Distance', fontsize=10)
g.fig.suptitle(f'{type}', fontsize=16, y=1.02)

g.ax_joint.set_facecolor('#f0f0f0')
g.ax_joint.grid(True, linestyle='--', alpha=0.3)

plt.subplots_adjust(right=0.9)
plt.savefig(f'Figures/SNE_func_{type}.pdf', format='pdf', dpi=300, bbox_inches='tight', pad_inches=0.1)
