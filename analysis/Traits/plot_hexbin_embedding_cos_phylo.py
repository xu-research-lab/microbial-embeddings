import json
import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib.pyplot as plt
from phylodm import PhyloDM
from scipy.stats import pearsonr

def load_pretrained_embeddings(vectors_file):
    with open(vectors_file, 'r') as f:
        vectors = {vals[0]: [float(x) for x in vals[1:]] 
                 for line in f if (vals := line.rstrip().split(' '))}
    return vectors

embedding_dict = load_pretrained_embeddings('/softerware/glove_embedding_new/new_data/p80_0/result/embeddings_100.txt')

# 获取有效微生物列表
microbes = [k for k in embedding_dict.keys() if k != '<unk>']
embeddings = np.array([embedding_dict[k] for k in microbes])

# 计算全量相似性矩阵
similarity_matrix = cosine_similarity(embeddings)
df_cos_sim = pd.DataFrame(similarity_matrix, index=microbes, columns=microbes)

# 加载进化距离矩阵
pdm = PhyloDM.load_from_newick_path('/softerware/tree_14k.nwk')
mat = pdm.dm(norm=False)
tree_labels = pdm.taxa()
df_phylo_dis = pd.DataFrame(mat, index=tree_labels, columns=tree_labels)



# 获取共同存在的微生物
common_microbes = sorted(set(microbes).intersection(df_phylo_dis.index.tolist()))
n = len(common_microbes)
print(f"Valid microbe pairs: {n*(n-1)//2}")

# 生成唯一组合对
def get_upper_triangle_values(df):
    """高效获取矩阵上三角值"""
    matrix = df.values
    return matrix[np.triu_indices_from(matrix, k=1)]

# 对齐并提取数据
x = get_upper_triangle_values(df_cos_sim.loc[common_microbes, common_microbes])
y = get_upper_triangle_values(df_phylo_dis.loc[common_microbes, common_microbes])

non_nan_idx = ~np.isnan(x) & ~np.isnan(y)
x = x[non_nan_idx]
y = y[non_nan_idx]

corr_coef, p_value = pearsonr(x, y)


# # 可视化优化
# plt.figure(figsize=(4, 3))
# hb = plt.hexbin(x, y,
#                 gridsize=100,  # 增加网格密度
#                 cmap='viridis', 
#                 bins='log',
#                 mincnt=5,      # 提高噪声过滤阈值
#                 edgecolors='none',
#                 alpha=0.85)
# plt.colorbar(hb, label='Data Density')
# plt.xlabel('SNE similarity', fontsize=14)
# plt.ylabel('Phylogenetic distance', fontsize=14)
# # plt.title('Pairwise Microbial Relationships Analysis', fontsize=16, pad=20)
# # plt.gca().set_facecolor('#f0f0f0')  # 添加背景色
# plt.grid(True, linestyle='--', alpha=0.3)
# plt.tight_layout()
# plt.show()


import seaborn as sns
from scipy.stats import gaussian_kde

# 创建画布和子图布局
fig = plt.figure(figsize=(4, 4))
gs = fig.add_gridspec(2, 2, 
                     width_ratios=[4, 1],
                     height_ratios=[1, 4],
                     hspace=0, wspace=0)
ax_main = fig.add_subplot(gs[1, 0])
ax_histx = fig.add_subplot(gs[0, 0], sharex=ax_main)
ax_histy = fig.add_subplot(gs[1, 1], sharey=ax_main)

# 主图：密度等高线
sns.kdeplot(
    x=x, y=y,
    ax=ax_main,
    cmap='Spectral_r',
    fill=True,
    thresh=0.05,
    levels=15,
    alpha=0.8
)

# 添加趋势线（OLS回归）
sns.regplot(
    x=x, y=y,
    ax=ax_main,
    scatter=False,
    color='black',
    line_kws={'lw':1.5, 'ls':'--'}
)

stats_text = (
    f"R={corr_coef:.3f}\n"
    f"P={p_value:.2e}" if p_value < 0.0001 else f"P-value = {p_value:.4f}"
)

# 边缘直方图（带核密度估计）
sns.kdeplot(
    x=x, 
    ax=ax_histx,
)

sns.kdeplot(
    y=y,
    ax=ax_histy,
)

# 美化图形
ax_main.set(
    xlabel='SNE similarity',
    ylabel='Phylogenetic distance',
    xlim=(-0.5, 1.1),
    ylim=(0, y.max()*1.05)
)
# ax_histx.set(title='')
# ax_histy.set(title='')

plt.tight_layout()
plt.show()

