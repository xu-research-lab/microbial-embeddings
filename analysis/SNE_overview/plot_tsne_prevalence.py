import json
import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from matplotlib.colors import LogNorm
import seaborn as sns
from collections import Counter

def load_pretrained_embeddings(vectors_file):
    with open(vectors_file, 'r') as f:
        vectors = {}
        for line in f:
            vals = line.rstrip().split(' ')
            if vals[0] != '<unk>':
                vectors[vals[0]] = [float(x) for x in vals[1:]]
    return vectors
embedding_dict = load_pretrained_embeddings('../../data/social_niche_embedding_100.txt')

def load_prevalence_dict(json_file):
    with open(json_file, 'r') as file:
        prevalence_dict = json.load(file)
    return prevalence_dict

prevalence_dict = load_prevalence_dict('prevalence_dict.json')
embedding_array = np.array(list(embedding_dict.values()))
prev_list = [prevalence_dict[fid] for fid in embedding_dict.keys()]

tsne = TSNE(n_components=2, random_state=42)
embedding_2d = tsne.fit_transform(embedding_array)
# 创建带颜色条的图形
plt.figure(figsize=(4, 3))
cmap = plt.get_cmap('viridis')  # 可以使用其他色谱如 'plasma', 'magma'
# 绘制散点图（颜色映射到prev_log值）
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
# 添加颜色条
cbar = plt.colorbar(scatter, pad=0.03)
cbar.set_label('Prevalence', fontsize=12)
cbar.ax.tick_params(labelsize=10)
# 添加标题和标签
# plt.title('t-SNE Visualization Colored by Prevalence', fontsize=14, pad=20)
plt.xlabel('t-SNE 1', fontsize=12)
plt.ylabel('t-SNE 2', fontsize=12)
plt.grid(alpha=0.2, linestyle='--')
# 优化显示
plt.tight_layout()
plt.show()






















