import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.metrics.pairwise import cosine_similarity
from skbio.stats.distance import mantel
import matplotlib.pyplot as plt
import os # 用于文件路径操作

# --- 1. 加载基准数据 (E_global_v23.csv) ---

df_randomE = pd.read_csv('original_data_files/E_global_v46.csv', index_col=0)
species_randomE_dict = {}
for microbe_name, row_data in df_randomE.iterrows():
    resource_vector = row_data.tolist()
    species_randomE_dict[microbe_name] = resource_vector

# --- 2. 定义辅助函数 ---
def load_pretrained_embeddings(vectors_file):
    """加载 GloVe embedding 文件"""
    if not os.path.exists(vectors_file):
        # print(f"警告: Embedding 文件未找到: {vectors_file}")
        return None
    with open(vectors_file, 'r', encoding='utf-8') as f:
        vectors = {}
        for line in f:
            vals = line.rstrip().split(' ')
            # 确保行数据有效 (至少包含一个词和一个向量值)
            if len(vals) > 1 and vals[0] != '<unk>': # 排除未知词 <unk>
                try:
                    vectors[vals[0]] = [float(x) for x in vals[1:]]
                except ValueError:
                    # print(f"警告: 无法解析文件 {vectors_file} 中的行: {line.strip()}")
                    continue # 跳过无法解析的行
    return vectors

def get_upper_triangle_values(df):
    """高效获取DataFrame对应矩阵的上三角值 (不包括对角线)"""
    matrix = df.values
    # k=1 表示不包括对角线
    return matrix[np.triu_indices_from(matrix, k=1)]

# --- 3. 定义数据集大小、重复次数和文件路径模板 ---
dataset_sizes = [1000, 2000, 5000, 10000, 20000, 40000, 80000, 160000, 200000]
repetitions = [1, 2, 3, 4, 5]

# !! 修改这里 !! 指定您的 embedding 文件所在的目录和命名模板
# 假设文件在当前脚本运行目录下的 'downloaded_embeddings' 子目录中
# 如果文件就在当前目录，可以将 EMBEDDING_DIR 改为 "."
# EMBEDDING_DIR = "/softerware/analysis_datasize_simul_V4" # 或者例如 "downloaded_embeddings"
# embedding_file_template = os.path.join(EMBEDDING_DIR, 'processed_size_{size}', "subset_{rep}/result/embeddings_100.txt") 
embedding_file_template = os.path.join("/softerware/datasize_embeddings/embeddings_100_size{size}_{rep}.txt")


# --- 4. 计算每个数据集大小和重复实验的相关性 ---
all_correlations_by_size = {size: [] for size in dataset_sizes} # 用于存储每个size下的5个相关性值

print("\n开始计算相关性...")
for size_val in tqdm(dataset_sizes, desc="处理数据集大小"):
    correlations_for_current_size = []
    for rep_val in repetitions:
        current_embedding_file = embedding_file_template.format(size=size_val, rep=rep_val)
        
        embedding_dict_current = load_pretrained_embeddings(current_embedding_file)
        
        if embedding_dict_current is None or not embedding_dict_current:
            print(f"警告: 未能加载或为空 {current_embedding_file}, 跳过此文件。")
            correlations_for_current_size.append(np.nan) # 添加 NaN 表示缺失数据
            continue
            
        # 找出基准字典和当前 embedding 字典共有的键 (物种名)
        common_keys = sorted(list(set(species_randomE_dict.keys()).intersection(embedding_dict_current.keys())))
        
        if len(common_keys) < 2: # Pearson 相关性至少需要两个数据点（即至少两个物种的两个成对比较）
            print(f"警告: {current_embedding_file} 与基准数据的共同物种少于2个 ({len(common_keys)}个)，无法计算相关性。跳过。")
            correlations_for_current_size.append(np.nan)
            continue
            
        # 根据共同键过滤和排序向量，以确保对齐
        filtered_randomE_vectors_list = []
        for key in common_keys:
            filtered_randomE_vectors_list.append(species_randomE_dict[key])
        
        filtered_embedding_vectors_list = []
        for key in common_keys:
            filtered_embedding_vectors_list.append(embedding_dict_current[key])

        # 转换为 NumPy 数组以计算余弦相似度
        filtered_randomE_np = np.array(filtered_randomE_vectors_list)
        filtered_embedding_np = np.array(filtered_embedding_vectors_list)

        # 计算余弦相似度矩阵
        df_sim_randomE_current = pd.DataFrame(cosine_similarity(filtered_randomE_np), index=common_keys, columns=common_keys)
        df_sim_embedding_current = pd.DataFrame(cosine_similarity(filtered_embedding_np), index=common_keys, columns=common_keys)
        
        # 计算相关性
        try:
            array_sim_randomE = df_sim_randomE_current.to_numpy()
            np.fill_diagonal(array_sim_randomE, 0)
            array_sim_embedding = df_sim_embedding_current.to_numpy()
            np.fill_diagonal(array_sim_embedding, 0)
            corr, p_value, n = mantel(array_sim_randomE, array_sim_embedding, method='pearson', permutations=999)
            correlations_for_current_size.append(corr)
        except ValueError as e:
            print(f"警告: 计算 {current_embedding_file} 的 Pearson 相关性时出错: {e}。可能是由于输入数组方差为0。跳过。")
            correlations_for_current_size.append(np.nan)
            
    all_correlations_by_size[size_val] = correlations_for_current_size

# --- 5. 绘制折线图 + Error Bar ---
print("\n开始绘制折线图...")

# 准备绘图数据
plot_sizes = []
mean_correlations = []
std_dev_correlations = []

for size_val in dataset_sizes:
    corrs_for_size = all_correlations_by_size[size_val]
    # 过滤掉 NaN 值进行统计计算
    valid_corrs = [c for c in corrs_for_size if not np.isnan(c)]
    
    if valid_corrs: # 只有当该 size 下有有效的相关性值时才进行计算和绘图
        plot_sizes.append(size_val)
        mean_correlations.append(np.mean(valid_corrs))
        std_dev_correlations.append(np.std(valid_corrs))
    else:
        print(f"警告: 数据集大小 {size_val} 没有有效的相关性数据，将不在图表中显示。")

df_plot = pd.DataFrame({'size': plot_sizes, 'mean_correlation': mean_correlations, 'std_dev': std_dev_correlations})
df_plot.to_csv('df_plot_datasetsize.csv', index=False)
if not plot_sizes:
    print("错误：没有足够的数据来绘制折线图。请检查您的 embedding 文件和路径。")
else:
    fig, ax = plt.subplots(figsize=(4, 3)) # 调整图像大小

    # 绘制带误差线的折线图
    # yerr 参数接收标准差作为误差线的长度
    # capsize 控制误差线末端横线的长度
    # fmt 用于设置线条和标记的样式，例如 '-o' 表示实线和圆点标记
    ax.errorbar(plot_sizes, mean_correlations, yerr=std_dev_correlations, ecolor='black', elinewidth=2,
                fmt='-', color='black', alpha=0.5, capsize=3)


    ax.set_xlabel("Dataset size", fontsize=14)
    ax.set_ylabel("R", fontsize=14)
    # ax.legend(fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.6)

    ax.set_xscale('log') # 取消注释以使用对数刻度
    # ['1,000', '2,000', '5,000', '10,000', '20,000', '40,000', '80,000', '160,000', '200,000']
    custom_xtick_labels = ['1,000', '2,000', '5,000', '10,000', '20,000', '40,000', '80,000', '160,000', '200,000']
    ax.set_xticks(plot_sizes) # 尝试在每个数据点位置设置刻度
    ax.set_xticklabels(custom_xtick_labels, rotation=45, ha="right", fontsize=10)
    plt.xticks(fontsize=9, rotation=45, ha="right")

    plt.yticks(fontsize=12)
    
    # 自动调整Y轴范围以更好地显示数据
    if mean_correlations: # 如果有有效数据点
        min_val = min(m - s for m, s in zip(mean_correlations, std_dev_correlations) if not (np.isnan(m) or np.isnan(s)))
        max_val = max(m + s for m, s in zip(mean_correlations, std_dev_correlations) if not (np.isnan(m) or np.isnan(s)))
        if not (np.isnan(min_val) or np.isnan(max_val)):
            padding = (max_val - min_val) * 0.1 # 上下留10%的边距
            ax.set_ylim(min_val - padding, max_val + padding)


    plt.tight_layout() # 自动调整子图参数，使之填充整个图像区域
    plt.savefig('correlation_by_datasetsize.png', dpi=300)

print("\n脚本执行完毕。")
