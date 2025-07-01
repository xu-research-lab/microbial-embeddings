import os
import glob
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
from concurrent.futures import ProcessPoolExecutor
from functools import partial

# --- 配置 ---

# 1. 输入文件夹
INPUT_FOLDER = '/home/cjj/projects/membed_local/membed/analysis_distance_true_datasize/embeding_filter'

# 2. 输出文件夹
OUTPUT_FOLDER = os.path.join('./', 'similarity_matrix')

# 3. 自定义使用的CPU核心数  # <--- 关键改动 1: 添加新配置项
#    - 设置为您想使用的具体数字，例如: 4 或 8
#    - 如果设置为 None，它将默认使用您计算机上所有可用的核心
#    - 推荐做法: os.cpu_count() - 1，可以为系统保留一个核心，防止电脑卡顿
NUM_CORES_TO_USE = 24


def calculate_and_save_similarity(filepath, output_dir):
    """
    为单个嵌入文件计算并保存余弦相似度矩阵。
    (此函数内容保持不变)
    """
    try:
        basename = os.path.basename(filepath)
        # if basename !='russell_rao_1_100.txt':
        #     return f"跳过: {basename}"
        
        output_filename = f"similarity_{os.path.splitext(basename)[0]}.csv"
        output_path = os.path.join(output_dir, output_filename)
        
        df = pd.read_csv(filepath, sep=' ', header=None)
        df.dropna(inplace=True)
        
        if df.empty:
            return f"文件为空，已跳过: {basename}"
            
        otu_ids = df.iloc[:, 0].values
        embeddings = df.iloc[:, 1:].values
        
        similarity_matrix = cosine_similarity(embeddings)
        
        result_df = pd.DataFrame(similarity_matrix, index=otu_ids, columns=otu_ids)
        
        result_df.to_csv(output_path)
        
        return f"处理完成: {basename} -> {output_filename}"
        
    except Exception as e:
        return f"处理失败: {os.path.basename(filepath)} - 错误: {e}"


def main():
    """
    主函数，负责扫描文件和分派任务
    """
    print(f"输入文件夹: {INPUT_FOLDER}")
    
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    print(f"结果将保存在: {OUTPUT_FOLDER}")
    
    embedding_files = glob.glob(os.path.join(INPUT_FOLDER, '*.txt'))
    
    if not embedding_files:
        print("\n错误：在指定文件夹下未找到任何 .txt 文件。请检查路径。")
        return
        
    print(f"\n找到 {len(embedding_files)} 个嵌入文件。")
    # 打印将要使用的核心数
    if NUM_CORES_TO_USE is None:
        print(f"将使用所有可用的CPU核心进行计算...")
    else:
        print(f"将使用 {NUM_CORES_TO_USE} 个CPU核心进行计算...")
    
    task_function = partial(calculate_and_save_similarity, output_dir=OUTPUT_FOLDER)
    
    # <--- 关键改动 2: 在此处传入 max_workers 参数
    with ProcessPoolExecutor(max_workers=NUM_CORES_TO_USE) as executor:
        results = executor.map(task_function, embedding_files)
        
        for result in results:
            print(result)
            
    print("\n所有文件处理完毕！")


if __name__ == "__main__":
    main()