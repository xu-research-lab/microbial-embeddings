import re
import os
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

def parse_glove_log(file_path):
    """
    Parses a GloVe training log file to extract cost per iteration.

    Args:
        file_path (str): The path to the log file.

    Returns:
        list: A list of 100 cost values (floats). Returns empty if parsing fails
              or fewer than 100 costs are found.
    """
    costs = []
    # Regex to find lines with iteration and cost
    # Example line: 05/17/25 - 09:56.07PM, iter: 100, cost: 0.048117
    cost_line_regex = re.compile(r"iter:\s*\d+,\s*cost:\s*(\d+\.\d+)")

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                match = cost_line_regex.search(line)
                if match:
                    costs.append(float(match.group(1)))
        # Ensure we have exactly 100 costs, otherwise it might indicate an issue
        if len(costs) == 100:
            return costs
        else:
            print(f"Warning: Found {len(costs)} cost entries in {file_path}, expected 100.")
            # If you want to be strict and only accept perfect files:
            # return []
            # If you want to use partial data (up to 100):
            return costs[:100] if len(costs) > 0 else []
    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        return []
    except Exception as e:
        print(f"Error parsing file {file_path}: {e}")
        return []

def main(output_csv_filename="plot_dimension.csv"): # 新增 output_csv_filename 参数
    # Define dimensions and number of runs
    dimensions = [25, 50, 100, 150, 200]
    num_runs = 5
    num_iterations = 100

    log_file_directory = "/softerware/glove_embedding_new/embedding_dim" # 请确保这是正确的路径

    all_costs_per_dimension = defaultdict(list)

    # 1. Extract costs from each log file
    print("Starting log file parsing...")
    for dim in dimensions:
        for run_num in range(1, num_runs + 1):
            file_name = f"glove_dim{dim}_{run_num}.log"
            file_path = os.path.join(log_file_directory, file_name)
            
            print(f"Processing {file_path}...")
            costs_for_run = parse_glove_log(file_path)
            
            if costs_for_run: # 只有当成功解析出成本时才添加
                all_costs_per_dimension[dim].append(costs_for_run)
            # 注意：如果costs_for_run为空列表，这里不会添加，后续的长度检查会处理

    # Filter out dimensions that didn't have all successful runs with consistent iterations
    valid_dimensions_data = {}
    for dim, runs_data in all_costs_per_dimension.items():
        if len(runs_data) == num_runs and all(len(run) == num_iterations for run in runs_data):
            valid_dimensions_data[dim] = np.array(runs_data)
        else:
            print(f"警告: 维度 {dim} 没有 {num_runs} 个包含 {num_iterations} 次迭代的完整运行。将跳过此维度。")
            if len(runs_data) > 0:
                for i, run_d in enumerate(runs_data):
                    print(f"  维度 {dim}, 运行 {i+1} 有 {len(run_d)} 次迭代。")

    if not valid_dimensions_data:
        print("解析后没有有效数据可供绘图或保存。请检查日志文件和路径。")
        return None # 或者返回一个空字典，取决于后续如何处理

    # 2. Calculate mean and standard deviation for each dimension
    mean_costs = {}
    std_dev_costs = {}
    sorted_dims = sorted(valid_dimensions_data.keys())

    for dim in sorted_dims:
        runs_array = valid_dimensions_data[dim]
        mean_costs[dim] = np.mean(runs_array, axis=0)
        std_dev_costs[dim] = np.std(runs_array, axis=0)

    # --- 新增：准备数据并写入CSV文件 ---
    print("\n准备数据用于CSV导出...")
    plotting_data_for_csv = []
    iterations_range = np.arange(1, num_iterations + 1) # 迭代次数从1到100

    for dim in sorted_dims:
        means = mean_costs[dim]
        stds = std_dev_costs[dim]
        for i in range(num_iterations):
            plotting_data_for_csv.append({
                'Dimension': dim,
                'Iteration': iterations_range[i],
                'MeanCost': means[i],
                'StdDevCost': stds[i]
            })
    
    df_plot_summary = pd.DataFrame(plotting_data_for_csv)
    
    try:
        df_plot_summary.to_csv(output_csv_filename, index=False)
        print(f"绘图数据已成功保存到: {output_csv_filename}")
    except Exception as e:
        print(f"保存数据到CSV时发生错误: {e}")
    # --- CSV写入结束 ---

    # 3. Plot the results
    plt.figure(figsize=(4, 3)) # 原始尺寸，如果需要可以调整
    # iterations 变量在上面已定义为 iterations_range
    
    error_every_n_points = 20 
    error_bar_indices = np.array([0,19,39,59,79,99]) # 确保这些索引在0到num_iterations-1范围内

    colors = ['#C0C0BFFF', '#FFCD44FF', '#EE7C7AFF', '#4589C8FF', '#008F91FF']
    if len(colors) < len(sorted_dims): # 如果颜色不够，进行扩展或提示
        print(f"警告: 颜色数量 ({len(colors)}) 少于维度数量 ({len(sorted_dims)})。部分维度可能共享颜色或出错。")
        # 可以简单地循环使用颜色
        num_colors_needed = len(sorted_dims)
        colors = [colors[j % len(colors)] for j in range(num_colors_needed)]


    for i,dim in enumerate(sorted_dims):
        means = mean_costs[dim]
        stds = std_dev_costs[dim]
        
        line, = plt.plot(iterations_range, means, label=f'Dim {dim}', color=colors[i])
        
        # 确保 error_bar_indices 不会超出 means/stds 数组的界限
        valid_error_bar_indices = error_bar_indices[error_bar_indices < len(means)]

        if len(valid_error_bar_indices) > 0:
            plt.errorbar(iterations_range[valid_error_bar_indices], 
                         means[valid_error_bar_indices], 
                         yerr=stds[valid_error_bar_indices], 
                         fmt='none',
                         ecolor='black', # 之前是line.get_color()，改为固定黑色以确保可见性      
                         alpha=0.5,
                         capsize=3, 
                         # label=f'_Dim {dim} error' # 通常不需要单独为error bar加图例
                        )

    plt.xlabel("Iteration", fontsize=14)
    plt.ylabel("Loss", fontsize=14)
    plt.xticks(fontsize=12) 
    plt.yticks(fontsize=12) 
    plt.legend(loc='upper right') 
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xlim(1, num_iterations) 
    
    all_mean_costs_flat = np.concatenate([mean_costs[dim] for dim in sorted_dims if dim in mean_costs])
    if all_mean_costs_flat.size > 0 : # 确保数组不为空
        min_positive_cost = np.min(all_mean_costs_flat[all_mean_costs_flat > 0]) if np.any(all_mean_costs_flat > 0) else 0.001 # 避免0
        if np.max(all_mean_costs_flat) < 0.1: 
            plt.ylim(bottom=min_positive_cost * 0.9) 

    plt.tight_layout()
    # 保存图像 (如果需要，可以取消注释并指定文件名)
    # output_plot_filename = "glove_cost_plot.png"
    # plt.savefig(output_plot_filename, dpi=300)
    # print(f"图像已保存到: {output_plot_filename}")
    plt.show()

    return df_plot_summary


if __name__ == "__main__":
    import pandas as pd
    df_plot_summary = main()