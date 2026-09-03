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

def main(output_csv_filename="results/plot_dimension.csv"):
    # Define dimensions and number of runs
    dimensions = [25, 50, 100, 150, 200]
    num_runs = 5
    num_iterations = 100

    log_file_directory = "results/embedding_dim_log/embedding_dim_log"

    all_costs_per_dimension = defaultdict(list)

    # 1. Extract costs from each log file
    print("Starting log file parsing...")
    for dim in dimensions:
        for run_num in range(1, num_runs + 1):
            file_name = f"glove_dim{dim}_{run_num}.log"
            file_path = os.path.join(log_file_directory, file_name)
            
            print(f"Processing {file_path}...")
            costs_for_run = parse_glove_log(file_path)
            
            if costs_for_run: # Add only successfully parsed runs.
                all_costs_per_dimension[dim].append(costs_for_run)
            # Empty runs are handled by the completeness check below.

    # Filter out dimensions that didn't have all successful runs with consistent iterations
    valid_dimensions_data = {}
    for dim, runs_data in all_costs_per_dimension.items():
        if len(runs_data) == num_runs and all(len(run) == num_iterations for run in runs_data):
            valid_dimensions_data[dim] = np.array(runs_data)
        else:
            print(f"Warning: Dimension {dim} does not have {num_runs} complete runs of {num_iterations} iterations and will be skipped.")
            if len(runs_data) > 0:
                for i, run_d in enumerate(runs_data):
                    print(f"  Dimension {dim}, run {i+1}: {len(run_d)} iterations.")

    if not valid_dimensions_data:
        print("No valid data were parsed. Check the log files and paths.")
        return None

    # 2. Calculate mean and standard deviation for each dimension
    mean_costs = {}
    std_dev_costs = {}
    sorted_dims = sorted(valid_dimensions_data.keys())

    for dim in sorted_dims:
        runs_array = valid_dimensions_data[dim]
        mean_costs[dim] = np.mean(runs_array, axis=0)
        std_dev_costs[dim] = np.std(runs_array, axis=0)

    # Prepare the plotting data for CSV export.
    print("\nPreparing data for CSV export...")
    plotting_data_for_csv = []
    iterations_range = np.arange(1, num_iterations + 1)

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
        print(f"Plotting data saved to: {output_csv_filename}")
    except Exception as e:
        print(f"Error saving data to CSV: {e}")

    # 3. Plot the results
    plt.figure(figsize=(4, 3))
    
    
    error_bar_indices = np.array([0,19,39,59,79,99])

    colors = ['#C0C0BFFF', '#FFCD44FF', '#EE7C7AFF', '#4589C8FF', '#008F91FF']
    if len(colors) < len(sorted_dims):
        print(f"Warning: There are fewer colors ({len(colors)}) than dimensions ({len(sorted_dims)}); some dimensions may share colors.")
        # Reuse colors if more dimensions are added.
        num_colors_needed = len(sorted_dims)
        colors = [colors[j % len(colors)] for j in range(num_colors_needed)]


    for i,dim in enumerate(sorted_dims):
        means = mean_costs[dim]
        stds = std_dev_costs[dim]
        
        line, = plt.plot(iterations_range, means, label=f'Dim {dim}', color=colors[i])
        
        # Keep error-bar indices within the available iterations.
        valid_error_bar_indices = error_bar_indices[error_bar_indices < len(means)]

        if len(valid_error_bar_indices) > 0:
            plt.errorbar(iterations_range[valid_error_bar_indices], 
                         means[valid_error_bar_indices], 
                         yerr=stds[valid_error_bar_indices], 
                         fmt='none',
                         ecolor='black',
                         alpha=0.5,
                         capsize=3, 
                         # label=f'_Dim {dim} error'
                        )

    plt.xlabel("Iteration", fontsize=14)
    plt.ylabel("Loss", fontsize=14)
    plt.xticks(fontsize=12) 
    plt.yticks(fontsize=12) 
    plt.legend(loc='upper right') 
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xlim(1, num_iterations) 
    
    all_mean_costs_flat = np.concatenate([mean_costs[dim] for dim in sorted_dims if dim in mean_costs])
    if all_mean_costs_flat.size > 0:
        min_positive_cost = np.min(all_mean_costs_flat[all_mean_costs_flat > 0]) if np.any(all_mean_costs_flat > 0) else 0.001 # Avoid zero.
        if np.max(all_mean_costs_flat) < 0.1: 
            plt.ylim(bottom=min_positive_cost * 0.9) 

    plt.tight_layout()
    # Optional plot export.
    # output_plot_filename = "glove_cost_plot.png"
    # plt.savefig(output_plot_filename, dpi=300)
    plt.show()

    return df_plot_summary


if __name__ == "__main__":
    import pandas as pd
    df_plot_summary = main()
