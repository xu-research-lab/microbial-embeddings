# Load necessary libraries
# 加载所需的库
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tools) # For file path manipulation
library(scales) # For log scale formatting

# Specify the input directory containing the CSV files
# 指定包含CSV文件的输入目录
input_dir <- "v46/plot_csvs"

# Get a list of all CSV files in the directory
# 获取目录中所有CSV文件的列表
csv_files <- list.files(path = input_dir, pattern = "\\.csv$", full.names = TRUE)

# --- Pre-computation Step: Find the global y-axis range across all files ---
# --- 预计算步骤：查找所有文件的全局y轴范围 ---
global_min_ratio <- Inf
global_max_ratio <- -Inf

# Loop through each file just to determine the overall range
# 循环遍历每个文件以确定总体范围
for (file_path in csv_files) {
  # Load and preprocess data, filtering for valid values for log scale
  # 加载并预处理数据，筛选适用于对数刻度的有效值
  data <- read.csv(file_path) %>%
    filter(ratio != 0 & is.finite(ratio)) # Log scale cannot handle zero or negative values

  # Clean data by removing outliers
  # 通过移除异常值来清理数据
  if (nrow(data) > 0) {
    data_clean <- data
    # %>%
    #   group_by(category) %>%
    #   mutate(
    #     Q1 = quantile(ratio, 0.25, na.rm = TRUE),
    #     Q3 = quantile(ratio, 0.75, na.rm = TRUE),
    #     IQR = Q3 - Q1,
    #     lower = Q1 - 1.5 * IQR,
    #     upper = Q3 + 1.5 * IQR
    #   ) %>%
    #   filter(ratio >= lower & ratio <= upper) %>%
    #   ungroup()

    # Update global min and max if the cleaned data is not empty
    # 如果清理后的数据不为空，则更新全局最小值和最大值
    if (nrow(data_clean) > 0) {
      current_min <- min(data_clean$ratio, na.rm = TRUE)
      current_max <- max(data_clean$ratio, na.rm = TRUE)
      if (current_min < global_min_ratio) global_min_ratio <- current_min
      if (current_max > global_max_ratio) global_max_ratio <- current_max
    }
  }
}

# Define the unified y-axis limits and breaks for the log scale
# 为对数刻度定义统一的y轴限制和刻度
# Round down to the nearest power of 10 for the lower limit, and up for the upper limit
y_limits <- c(10^floor(log10(global_min_ratio)), 10^ceiling(log10(global_max_ratio)))
y_breaks <- 10^seq(floor(log10(y_limits[1])), ceiling(log10(y_limits[2])))
y_temp <- y_breaks
y_breaks <- c(y_breaks[1], y_breaks[3], y_breaks[5], y_breaks[7], y_breaks[9], y_breaks[11], y_breaks[13], y_breaks[15])

# --- Main Plotting Loop ---
# --- 主绘图循环 ---
# Loop through each CSV file to process and plot
# 循环处理每个CSV文件并绘图
for (file_path in csv_files) {
  # --- 1. Load and preprocess the data ---
  # --- 1. 读取并预处理数据 ---
  data <- read.csv(file_path) %>%
    mutate(
      category = recode(category,
        "Sim_ik_Sim_jk" = "Sim.(i,k)/Sim.(j,k)",
        "Sim_ik_Unsim_jk" = "Sim.(i,k)/Dissim.(j,k)",
        "Unsim_ik_Sim_jk" = "Dissim.(i,k)/Sim.(j,k)",
        "Unsim_ik_Unsim_jk" = "Dissim.(i,k)/Dissim.(j,k)"
      )
    ) %>%
    # Filter out values that are not suitable for a log scale
    filter(ratio != 0 & is.finite(ratio))

  # --- 2. Clean data by removing outliers ---
  # --- 2. 通过移除异常值来清理数据 ---
  data_clean <- data
  # %>%
  #   group_by(category) %>%
  #   mutate(
  #     Q1 = quantile(ratio, 0.25, na.rm = TRUE),
  #     Q3 = quantile(ratio, 0.75, na.rm = TRUE),
  #     IQR = Q3 - Q1,
  #     lower = Q1 - 1.5 * IQR,
  #     upper = Q3 + 1.5 * IQR
  #   ) %>%
  #   filter(ratio >= lower & ratio <= upper) %>%
  #   ungroup()

  # Set the order of categories for the x-axis
  # 设置x轴类别的顺序
  data_clean$category <- factor(data_clean$category, levels = c(
    "Sim.(i,k)/Sim.(j,k)",
    "Sim.(i,k)/Dissim.(j,k)",
    "Dissim.(i,k)/Sim.(j,k)",
    "Dissim.(i,k)/Dissim.(j,k)"
  ))

  # --- 3. Identify significantly different pairs ---
  # --- 3. 筛选显著差异的组别对 ---
  categories <- unique(data_clean$category)
  pairs <- combn(categories, 2, simplify = FALSE)

  signif_pairs <- list()
  for (pair in pairs) {
    excluded_pair <- c("Sim.(i,k)/Sim.(j,k)", "Dissim.(i,k)/Dissim.(j,k)")
    if (setequal(pair, excluded_pair)) {
      next
    }
    if (all(pair %in% levels(data_clean$category))) {
      group1 <- filter(data_clean, category == pair[1]) %>% dplyr::select(ratio)
      group2 <- filter(data_clean, category == pair[2]) %>% dplyr::select(ratio)
      if (nrow(group1) > 1 && nrow(group2) > 1) {
        test_result <- tryCatch(t.test(group1, group2), error = function(e) NULL)
        if (!is.null(test_result) && test_result$p.value < 0.05) {
          signif_pairs <- c(signif_pairs, list(as.character(pair)))
        }
      }
    }
  }
  signif_pairs <- list(
    c("Sim.(i,k)/Sim.(j,k)", "Sim.(i,k)/Dissim.(j,k)"),
    c("Sim.(i,k)/Dissim.(j,k)", "Dissim.(i,k)/Sim.(j,k)"),
    c("Sim.(i,k)/Dissim.(j,k)", "Dissim.(i,k)/Dissim.(j,k)"),
    c("Dissim.(i,k)/Sim.(j,k)", "Dissim.(i,k)/Dissim.(j,k)"),
    c("Dissim.(i,k)/Sim.(j,k)", "Sim.(i,k)/Sim.(j,k)")
  )

  # --- 4. Plot the graph ---
  # --- 4. 绘制图形 ---
  p1 <- ggplot(data_clean, aes(x = category, y = ratio)) +
    geom_violin(weight = 0.1, alpha = 0.5, fill = "#0099FF") +
    geom_boxplot(width = 0.2, fatten = 1.5, outlier.shape = NA) +
    {
      if (length(signif_pairs) > 0) {
        stat_compare_means(
          comparisons = signif_pairs,
          method = "t.test",
          label = "p.signif",
          tip.length = 0.02,
          step.increase = 0.05,
          vjust = 0.5,
          size = 2.5
        )
      }
    } +
    theme_bw() +
    geom_hline(yintercept = 1, colour = "red", linewidth = 0.3) + # Fixed color specification
    # Apply log scale with unified limits and smart decimal labels
    # 应用具有统一限制和智能十进制标签的对数刻度
    scale_y_log10(
      limits = c(y_limits[1], NA),
      breaks = y_breaks,
      labels = function(x) format(x, scientific = T, drop0trailing = TRUE),
      expand = expansion(mult = c(0, 0.07))
    ) +
    theme(
      axis.text.x = element_text(
        angle = 30,
        hjust = 1
      ),
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 15, color = "black"),
      panel.spacing.x = unit(0.01, "lines"),
      axis.title.x = element_text(hjust = 1),
      panel.border = element_rect(color = "black", fill = NA)
    ) +
    labs(y = "Co-occur. prob. ratio", x = "Semantic relationship categories") +
    scale_x_discrete(expand = expansion(mult = 0.1))

  # --- 5. Save the plot with a modified filename ---
  # --- 5. 使用修改后的文件名保存图形 ---
  base_name <- basename(file_path)
  new_name <- sub("df_plot_ratio", "boxplot", base_name)
  new_name <- sub("\\.csv$", "_log10.png", new_name)

  output_path <- file.path("v46", new_name)

  ggsave(p1, filename = output_path, width = 3.8, height = 3.8, dpi = 300)

  cat("Processed:", base_name, "-> Saved as:", new_name, "\n")
}
