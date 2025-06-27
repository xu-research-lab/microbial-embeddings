# load biom
library(biomformat)
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggraph)
library(tidygraph)
library(viridis)
micro_abundance <- read_biom("Data/AGP_test.biom")
otu <- biom_data(micro_abundance) %>% as.matrix()
# random select sample
set.seed(1234)
sample_index <- sample(1:ncol(otu), 1)
sample_data <- otu[, sample_index, drop = FALSE]
sample_data <- sample_data[sample_data[,1]> 0,,drop = F]
sample_data_tss <- sample_data / sum(sample_data)
# 分层抽样 20个otus 
quantiles <- quantile(sample_data_tss, probs = c(0, 0.25, 0.5, 0.75, 0.99))
strata <- cut(sample_data_tss, 
              breaks = quantiles, 
              include.lowest = TRUE, 
              labels = c("low", "medium", "high", "very_high"))  # 修正为四层标签

# 设定四层抽样数量（总和20）
samples_per_stratum <- c(low = 5, medium = 5, high = 5, very_high = 5)  # 每层5个

# 检查并抽样（按实际存在的分层标签）
selected_indices <- lapply(names(samples_per_stratum), function(s) {
  indices <- which(strata == s)
  if (length(indices) == 0) return(NULL)
  sample_size <- min(samples_per_stratum[s], length(indices))
  sample(indices, size = sample_size)
}) %>% unlist()


# 提取数据
sample_data <- sample_data[selected_indices, , drop = FALSE]
#sample_data_tss <- sample_data_tss[selected_indices, , drop = FALSE]
sample_data_tss <- sample_data / sum(sample_data)

remove(micro_abundance, otu)
# Russell rao (a/n)
cal_russell_rao_table <- function(otu_table) {
  # 转换为二进制数据（存在/不存在）
  binary_data <- as.matrix(otu_table > 0)
  
  # 计算样本间共同存在物种数
  common_counts <- binary_data %*% t(binary_data)
  
  # 总样本数
  total_species <- ncol(binary_data)
  
  # 计算 Russell-Rao 相似性矩阵
  cooccurrence_matrix <- common_counts / total_species
  
  # 提取上三角（不包含对角线）的索引
  upper_indices <- which(upper.tri(cooccurrence_matrix, diag = FALSE), arr.ind = TRUE)
  
  # 构建三列数据框
  result <- data.frame(
    micro_1 = rownames(cooccurrence_matrix)[upper_indices[, 1]],
    micro_2 = colnames(cooccurrence_matrix)[upper_indices[, 2]],
    Co_occurrence = cooccurrence_matrix[upper_indices]
  )
  
  return(result)
}
Russell_rao <- cal_russell_rao_table(sample_data)

  #plot 
  # 节点数据准备
  sample_data_temp <- sample_data_tss %>% 
    as.data.frame() %>% 
    setNames("Relative_Abundance") %>% 
    rownames_to_column("point")
  
  # 边数据校验
  Russell_rao_edges <- as.data.frame(Russell_rao)
  colnames(Russell_rao_edges) <- c("from", "to", "co_value")
  
  # 创建图形对象
  Russell_rao_graph <- tbl_graph(
    nodes = sample_data_temp,
    edges = Russell_rao_edges,
    directed = FALSE
  ) %>% 
    activate(nodes) %>% 
    arrange(desc(Relative_Abundance))
  
  Russell_rao_plot <- ggraph(Russell_rao_graph, layout = "linear", circular = TRUE) +
    geom_edge_fan(
      aes(color = co_value),
      alpha = 0.4,
      edge_width = 0.5
    ) +
    scale_edge_color_gradientn(
      colors = c("#d7191c")
    ) +
    geom_node_point(
      aes(size = Relative_Abundance),
      color = "grey30"
    ) +
    scale_size_continuous(range = c(2, 8)) +
    labs(title = "Russell rao") +
    theme_void() +
    coord_fixed() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold") ,
      legend.position = "none"
    )
  
  Russell_rao_plot

# Russeellrao weight ()
  cooccurrence_Russeellrao <- matrix(0, nrow = nrow(sample_data), ncol = nrow(sample_data))
  row_order <- order(sample_data, decreasing = TRUE)
  
  # 按行和排序（保持列顺序不变）
  sample_data <- sample_data[row_order, ,drop = FALSE]
  rownames(cooccurrence_Russeellrao) <- rownames(sample_data)
  colnames(cooccurrence_Russeellrao) <- rownames(sample_data)
  
  for (i in 1:(nrow(sample_data)-1)) {
    for (j in (i+1):nrow(sample_data)) {
      cooccurrence_Russeellrao[i, j] <- 1/(abs(i-j))
    }
  }
  # 提取上三角（不包含对角线）的索引
  upper_indices <- which(upper.tri(cooccurrence_Russeellrao, diag = FALSE), arr.ind = TRUE)
  # 构建三列数据框
  Russell_rao_weight <- data.frame(
    micro_1 = rownames(cooccurrence_Russeellrao)[upper_indices[, 1]],
    micro_2 = colnames(cooccurrence_Russeellrao)[upper_indices[, 2]],
    Co_occurrence = cooccurrence_Russeellrao[upper_indices]
  )
  #plot
  #边数据
  Russell_rao_weight_edges <- as.data.frame(Russell_rao_weight)
  colnames(Russell_rao_weight_edges) <- c("from", "to", "co_value")
  
  # 创建图形对象
  Russell_rao_weight_graph <- tbl_graph(
    nodes = sample_data_temp,
    edges = Russell_rao_weight_edges,
    directed = FALSE
  ) %>% 
    activate(nodes) %>% 
    arrange(desc(Relative_Abundance))
  
  Russell_rao_weight_plot <- ggraph(Russell_rao_weight_graph, layout = "linear", circular = TRUE) +
    geom_edge_fan(
      aes(color = co_value),
      alpha = 0.7,
      edge_width = 0.5
    ) +
    scale_edge_color_gradientn(
      colors = c("#000080", "#2c7bb6", "#7bc8f6", "#ffffbf", "#d7191c"),
      name = NULL,
      guide = guide_edge_colorbar(  # 将参数移到 guide 内部
        breaks = c(0, 1),          # 断点位置
        labels = c("low", "high"), # 对应标签
        show.limits = TRUE         # 强制显示端点标签
      )
    ) +
    geom_node_point(
      aes(size = Relative_Abundance),
      color = "grey30"
    ) +
    scale_size_continuous(range = c(2, 8)) +
    labs(title = "Russellrao weight") +
    theme_void() +
    coord_fixed() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
  Russell_rao_weight_plot
  


# Jaccard (a/(n - c))
  # 计算 Jaccard 距离矩阵
  dist_binary <- vegdist(sample_data, method = "jaccard", binary = TRUE)
  
  # 转换为 Jaccard 相似性矩阵
  sim_matrix <- 1 - as.matrix(dist_binary)
  
  # 提取上三角（不包含对角线）的索引
  upper_indices <- which(upper.tri(sim_matrix, diag = FALSE), arr.ind = TRUE)
  
  # 构建三列数据框
  jaccard_res <- data.frame(
    micro_1 = rownames(sim_matrix)[upper_indices[, 1]],
    micro_2 = colnames(sim_matrix)[upper_indices[, 2]],
    Co_occurrence = sim_matrix[upper_indices]
  )
  # plot
  # 边数据
  jaccard_edges <- as.data.frame(jaccard_res)
  colnames(jaccard_edges) <- c("from", "to", "co_value")
  
  # 创建图形对象
  jaccard_graph <- tbl_graph(
    nodes = sample_data_temp,
    edges = jaccard_edges,
    directed = FALSE
  ) %>% 
    activate(nodes) %>% 
    arrange(desc(Relative_Abundance))
  
  jaccard_plot <- ggraph(jaccard_graph, layout = "linear", circular = TRUE) +
    geom_edge_fan(
      aes(color = co_value),
      alpha = 0.4,
      edge_width = 0.5
    ) +
    scale_edge_color_gradientn(
      colors = c("#d7191c")
    ) +
    geom_node_point(
      aes(size = Relative_Abundance),
      color = "grey30"
    ) +
    scale_size_continuous(range = c(2, 8)) +
    labs(title = "Jaccard") +
    theme_void() +
    coord_fixed() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
  jaccard_plot
# Faith ((a + (c/n))/n)
  faith_res <- jaccard_res
  
  #plot 
  # 边数据
  faith_edges <- as.data.frame(faith_res)
  colnames(faith_edges) <- c("from", "to", "co_value")
  
  # 创建图形对象
  faith_graph <- tbl_graph(
    nodes = sample_data_temp,
    edges = faith_edges,
    directed = FALSE
  ) %>% 
    activate(nodes) %>% 
    arrange(desc(Relative_Abundance))
  
  faith_plot <- ggraph(faith_graph, layout = "linear", circular = TRUE) +
    geom_edge_fan(
      aes(color = co_value),
      alpha = 0.4,
      edge_width = 0.5
    ) +
    scale_edge_color_gradientn(
      colors = c("#d7191c")
    ) +
    geom_node_point(
      aes(size = Relative_Abundance),
      color = "grey30"
    ) +
    scale_size_continuous(range = c(2, 8)) +
    labs(title = "Faith") +
    theme_void() +
    coord_fixed() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
  faith_plot

# Braycurtis totalsum
  cooccurrence_braycurtis <- matrix(0, nrow = nrow(sample_data), ncol = nrow(sample_data))
  rownames(cooccurrence_braycurtis) <- rownames(sample_data_tss)
  colnames(cooccurrence_braycurtis) <- rownames(sample_data_tss)
  for (i in 1:(nrow(sample_data)-1)) {
    for (j in (i+1):nrow(sample_data)) {
      cooccurrence_braycurtis[i, j] <- min(sample_data_tss[i], sample_data_tss[j])*2/(sample_data_tss[i]+sample_data_tss[j])
    }
  }
  # 提取上三角（不包含对角线）的索引
  upper_indices <- which(upper.tri(cooccurrence_braycurtis, diag = FALSE), arr.ind = TRUE)
  # 构建三列数据框
  braycurtis_res <- data.frame(
    micro_1 = rownames(cooccurrence_braycurtis)[upper_indices[, 1]],
    micro_2 = colnames(cooccurrence_braycurtis)[upper_indices[, 2]],
    Co_occurrence = cooccurrence_braycurtis[upper_indices]
  )
  
  # plot
  # 边数据
  braycurtis_edges <- as.data.frame(braycurtis_res)
  colnames(braycurtis_edges) <- c("from", "to", "co_value")
  
  # 创建图形对象
  braycurtis_graph <- tbl_graph(
    nodes = sample_data_temp,
    edges = braycurtis_edges,
    directed = FALSE
  ) %>% 
    activate(nodes) %>% 
    arrange(desc(Relative_Abundance))
  
  braycurtis_plot <- ggraph(braycurtis_graph, layout = "linear", circular = TRUE) +
    geom_edge_fan(
      aes(color = co_value),
      alpha = 0.7,
      edge_width = 0.5
    ) +
    scale_edge_color_gradientn(
      colours = c("#000080", "#2c7bb6", "#7bc8f6", "#ffffbf", "#d7191c"),
  
    ) +
    geom_node_point(
      aes(size = Relative_Abundance),
      color = "grey30"
    ) +
    scale_size_continuous(range = c(2, 8)) +
    labs(title = "Braycurtis totalsum") +
    theme_void() +
    coord_fixed() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
  braycurtis_plot


# Braycurtis percentile
  # 定义归一化函数
  rank_normalize <- function(x) {
    ranks <- rank(x, ties.method = "average")  # 处理并列排名（可选参数：min/max）
    normalized <- ranks / max(ranks)           # 缩放到 [0,1]
    return(normalized)
  }
  sample_data_pct <- apply(sample_data, 2, rank_normalize)
  
  cooccurrence_braycurtis_pct <- matrix(0, nrow = nrow(sample_data), ncol = nrow(sample_data))
  rownames(cooccurrence_braycurtis_pct) <- rownames(sample_data)  
  colnames(cooccurrence_braycurtis_pct) <- rownames(sample_data)  
  for (i in 1:(nrow(sample_data)-1)) {
    for (j in (i+1):nrow(sample_data)) {
      cooccurrence_braycurtis_pct[i, j] <- min(sample_data_pct[i], sample_data_pct[j])*2/(sample_data_pct[i]+sample_data_pct[j])
    }
  }  
  # 提取上三角（不包含对角线）的索引
  upper_indices <- which(upper.tri(cooccurrence_braycurtis_pct, diag = FALSE), arr.ind = TRUE)
  # 构建三列数据框
  braycurtis_pct_res <- data.frame(
    micro_1 = rownames(cooccurrence_braycurtis_pct)[upper_indices[, 1]],
    micro_2 = colnames(cooccurrence_braycurtis_pct)[upper_indices[, 2]],
    Co_occurrence = cooccurrence_braycurtis_pct[upper_indices]
  )
  
  # plot
  # 边数据
  braycurtis_pct_edges <- as.data.frame(braycurtis_pct_res)
  colnames(braycurtis_pct_edges) <- c("from", "to", "co_value")
  
  # 创建图形对象
  braycurtis_pct_graph <- tbl_graph(
    nodes = sample_data_temp,
    edges = braycurtis_pct_edges,
    directed = FALSE
  ) %>% 
    activate(nodes) %>% 
    arrange(desc(Relative_Abundance))
  
  braycurtis_pct_plot <- ggraph(braycurtis_pct_graph, layout = "linear", circular = TRUE) +
    geom_edge_fan(
      aes(color = co_value),
      alpha = 0.7,
      edge_width = 0.5
    ) +
    scale_edge_color_gradientn(
      colours = c("#000080", "#2c7bb6", "#7bc8f6", "#ffffbf", "#d7191c"),
    ) +
    geom_node_point(
      aes(size = Relative_Abundance),
      color = "grey30"
    ) +
    scale_size_continuous(range = c(2, 8)) +
    labs(title = "Braycurtis percentile") +
    theme_void() +
    coord_fixed() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
  braycurtis_pct_plot
 

  
# Abundance totalsum
  cooccurrence_Abundance_Tss <- matrix(0, nrow = nrow(sample_data), ncol = nrow(sample_data))
  rownames(cooccurrence_Abundance_Tss) <- rownames(sample_data_tss)
  colnames(cooccurrence_Abundance_Tss) <- rownames(sample_data_tss)  
  for (i in 1:(nrow(sample_data)-1)) {
    for (j in (i+1):nrow(sample_data)) {
      cooccurrence_Abundance_Tss[i, j] <- min(sample_data_tss[i], sample_data_tss[j])*(1 - abs(sample_data_tss[i] - sample_data_tss[j]))
    }
  }
  # 提取上三角（不包含对角线）的索引  
  upper_indices <- which(upper.tri(cooccurrence_Abundance_Tss, diag = FALSE), arr.ind = TRUE)
  # 构建三列数据框
  Abundance_Tss_res <- data.frame(
    micro_1 = rownames(cooccurrence_Abundance_Tss)[upper_indices[, 1]],
    micro_2 = colnames(cooccurrence_Abundance_Tss)[upper_indices[, 2]],
    Co_occurrence = cooccurrence_Abundance_Tss[upper_indices]
  )
  
  # plot
  # 边数据
  Abundance_Tss_edges <- as.data.frame(Abundance_Tss_res)
  colnames(Abundance_Tss_edges) <- c("from", "to", "co_value")
  # 创建图形对象
  Abundance_Tss_graph <- tbl_graph(
    nodes = sample_data_temp,
    edges = Abundance_Tss_edges,
    directed = FALSE
  ) %>% 
    activate(nodes) %>% 
    arrange(desc(Relative_Abundance))
  
  Abundance_Tss_plot <- ggraph(Abundance_Tss_graph, layout = "linear", circular = TRUE) +
    geom_edge_fan(
      aes(color = co_value),
      alpha = 0.7,
      edge_width = 0.5
    ) +
    scale_edge_color_gradientn(
      colours = c("#000080", "#2c7bb6", "#7bc8f6", "#ffffbf", "#d7191c"),
      values = scales::rescale(
        quantile(Abundance_Tss_res$Co_occurrence,  # 关键修改点
                 probs = c(0, 0.25, 0.5, 0.75, 1))  # 四分位点设置
      )
    ) +
    geom_node_point(
      aes(size = Relative_Abundance),
      color = "grey30"
    ) +
    scale_size_continuous(range = c(2, 8)) +
    labs(title = "Abundance totalsum") +
    theme_void() +
    coord_fixed() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
  Abundance_Tss_plot


# Abundance percentile
  cooccurrence_Abundance_Pct <- matrix(0, nrow = nrow(sample_data), ncol = nrow(sample_data))
  rownames(cooccurrence_Abundance_Pct) <- rownames(sample_data)
  colnames(cooccurrence_Abundance_Pct) <- rownames(sample_data)  
  for (i in 1:(nrow(sample_data)-1)) {
    for (j in (i+1):nrow(sample_data)) {
      cooccurrence_Abundance_Pct[i, j] <- min(sample_data_pct[i], sample_data_pct[j])*(1 - abs(sample_data_pct[i] - sample_data_pct[j]))
    }
  }
  # 提取上三角（不包含对角线）的索引  
  upper_indices <- which(upper.tri(cooccurrence_Abundance_Pct, diag = FALSE), arr.ind = TRUE)
  # 构建三列数据框
  Abundance_Pct_res <- data.frame(
    micro_1 = rownames(cooccurrence_Abundance_Pct)[upper_indices[, 1]],
    micro_2 = colnames(cooccurrence_Abundance_Pct)[upper_indices[, 2]],
    Co_occurrence = cooccurrence_Abundance_Pct[upper_indices]
  )
  
  # plot
  # 边数据
  Abundance_Pct_edges <- as.data.frame(Abundance_Pct_res)
  colnames(Abundance_Pct_edges) <- c("from", "to", "co_value")
  # 创建图形对象
  Abundance_Pct_graph <- tbl_graph(
    nodes = sample_data_temp,
    edges = Abundance_Pct_edges,
    directed = FALSE
  ) %>% 
    activate(nodes) %>% 
    arrange(desc(Relative_Abundance))
  
  
  Abundance_Pct_plot <- ggraph(Abundance_Pct_graph, layout = "linear", circular = TRUE) +
    geom_edge_fan(
      aes(color = co_value),
      alpha = 0.7,
      edge_width = 0.5
    ) +
    scale_edge_color_gradientn(
      colours = c("#000080", "#2c7bb6", "#7bc8f6", "#ffffbf", "#d7191c"),
      values = scales::rescale(
        quantile(Abundance_Pct_res$Co_occurrence,  # 关键修改点
                 probs = c(0, 0.25, 0.5, 0.75, 1))  # 四分位点设置
      )
    ) +
    geom_node_point(
      aes(size = Relative_Abundance),
      color = "grey30"
    ) +
    scale_size_continuous(range = c(2, 8)) +
    labs(title = "Abundance percentile") +
    theme_void() +
    coord_fixed() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
  Abundance_Pct_plot
  
  # 图例绘制
  

  

# 网络图绘制(cycle)
  library(ggplot2)
  library(ggraph)
  library(tidygraph)
  library(viridis)
  library(patchwork)
  
  combined_plot <- ggpubr::ggarrange(
    Russell_rao_plot,
    Russell_rao_weight_plot,
    jaccard_plot,
    faith_plot,
    braycurtis_plot,
    braycurtis_pct_plot,
    Abundance_Tss_plot,
    Abundance_Pct_plot,
    ncol = 4,
    nrow = 2
  )
  combined_plot
  
  # 假设数据中的变量范围（替换为实际数据）
  co_occurrence_values <- Abundance_Pct_res$Co_occurrence
  rel_abundance_range <- range(sample_data_tss)
  
  # 1. 生成颜色图例的虚拟图形（标注 Low/High）
  color_legend_plot <- ggplot() +
    geom_segment(
      aes(x = 0, xend = 1, y = 0, yend = 0, color = co_occurrence_values),
      data = data.frame(co_occurrence_values = co_occurrence_values)
    ) +
    scale_color_gradientn(
      name = NULL,  # 不显示图例标题
      colours = c("#000080", "#2c7bb6", "#7bc8f6", "#ffffbf", "#d7191c"),
      values = scales::rescale(
        quantile(co_occurrence_values, probs = c(0, 0.25, 0.5, 0.75, 1))
      ),
      breaks = quantile(co_occurrence_values, probs = c(0, 1)),  # 仅显示最小和最大值
      labels = c("Low", "High")  # 标签强制设为 Low 和 High
    ) +
    theme_void() +
    theme(
      legend.key.width = unit(1, "cm"),  # 调整颜色条宽度
      legend.text = element_text(size = 10)
    )
  
  # 提取颜色图例
  color_legend <- cowplot::get_legend(color_legend_plot)
  
  # 2. 生成点大小图例的虚拟图形（保持数值标签）
  size_legend_plot <- ggplot() +
    geom_point(
      aes(x = 0, y = 0, size = rel_abundance_range),
      data = data.frame(rel_abundance_range = rel_abundance_range),
      color = "grey30"
    ) +
    scale_size_continuous(
      name = "Rel_abundance",  # 不显示图例标题
      range = c(2, 8),
      breaks = scales::pretty_breaks(n = 3)  # 自动生成合理刻度（如 0, 50, 100）
    ) +
    theme_void() +
    theme(
      legend.text = element_text(size = 10)
    )
  
  # 提取大小图例
  size_legend <- cowplot::get_legend(size_legend_plot)
  
  # 3. 组合图例
  combined_legends <- gridExtra::grid.arrange(
    color_legend,
    size_legend,
    ncol = 1
  )
  
  # 显示图例
  grid::grid.draw(combined_legends)
  
  combined_plot_final <- combined_plot + combined_legends + plot_layout(widths = c(8, 1))
  combined_plot_final
  