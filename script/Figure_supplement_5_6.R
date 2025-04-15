# 加载所需库
library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(aplot)
library(patchwork)
library(gridExtra)
library(grid)

#获取min-max值
# 初始化 min 和 max 值
min_data_fc <- Inf
max_data_fc <- -Inf
min_data_model <- Inf
max_data_model <- -Inf
min_data_model_log <- Inf
max_data_model_log <- -Inf


# 读取 CRC 和 IBD 文件夹下所有 CSV 文件
for (folder in c("CRC", "IBD")) {
    for (i in 1:7) {
        
        if (folder == "CRC"){
            file_path = paste0("/home/dongbiao/word_embedding_microbiome/programe_test/CRC/explain_model")
        }
        else{
            file_path = paste0("/home/dongbiao/word_embedding_microbiome/programe_test/IBD/explain")
        }
        
        # 读取 fc 文件
        data_fc <- read.csv(paste0(file_path, "/fc_", i, ".csv")) %>% column_to_rownames("X")
        min_data_fc <- min(min_data_fc, min(data_fc, na.rm = TRUE))
        max_data_fc <- max(max_data_fc, max(data_fc, na.rm = TRUE))
        
        # 读取 model 文件
        data_model <- read.csv(paste0(file_path, "/Model_", i, ".csv")) %>% column_to_rownames("X")
        min_data_model <- min(min_data_model, min(data_model, na.rm = TRUE))
        max_data_model <- max(max_data_model, max(data_model, na.rm = TRUE))
        
        #model_log
        # data_model <- log10(data_model+min(data_model[data_model != 0]))
        min_data_model_log <- min(min_data_model_log, min(data_model, na.rm = TRUE))
        max_data_model_log <- max(max_data_model_log, max(data_model, na.rm = TRUE))
      }
}

# 输出最终的min和max值

print(min_data_fc)
print(max_data_fc)
print(min_data_model)
print(max_data_model)


plot_heatmap <- function(data, metadata, title, xlab_text, fill_colors, min_data, max_data, col_cluster = TRUE,row_cluster = TRUE) {
  # 数据处理
  sample_names <- rownames(data)
  dimension_names <- colnames(data)
  
  # 层次聚类
  if (row_cluster) {
    cluster_sample_order <- hclust(dist(data, method = "euclidean"), method = "complete")$order
  } else {
    cluster_sample_order <- 1:nrow(data)
  }
  
  if (col_cluster) {
    cluster_dimension_order <- hclust(dist(t(data), method = "euclidean"), method = "complete")$order
  } else {
    cluster_dimension_order <- 1:ncol(data)
  }
  
  # 数据转换为长格式
  melted_data <- data %>%
    rownames_to_column("sample") %>%
    melt(id.vars = "sample") %>%
    mutate(
      sample = factor(sample, levels = sample_names[cluster_sample_order]),
      variable = factor(variable, levels = dimension_names[cluster_dimension_order])
    )
  
  # 创建热图
  p_heatmap <- ggplot(melted_data, aes(x = variable, y = sample)) +
    theme_classic() +
    scale_fill_gradientn(colors = fill_colors, limits = c(min_data, max_data)) +
    geom_tile(aes(fill = value), show.legend = FALSE) +
    labs(x = xlab_text, y = NULL, fill = NULL) +
    ggtitle(title) +
    theme(
      text = element_text(face = "bold"),
      axis.title = element_text(size = 12, color = "black", face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line = element_blank(),
      plot.title = element_text(hjust = 0.5)  # 标题居中
    )
  
  return(p_heatmap)
}

plot_row_annotation <- function(data, metadata, fill_values) {
  # 数据处理
  annotation_data <- data %>%
    rownames_to_column("sample") %>%
    mutate(Group = metadata$Group[match(sample, rownames(metadata))])
  
  # 创建行注释图
  p_row_annotation <- ggplot(annotation_data, aes(x = 1, y = sample)) +
    theme_classic() +
    geom_tile(aes(fill = Group), show.legend = FALSE) +
    scale_fill_manual(values = fill_values) +
    labs(y = "Samples", x = NULL, fill = NULL) +
    theme(
      text = element_text(face = "bold"),
      axis.title.y = element_text(size = 12, color = "black", face = "bold",hjust = 0.5, vjust = 0),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line = element_blank()
    )
  
  return(p_row_annotation)
}

# CRC函数
CRC_fc <- function(fc_file, metadata_file, title) {
  data <- read.csv(fc_file) %>% column_to_rownames("X")
  metadata <- read.csv(metadata_file, sep = "\t") %>%
    column_to_rownames("Run") %>%
    mutate(Group = if_else(group == 1, "CRC", "Control"))
  heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)
  annotation_colors <- c("Control" = "#2878B5", "CRC" = "#C82423")
  
  p_heatmap <- plot_heatmap(data, metadata, title, "Attention pooling", heatmap_colors,min_data = min_data_fc,max_data = max_data_fc)
  p_row_annotation <- plot_row_annotation(data, metadata, annotation_colors)
  
  p_combined <- p_heatmap %>% insert_left(p_row_annotation, width = 0.05)
  
  return(p_combined)
}

# CRC模型函数
CRC_model <- function(model_file, metadata_file, title) {
  data <- read.csv(model_file) %>% column_to_rownames("X")
  metadata <- read.csv(metadata_file, sep = "\t") %>%
    column_to_rownames("Run") %>%
    mutate(Group = if_else(group == 1, "CRC", "Control"))
  my_color_palette <- colorRampPalette(brewer.pal(9, "OrRd"))
  heatmap_colors <- my_color_palette(100)
  annotation_colors <- c("Control" = "#2878B5", "CRC" = "#C82423")
  
  p_heatmap <- plot_heatmap(data, metadata, title, "Attention weight", heatmap_colors,min_data = min_data_model,max_data = max_data_model,col_cluster = FALSE)
  p_row_annotation <- plot_row_annotation(data, metadata, annotation_colors)
  
  p_combined <- p_heatmap %>% insert_left(p_row_annotation, width = 0.05)
  
  return(p_combined)
}

# CRC模型函数_log
CRC_model_log <- function(model_file, metadata_file, title) {
  data <- read.csv(model_file) %>% column_to_rownames("X")
  data <- log10(data + min(data[data!=0]))
  metadata <- read.csv(metadata_file, sep = "\t") %>%
    column_to_rownames("Run") %>%
    mutate(Group = if_else(group == 1, "CRC", "Control"))
  my_color_palette <- colorRampPalette(brewer.pal(9, "OrRd"))
  heatmap_colors <- my_color_palette(100)
  annotation_colors <- c("Control" = "#2878B5", "CRC" = "#C82423")
  
  p_heatmap <- plot_heatmap(data, metadata, title, "Attention weight", heatmap_colors,min_data = min_data_model_log,max_data = max_data_model_log,col_cluster = FALSE)
  p_row_annotation <- plot_row_annotation(data, metadata, annotation_colors)
  
  p_combined <- p_heatmap %>% insert_left(p_row_annotation, width = 0.05)
  
  return(p_combined)
}

# IBD函数
IBD_fc <- function(fc_file, metadata_file, title) {
  data <- read.csv(fc_file) %>% column_to_rownames("X")
  metadata <- read.csv(metadata_file, sep = "\t") %>%
    column_to_rownames("sample") %>%
    mutate(Group = if_else(group == 1, "IBD", "Control"))
  
  heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)
  annotation_colors <- c("Control" = "#2878B5", "IBD" = "#FF6F61")
  p_heatmap <- plot_heatmap(data, metadata, title, "Attention pooling", heatmap_colors,min_data = min_data_fc,max_data = max_data_fc)
  p_row_annotation <- plot_row_annotation(data, metadata, annotation_colors)
  
  p_combined <- p_heatmap %>% insert_left(p_row_annotation, width = 0.05)
  
  return(p_combined)
}

# IBD模型函数
IBD_model <- function(model_file, metadata_file, title) {
  data <- read.csv(model_file) %>% column_to_rownames("X")
  metadata <- read.csv(metadata_file, sep = "\t") %>%
    column_to_rownames("sample") %>%
    mutate(Group = if_else(group == 1, "IBD", "Control"))
  my_color_palette <- colorRampPalette(brewer.pal(9, "OrRd"))
  heatmap_colors <- my_color_palette(100)
  annotation_colors <- c("Control" = "#2878B5", "IBD" = "#FF6F61")
  
  p_heatmap <- plot_heatmap(data, metadata, title, "Attention weight", heatmap_colors,min_data = min_data_model,max_data = max_data_model,col_cluster = FALSE)
  p_row_annotation <- plot_row_annotation(data, metadata, annotation_colors)
  
  p_combined <- p_heatmap %>% insert_left(p_row_annotation, width = 0.05)
  
  return(p_combined)
}

# IBD模型函数_log
IBD_model_log <- function(model_file, metadata_file, title) {
  data <- read.csv(model_file) %>% column_to_rownames("X")
  data <- log10(data + min(data[data!=0]))
  metadata <- read.csv(metadata_file, sep = "\t") %>%
    column_to_rownames("sample") %>%
    mutate(Group = if_else(group == 1, "IBD", "Control"))
  my_color_palette <- colorRampPalette(brewer.pal(9, "OrRd"))
  heatmap_colors <- my_color_palette(100)
  annotation_colors <- c("Control" = "#2878B5", "IBD" = "#FF6F61")
  
  p_heatmap <- plot_heatmap(data, metadata, title, "Attention weight", heatmap_colors,min_data = min_data_model_log,max_data = max_data_model_log,col_cluster = FALSE)
  p_row_annotation <- plot_row_annotation(data, metadata, annotation_colors)
  
  p_combined <- p_heatmap %>% insert_left(p_row_annotation, width = 0.05)
  
  return(p_combined)
}

# 画图
# 定义空列表来存储图
fc_plot_list <- list()
model_plot_list <- list()
model_plot_list_log <- list()

file_path = paste0("/home/dongbiao/word_embedding_microbiome/programe_test/IBD/explain")
for (i in 1:7) {
  #  IBD_fc 
  fc_plot_list[[i]] <- IBD_fc(paste0(file_path, "/fc_", i, ".csv"), "/home/dongbiao/word_embedding_microbiome/programe_test/IBD/data/metadata.txt", 
                              title = paste0("Model_", i)) %>% as.patchwork()
  #  IBD_model 
  # model_plot_list[[i]] <- IBD_model(paste0("IBD/Model_", i, ".csv"), "/home/dongbiao/word_embedding_microbiome/programe_test/IBD/data/metadata.txt", 
  #                                   title = paste0("Model_", i)) %>% as.patchwork()
  # #  IBD_model_log 
  # model_plot_list_log[[i]] <- IBD_model_log(paste0("IBD/Model_", i, ".csv"), "/home/dongbiao/word_embedding_microbiome/programe_test/IBD/data/metadata.txt", 
  #                                           title = paste0("Model_", i)) %>% as.patchwork()
}


for (i in 1:7) {
  # CRC_fc 
  fc_plot_list[[i + 7]] <- CRC_fc(paste0("CRC/fc_", i, ".csv"), "CRC/metadata.txt", title = paste0("Model_", i)) %>% as.patchwork()
  # CRC_model
  model_plot_list[[i + 7]] <- CRC_model(paste0("CRC/Model_", i, ".csv"), "CRC/metadata.txt", title = paste0("Model_", i)) %>% as.patchwork()
  # CRC_model_log
  model_plot_list_log[[i + 7]] <- CRC_model_log(paste0("CRC/Model_", i, ".csv"), "CRC/metadata.txt", title = paste0("Model_", i)) %>% as.patchwork()
}


# 创建图例
create_combined_legend <- function(min_data = -1, max_data = 1,heatmap_color) {
  # 创建分类图例
  create_color_legend <- function() {
    legend_data <- data.frame(
      Group = c("Control", "CRC", "IBD"),
      Color = c("#2878B5", "#C82423", "#FF6F61")
    )
    
    p_color_legend <- ggplot(legend_data, aes(x = Group, y = 1, fill = Group)) +
      geom_tile() +
      scale_fill_manual(values = legend_data$Color) +
      labs(fill = "Group") +
      theme_void() +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)
      )
    
    return(p_color_legend)
  }
  
  # 创建渐变图例
  create_gradient_legend <- function() {
    gradient_data <- data.frame(
      value = c(-1, 0, 1)
    )
    
    p_gradient_legend <- ggplot(gradient_data, aes(x = value, y = 1, fill = value)) +
      geom_tile() +
      scale_fill_gradientn(
        colors = heatmap_color,
        limits = c(min_data, max_data)
      ) +
      labs(fill = "Value") +
      theme_void() +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)
      )
    
    return(p_gradient_legend)
  }
  
  # 提取图例的辅助函数
  extract_legend <- function(plot) {
    g <- ggplotGrob(plot)
    legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
    return(legend)
  }
  
  # 创建并提取图例
  color_legend <- extract_legend(create_color_legend())
  gradient_legend <- extract_legend(create_gradient_legend())
  
  # 合并图例
  combined_legend <- grid.arrange(
    color_legend,
    gradient_legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") * 0.5, unit(1, "npc") * 0.5)
  )
  
  return(combined_legend)
}

# 调用函数显示合并的图例
my_color_palette <- colorRampPalette(brewer.pal(9, "OrRd"))

heatmap_colors <- my_color_palette(100)
heatmap_colors_fc <- colorRampPalette(c("blue", "white", "red"))(100)

legends_fc <- create_combined_legend(min_data = min_data_fc, max_data = max_data_fc, heatmap_color = heatmap_colors_fc)
legends_model <- create_combined_legend(min_data = min_data_model, max_data = max_data_model,heatmap_color = heatmap_colors)
legends_model_log <- create_combined_legend(min_data = min_data_model_log, max_data = max_data_model_log, heatmap_color = heatmap_colors)


# 定义 fc 部分的布局和宽度
fc_layout <- (fc_plot_list[[1]] | fc_plot_list[[2]] | fc_plot_list[[3]] | fc_plot_list[[4]]) / 
  (fc_plot_list[[5]] | fc_plot_list[[6]] | fc_plot_list[[7]] | plot_spacer()) / 
  (fc_plot_list[[8]] | fc_plot_list[[9]] | fc_plot_list[[10]] | fc_plot_list[[11]]) / 
  (fc_plot_list[[12]] | fc_plot_list[[13]] | fc_plot_list[[14]] | plot_spacer())


# 合并图和图例
combined_plot_fc <- fc_layout | legends_fc
combined_plot_fc <- combined_plot_fc + plot_layout(widths = c(6,1))

# 定义 model 部分的布局和宽度
model_layout <- (model_plot_list[[1]] | model_plot_list[[2]] | model_plot_list[[3]] | model_plot_list[[4]]) / 
  (model_plot_list[[5]] | model_plot_list[[6]] | model_plot_list[[7]] | plot_spacer()) / 
  (model_plot_list[[8]] | model_plot_list[[9]] | model_plot_list[[10]] | model_plot_list[[11]]) / 
  (model_plot_list[[12]] | model_plot_list[[13]] | model_plot_list[[14]] | plot_spacer())

combined_plot_model <- model_layout | legends_model

combined_plot_model <- combined_plot_model + plot_layout(
  widths = c(6,1)
)

model_layout_log <- (model_plot_list_log[[1]] | model_plot_list_log[[2]] | model_plot_list_log[[3]] | model_plot_list_log[[4]]) / 
  (model_plot_list_log[[5]] | model_plot_list_log[[6]] | model_plot_list_log[[7]] | plot_spacer()) / 
  (model_plot_list_log[[8]] | model_plot_list_log[[9]] | model_plot_list_log[[10]] | model_plot_list_log[[11]]) / 
  (model_plot_list_log[[12]] | model_plot_list_log[[13]] | model_plot_list_log[[14]] | plot_spacer())

combined_plot_model_log <- model_layout_log | legends_model_log

combined_plot_model_log <- combined_plot_model_log + p_color_legend


# 保存图
ggsave("heatmap_fc.png", combined_plot_fc, width = 8, height = 6, dpi = 300)
ggsave("heatmap_model.png", combined_plot_model, width = 8, height = 6, dpi = 300)
ggsave("heatmap_model_log.png", combined_plot_model_log, width = 8, height = 6, dpi = 300)

