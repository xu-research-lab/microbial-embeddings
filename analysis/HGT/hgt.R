library(tidyverse)
library(ggpubr)
library(readr)
library(biomformat)
library(cowplot)
library(pROC)  
library(scales)
library(ggrepel)
library(patchwork)
library(ggpubr) 
library(ggsci)


# 将 tax 列转换为有序因子
table_group <- read.csv("Data/tax_group.csv")
table_group$tax <- factor(table_group$tax, levels = c("Phylum", "Class", "Order", "Family", "Genus"), ordered = TRUE)

# 计算 Pearson 相关系数和趋势线
table_group <- table_group %>%
    group_by(tax) %>%
    mutate(
        r = cor(embed_cosine_mean, phy_dist_mean), # 计算 Pearson 相关系数
        slope = lm(phy_dist_mean ~ embed_cosine_mean)$coefficients[2], # 计算回归斜率
        intercept = lm(phy_dist_mean ~ embed_cosine_mean)$coefficients[1] # 计算回归截距
    )

highlight_data <- table_group %>%
    filter(tax == "Genus") %>%
    filter(name %in% c("Bilophila", "Burkholderia-Caballeronia-Paraburkholderia", "Elusimicrobium", "Tatumella")) %>%
    mutate(
        # 为Burkholderia-Caballeronia-Paraburkholderia创建简写
        label = ifelse(name == "Burkholderia-Caballeronia-Paraburkholderia", "BCP", name)
    )

# 使用 ggplot2 绘制分面图
p1 <- ggplot(table_group, aes(x = embed_cosine_mean, y = phy_dist_mean)) +
    # 添加常规误差条 - 保持原始线宽
    geom_errorbar(aes(ymin = phy_dist_mean - phy_dist_sd, ymax = phy_dist_mean + phy_dist_sd),
                  width = 0.02, color = "gray", linewidth = 0.3
    ) +
    geom_errorbarh(aes(xmin = embed_cosine_mean - embed_cosine_sd, xmax = embed_cosine_mean + embed_cosine_sd),
                   height = 0.02, color = "gray", linewidth = 0.3
    ) +
    
    # 绘制常规散点图
    geom_point(size = 2, alpha = 0.4, color = "blue") +
    
    # 为Genus分面添加高亮物种的误差条（红色） - 不再加粗
    geom_errorbar(
        data = highlight_data,
        aes(ymin = phy_dist_mean - phy_dist_sd, ymax = phy_dist_mean + phy_dist_sd),
        width = 0.02, color = "red", linewidth = 0.3
    ) +
    geom_errorbarh(
        data = highlight_data,
        aes(xmin = embed_cosine_mean - embed_cosine_sd, xmax = embed_cosine_mean + embed_cosine_sd),
        height = 0.02, color = "red", linewidth = 0.3
    ) +
    
    # 为高亮物种添加点（红色）
    geom_point(data = highlight_data, size = 4, color = "red", shape = 16) +
    
    # 使用ggrepel添加标签 - 增加标签距离
    geom_text_repel(
        data = highlight_data,
        aes(label = label),
        color = "red",
        size = 4,
        box.padding = 1.0, # 增加标签周围的空间（原来是0.5）
        point.padding = 0.5, # 增加标签与点之间的距离（原来是0.5）
        min.segment.length = 0, # 确保总是有连接线
        segment.color = "red", # 连接线颜色
        segment.size = 0.6, # 连接线粗细（略微加粗以提高可见性）
        nudge_x = 0.1, # 水平方向额外偏移
        nudge_y = 0.01, # 垂直方向额外偏移
        direction = "both"
    ) + # 允许标签在任何方向移动
    
    # 添加趋势线
    geom_abline(aes(slope = slope, intercept = intercept), color = "red", linetype = "dashed", linewidth = 0.8) +
    
    # 添加 Pearson 相关系数文本
    geom_text(aes(x = 0.2, y = 1.5, label = paste0("r = ", round(r, 2))), size = 4, color = "black") +
    
    # 分面显示
    facet_wrap(~tax, nrow = 1) +
    
    # 添加标签和标题
    labs(
        x = "SNEsim",
        y = "Phylo distance",
        title = ""
    ) +
    
    # 美化主题
    theme_bw() +
    theme(
        text = element_text(face = "bold"), # 设置所有文本为加粗
        legend.position = "top",
        plot.title = element_text(hjust = 0.5), # 居中对齐标题
        axis.title = element_text(face = "bold"), # 加粗坐标轴标题
        axis.text = element_text(face = "bold")
    )

hgt_plot_res <- read.csv("Data/hgt_plot_res.csv")
hgt_plot_res$distance[hgt_plot_res$group == "SNE"] <- 1 - hgt_plot_res$distance[hgt_plot_res$group == "SNE"]

df <- hgt_plot_res %>% filter(group == "Phylo")
summary(lm(hgt_rate ~ distance, df))
# 获取 Phylogeny 组数据的 x 范围（主x轴基准）
phylogeny_data <- hgt_plot_res %>% filter(group != "SNE")
x_primary_range <- range(phylogeny_data$distance)

# 将 Co embedding 的 x 值（cosine）线性映射到主x轴范围
coembed_data <- hgt_plot_res %>% filter(group == "SNE") %>% 
    mutate(scaled_x = rescale(distance, to = x_primary_range))

# 创建基础图层（Phylogeny 组）
# 定义颜色方案
# 假设coembed_data包含原始cosine距离列
coembed_data$scaled_x <- rescale(coembed_data$distance, 
                                 to = range(phylogeny_data$distance))
primary_range <- range(phylogeny_data$distance)
secondary_range <- range(1 - coembed_data$distance)  # 假设需要显示1 - cosine
phylogeny_data$group <- factor(phylogeny_data$group, levels = c("Phylo", "Phylo SNEsim > 0.5", "Phylo SNEsim < 0.5"))
Phylo_paired <- phylogeny_data %>%
    dplyr::filter(group != "Phylo")

group1_data <- Phylo_paired$hgt_rate[Phylo_paired$group == "Phylo SNEsim > 0.5"]
group2_data <- Phylo_paired$hgt_rate[Phylo_paired$group == "Phylo SNEsim < 0.5"]

paired_test_res <- t.test(group2_data, group1_data, paired = TRUE)


# --- 主绘图对象 ---
p2 <- ggplot() +
    # Phylogeny数据层（三个子组）
    geom_point(
        data = phylogeny_data,
        aes(x = distance, y = hgt_rate, color = group, shape = group),
        size = 2
    ) +
    geom_smooth(
        data = phylogeny_data,
        aes(x = distance, y = hgt_rate, color = group, linetype = group),
        method = "loess", se = TRUE, size = 0.6
    ) +
    # 修正Phylo组的R²标签
    ggpubr::stat_cor(
        data = phylogeny_data[phylogeny_data$group == "Phylo", ],
        aes(x = distance, y = hgt_rate, label = ..rr.label..),
        label.x.npc = 0.15,  
        label.y.npc = 0.7,    # 标签y轴位置
        color = "#E41A1C",      # 强制指定颜色与Phylo组一致
        size = 3,
        show.legend = FALSE
    ) +
    # Co-embedding数据层
    geom_point(
        data = coembed_data,
        aes(x = scaled_x, y = hgt_rate, color = "SNE"),
        shape = 16, size = 2  
    ) +
    geom_smooth(
        data = coembed_data,
        aes(x = scaled_x, y = hgt_rate, color = "SNE"),
        method = "loess", se = TRUE, size = 0.6
    ) +
    # 修正SNE组的R²标签
    ggpubr::stat_cor(
        data = coembed_data,
        aes(x = scaled_x, y = hgt_rate, label = ..rr.label..),
        label.x.npc = 0.2,
        label.y.npc = 0.3,
        color = "#377EB8",       # 强制指定颜色与SNE组一致
        size = 3,
        show.legend = FALSE
    ) +
    
    # --- 双轴系统配置 ---
    scale_x_continuous(
        name = "Phylo distance",
        limits = range(phylogeny_data$distance),
        sec.axis = sec_axis(
            ~ rescale(., to = range(1 - coembed_data$distance)),
            name = "1 - SNESim",
            breaks = pretty_breaks(n = 5)
        )
    ) +
    
    # --- 视觉映射系统 ---
    scale_color_manual(
        name = "Group",
        values = c(
            "Phylo" = "#E41A1C",
            "Phylo SNEsim > 0.5" = "#E41A1C",
            "Phylo SNEsim < 0.5" = "#E41A1C",
            "SNE" = "#377EB8"
        ),
        # labels = c(
        #     "Phylogeny" = "Phylogeny (Baseline)",
        #     "Phylogeny_high" = "Phylogeny (High)",
        #     "Phylogeny_low" = "Phylogeny (Low)",
        #     "Co-embedding" = "Co-embedding"
        # )
    ) +
    scale_shape_manual(
        name = "Group",
        values = c(
            "Phylo" = 16,      # 圆形
            "Phylo SNEsim > 0.5" = 17, # 三角形
            "Phylo SNEsim < 0.5" = 15,  # 正方形
            "SNE" = 16    # 菱形
        ),
        guide = "none"  # 隐藏独立线型图例
    ) +
    scale_linetype_manual(
        name = "Group",
        values = c(
            "Phylo" = "solid",
            "Phylo SNEsim > 0.5" = "dashed",
            "Phylo SNEsim < 0.5" = "dotted",
            "SNE" = "solid"
        ),
        guide = "none"  # 隐藏独立线型图例
    ) +
    
    # --- 图例合并系统 ---
    guides(
        color = guide_legend(
            override.aes = list(
                linetype = c("solid", "dashed", "dotted", "solid"),
                shape = c(16, 17, 15, 16)
            )
        )
    ) +
    geom_text(
        data = data.frame(
            x = 0.8, 
            y = 28, 
            label = paste0("p-value: ", 
                           format(paired_test_res$p.value, format = "e", digits = 3, scientific=TRUE))
        ),
        aes(x = x, y = y, label = label),
        size = 3, 
        color = "black"
    )

# --- 主题优化 ---
p2 <- p2 + theme_bw(base_size = 12) +
    theme(
        legend.title = element_blank(),
        axis.title.x.bottom = element_text(color = "#E41A1C", face = "bold"),
        axis.text.x.bottom = element_text(color = "#E41A1C"),
        axis.title.x.top = element_text(color = "#377EB8", face = "bold"),
        axis.text.x.top = element_text(color = "#377EB8")
    ) +
    labs(y = "HGT Per 100 Comparisons")


### HGT predict
df <- read.csv("Data/hgt_predict_res_all.csv")
fpr_grid <- seq(0, 1, length.out = 100)  
# 计算每折的 ROC 曲线，并在固定 FPR 网格下插值出 TPR，同时计算 AUC  
roc_data <- df %>%   
    group_by(group, test) %>%             
    summarise(  
        roc_obj = list(roc(labels, proba, quiet = TRUE)),  
        auc = as.vector(roc(labels, proba, quiet = TRUE)$auc),  
        # 利用 pROC::coords 在固定 FPR 网格下获得 sensitivity (即 TPR)  
        tpr = coords(roc_obj[[1]], x = fpr_grid, input = "fpr", ret = "sensitivity", transpose = FALSE),  
        fpr = fpr_grid,  
        .groups = "drop"  
    )
roc_data <- roc_data[,c(1,2,4,5,6)]
roc_data <- roc_data %>%  
    unnest(cols = c(fpr, tpr))  
colnames(roc_data) <- c("group", "test", "auc", "tpr", "fpr")

# 计算统计量  
summary_data <- roc_data %>%  
    group_by(group, fpr) %>%  
    summarise(  
        mean_tpr = mean(tpr, na.rm = TRUE),  
        sd_tpr = sd(tpr, na.rm = TRUE),  
        .groups = "drop"  
    )  

auc_summary <- df %>%  
    group_by(test, group) %>%  
    summarise(auc = as.vector(roc(labels, proba)$auc), .groups = "drop") %>%  
    group_by(group) %>%  
    summarise(  
        mean_auc = mean(auc),  
        sd_auc = sd(auc),  
        label = sprintf("%s: %.3f ± %.3f",   
                        group, mean_auc, sd_auc)  
    )  

auc_data <- roc_data
auc_data$group1 <- paste0(auc_data$group, "_", auc_data$test)
auc_data <- auc_data[!duplicated(auc_data$group1), ]

p <- auc_data %>% ggplot(aes(x=group, y=auc, group=test)) +
    geom_point(size=3) +
    geom_line() +
    theme_bw(base_size = 14)+
    theme(text = element_text(face = "bold"), 
          legend.text = element_text(size = 14, color = "black", face = "bold"),
          legend.box.margin = margin(l = 0.1, unit = "in"),
          legend.title.align = 0.5 
    ) +
    labs(x="", y="AUC")