library(tidyverse)
library(ggpubr)
library(readr)
library(biomformat)
library(cowplot)
library(pROC)  
library(ggpubr) 
library(ggsci)

# otu_distance <- read.csv("/home/dongbiao/word_embedding_microbiome/HGT/results/otu_distance.csv")
# otu_distance$group <- factor(otu_distance$group, levels = c("Phylogeny", "Co-Em"))
# 
# p1 <- otu_distance %>% ggplot(aes(x = distance)) +
#     geom_density()+
#     facet_wrap(~group, scales = "free")+
#     theme_bw(base_size=14) +
#     labs(fill = NULL, x = "Distance", y = "Density") +
#     theme(text = element_text(face = "bold"), # 设置所有文本为加粗
#           plot.title = element_text(hjust = 0.5), # 居中对齐标题
#           axis.title = element_text(face = "bold"), # 加粗坐标轴标题
#           axis.text = element_text(face = "bold"))

table_group <- read.csv("/home/dongbiao/word_embedding_microbiome/HGT/results/tax_group.csv")

# 将 tax 列转换为有序因子
table_group$tax <- factor(table_group$tax, levels = c("Phylum", "Class", "Order", "Family", "Genus"), ordered = TRUE)

# 计算 Pearson 相关系数和趋势线
table_group <- table_group %>%
    group_by(tax) %>%
    mutate(
        r = cor(embed_cosine_mean, phy_dist_mean),  # 计算 Pearson 相关系数
        slope = lm(phy_dist_mean ~ embed_cosine_mean)$coefficients[2],  # 计算回归斜率
        intercept = lm(phy_dist_mean ~ embed_cosine_mean)$coefficients[1]  # 计算回归截距
    )

# 使用 ggplot2 绘制分面图
p1 <- ggplot(table_group, aes(x = embed_cosine_mean, y = phy_dist_mean)) +
    # 添加误差条
    geom_errorbar(aes(ymin = phy_dist_mean - phy_dist_sd, ymax = phy_dist_mean + phy_dist_sd), width = 0.02, color = "gray") +
    geom_errorbarh(aes(xmin = embed_cosine_mean - embed_cosine_sd, xmax = embed_cosine_mean + embed_cosine_sd), height = 0.02, color = "gray") +
    
    # 绘制散点图
    geom_point(size = 2, alpha = 0.4, color = "blue") +
    # 添加趋势线
    geom_abline(aes(slope = slope, intercept = intercept), color = "red", linetype = "dashed", linewidth = 0.8) +
    # 添加 Pearson 相关系数文本
    geom_text(aes(x = -0.1, y = 1, label = paste0("r = ", round(r, 2))), size = 4, color = "black") +
    # 分面显示
    facet_wrap(~ tax, nrow = 1) +
    # 添加标签和标题
    labs(
        x = "Embedding cosine",
        y = "Evolutionary distance",
        title = ""
    ) +
    # 美化主题
    theme_bw() +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          legend.position = "top",
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))

hgt_plot_res <- read.csv("/home/dongbiao/word_embedding_microbiome/HGT/results/hgt_plot_res.csv")

p2 <- hgt_plot_res %>% filter(group == "Phylogeny") %>% ggplot(aes(x = distance, y = hgt_rate)) +  
    geom_point() +  # 绘制散点图  
    geom_smooth(method = "loess", se = TRUE, color = "blue") +  # 添加线性回归趋势线  
    stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),   
             method = "pearson", label.x.npc = "left", label.y.npc = "top") +  # 标注R²和p值  
    theme_bw() +  # 使用简洁主题  
    labs(title = "", x = ("16S rRNA distance"), y = "HGT per 100 comparisons", color="")  +
    scale_color_brewer(palette = "Set1") + 
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          legend.position = "top",
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))

p3 <- hgt_plot_res %>% filter(group == "Co embedding") %>% ggplot(aes(x = distance, y = hgt_rate)) +  
    geom_point() +  # 绘制散点图  
    geom_smooth(method = "loess", se = TRUE, color = "blue") +  # 添加线性回归趋势线  
    stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),   
             method = "pearson", label.x.npc = "left", label.y.npc = "top") +  # 标注R²和p值  
    theme_bw() +  # 使用简洁主题  
    labs(title = "", x = ("Embedding cosine"), y = "HGT per 100 comparisons", color="")  +
    scale_color_brewer(palette = "Set1") + 
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          legend.position = "top",
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))

### with and without HGT
hgt_hight_low_sim <- read.csv("/home/dongbiao/word_embedding_microbiome/HGT/results/hgt_hight_low_sim.csv")
compire <- list(c('High EmbSim','Low EmbSim'))

p4 <- hgt_hight_low_sim %>%
    ggplot(aes(x = group, y = cosine_co, group = group)) +
    geom_jitter() + 
    theme_bw() + labs(x="", y="Clusters_111_abundance") +
    geom_signif(comparisons = compire, step_increase=0.5, map_signif_level = T, 
                test = wilcox.test, margin_top = 0.1, vjust = 1.5) +
    theme_bw(base_size=14) +
    labs(fill = "", x = "", y = "Embedding cosine", title = "") +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          legend.position = "None",
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))

p5 <- hgt_hight_low_sim %>%
    ggplot(aes(x = group, y = phy_dis, group = group)) +
    geom_jitter() + 
    theme_bw() + labs(x="", y="Clusters_111_abundance") +
    geom_signif(comparisons = compire, step_increase=0.5, map_signif_level = T, 
                test = wilcox.test, margin_top = 0.1, vjust = 1.5) +
    theme_bw(base_size=14) +
    labs(fill = "", x = "", y = "16S rRNA distance", title = "") +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          legend.position = "None",
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))

hgt_hight_sim <- hgt_hight_low_sim %>% filter(group == "High EmbSim")
hgt_low_sim <- hgt_hight_low_sim %>% filter(group == "Low EmbSim")

p_hgt_group <- plot_grid(p4, p5,
                         align="hv", labels = c('', ''),
                         nrow = 2, ncol=1, plot=FALSE, rel_heights = c(1, 1))
HGT_rate <- data.frame(value = c(sum(hgt_hight_sim$hgt > 0) / nrow(hgt_hight_sim) * 100, 
                                 sum(hgt_low_sim$hgt > 0) / nrow(hgt_low_sim) * 100),
                       group = c('High EmbSim','Low EmbSim'))


p6 <- HGT_rate %>% ggplot(aes(x = group, y = value)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw() + labs(x="", y="HGT per 100 comparisons", title = "") +
    theme_bw(base_size=14) +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          legend.position = "None",
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))


### HGT predict
df <- read.csv("/home/dongbiao/word_embedding_microbiome/HGT/results/hgt_predict_res.csv")
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

# 绘制图形  
p7 <- ggplot(summary_data, aes(x = fpr, y = mean_tpr, color = group)) +  
    geom_ribbon(aes(ymin = pmax(mean_tpr - sd_tpr, 0),  
                    ymax = pmin(mean_tpr + sd_tpr, 1),  
                    fill = group),  
                alpha = 0.2) +  
    geom_line(size = 1) +  
    geom_abline(slope = 1, intercept = 0,   
                linetype = "dashed", color = "gray50") +  
    scale_color_manual(values = c("#E69F00", "#56B4E9")) +  
    scale_fill_manual(values = c("#E69F00", "#56B4E9")) +  
    geom_text(data = auc_summary,  
              aes(x = 0.1, y = 0.1 - 0.07*(as.numeric(factor(group))-1),  
                  label = label),  
              hjust = 0, size = 4.5, show.legend = FALSE) +  
    labs(title = "HGT predict",  
         x = "False Positive Rate",  
         y = "True Positive Rate",  
         color = "",  
         fill = "") +  
    theme_bw(base_size = 14) +  
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          legend.position = "None",
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))




buttom_plot <- plot_grid(p2, p3, p_hgt_group, p6, p7, align="hv", 
                         labels = c('b', 'c', 'd', 'e', 'f'),
                         nrow = 1, ncol=5, plot=FALSE, rel_heights = c(1, 1, 1))

p <- plot_grid(p1, buttom_plot,
               align="hv", labels = c('a', ''),
               nrow = 2, ncol=1, plot=FALSE, rel_heights = c(1, 1))

ggsave("/home/dongbiao/word_embedding_microbiome/result/hgt.png", p,
       width = 38, height = 19, units = "cm")
