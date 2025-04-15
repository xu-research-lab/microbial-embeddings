library(dplyr)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(lemon)
library(cowplot)
library(reshape2)
library(tidyverse)
library(stringr)

library(ggdensity)
library(ggblanket)
library(ggsci)

library(gridExtra)
library(patchwork)
library(cowplot)
library(aplot)
library(ggplotify)

simulation_data_ratio <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/simulation_data/results/simulation_data_ratio.csv")
simulation_data_ratio <- simulation_data_ratio %>% filter(group2 == "abundance-percentile")

p1 <- simulation_data_ratio %>% ggplot(aes(x=log10(value), group=group1, color=group1)) +   
    geom_density(aes(y = ..density..), alpha = 0.5, size=1) +
    scale_color_brewer(palette = "Set2") +
    labs(title = "Probability Density Function", x = "Value", y = "Density") +  
    theme_bw()+
    labs(x="log10(Ratios of co-occurrence probabilities)", y="Density", title="")+
    theme(text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 0),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5))

dif_method_jasim <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/simulation_data/results/diff_method_jasim_cos.txt", 
                             sep = "\t")
dif_method_jasim <- dif_method_jasim %>% filter(methed == "abundance-percentile")
p2 <- dif_method_jasim %>% ggplot(aes(x = jac_smi, y = cosine)) +
    geom_point() +
    labs(x="Jaccard similarity", y="Cosine") +
    # facet_wrap(methed ~., scales = "free_y", nrow=2) +
    stat_smooth(method = "lm") +
    stat_cor(label.x = 0.35, label.y.npc = "top", color = "red", r.digits = ) +
    theme_bw(base_size = 14) +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))
upp_plot <- plot_grid(p1, p2,
                         labels = c('a', 'b'),
                         align="hv",
                         scale = c(1, 1),
                         nrow = 1, ncol=2, plot=FALSE, rel_widths = c(1, 1))

diff_dim <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/simulation_data/results/diff_embedsize_jasim_cosine.txt", sep = "\t")
diff_dim$embed_size <- paste("Embedding size", as.character(diff_dim$embed_size))
diff_dim$embed_size <- factor(diff_dim$embed_size, levels = c("Embedding size 50", "Embedding size 100", "Embedding size 150", "Embedding size 200"))
p3 <- diff_dim %>% ggplot(aes(x = jac_smi, y = cosine)) +
    geom_point() +
    labs(x="Jaccard similarity", y="Cosine") +
    facet_wrap(embed_size ~., scales = "free_y", nrow=1) +
    stat_smooth(method = "lm") +
    stat_cor(label.x = 0.35,
             label.y.npc = "top", color = "red", r.digits = ) +
    theme_bw(base_size = 14) +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))

diff_datasets <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/simulation_data/results/diff_datasize_r2_cosine.txt", sep = "\t") %>%
    column_to_rownames("data_size") %>% t() %>% as.data.frame() %>% melt(variable.name="data_size", value.name = "R")
p4 <- diff_datasets %>% ggplot(aes(x=data_size, y=R)) +
    geom_boxplot() +
    geom_smooth(level= 0.9) +
    geom_point() +
    theme_bw(base_size = 14) +
    labs(x="Data size", y="R") +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))




### k_means with pathway
k_means_pathway_embed <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/niche_func/k_means_func_embed.csv")
p_kmeans_func_embed <- ggplot(k_means_pathway_embed, aes(x = embed_sim, y = fun_sim, group = group, fill=group)) +  
    geom_hdr(show.legend = FALSE, method = "kde", n=100,
             probs = c(0.85, 0.8, 0.7, 0.6, 0.5, 0.4, 0.2, 0.25, 0.15, 0.1, 0.05)) +
    scale_fill_manual(values = c("Niche-equivalent" = "blue", "Niche-non-equivalent" = "red")) +  
    labs(x = "Embedding Similarity", y = "Pathway Similarity") +  
    theme_bw(base_size=14) +
    theme(legend.position = "right", text = element_text(face = "bold"))
# Embedding and Function Similarity by Jaccard Similarity Level

distribut_embed <- k_means_pathway_embed %>% ggplot(aes(x = embed_sim)) +
    geom_density(aes(fill = group), alpha = 0.4)+
    theme_bw(base_size=14) +
    labs(fill = NULL, x = NULL, y = NULL) +
    scale_fill_manual(values = c("Niche-equivalent" = "blue", "Niche-non-equivalent" = "red")) + 
    scale_color_brewer(palette = "Paired") +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          legend.position = "None",
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"), 
          axis.ticks = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_blank())

distribut_func <- k_means_pathway_embed %>% ggplot(aes(y = fun_sim)) +
    geom_density(aes(fill = group), alpha = 0.4)+
    scale_fill_manual(values = c("Niche-equivalent" = "blue", "Niche-non-equivalent" = "red")) +
    # facet_grid(. ~ Type, scales="free") + 
    theme_bw(base_size=14) +
    labs(fill = NULL, x = NULL, y = NULL) +
    scale_color_brewer(palette = "Paired") +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          legend.position = "None",
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"),axis.ticks = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_blank())
p5 <- p_kmeans_func_embed %>% insert_right(distribut_func, width=0.2) %>% 
    insert_top(distribut_embed, height = 0.2)
p5 <- as.ggplot(p5)
bottom_plot <- plot_grid(p4, p5,
                         labels = c('d', 'e'),
                         align="hv",
                         scale = c(0.5, 1),
                         nrow = 1, ncol=2, plot=FALSE, rel_widths = c(1, 1))

p <- plot_grid(upp_plot, p3, bottom_plot,
               labels = c('', "c", ""),
               align="hv",
               scale = c(1, 1, 1),
               rel_heights = c(0.8, 0.8, 1),
               nrow = 3, ncol=1, plot=FALSE)

ggsave("/home/dongbiao/word_embedding_microbiome/result/simulation_res.png", p,
       width = 30, height = 30, units = "cm")
