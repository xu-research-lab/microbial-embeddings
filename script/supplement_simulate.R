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

simulation_data_ratio <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/simulation_data/results/simulation_data_ratio.csv")
simulation_data_ratio$group2 <- factor(simulation_data_ratio$group2, levels = c("abundance-percentile", "abundance-totalsum", "russell_rao_weight", "russell_rao", 
                                                                                "braycurtis-percentile", "braycurtis-totalsum", "jaccard", "faith"))

p1 <- simulation_data_ratio %>% ggplot(aes(x=log10(value), group=group1, color=group1)) +   
    geom_density(aes(y = ..density..), alpha = 0.5, size=1) +
    scale_color_brewer(palette = "Set2") +
    facet_wrap(group2 ~., scales = "free_y", nrow=2) +
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
dif_method_jasim$methed <- factor(dif_method_jasim$methed, levels = c("abundance-percentile", "abundance-totalsum", "russell_rao_weight", "russell_rao", 
                                                                      "braycurtis-percentile", "braycurtis-totalsum", "jaccard", "faith"))
p2 <- dif_method_jasim %>% ggplot(aes(x = jac_smi, y = cosine)) +
    geom_point() +
    labs(x="Jaccard similarity", y="Cosine") +
    facet_wrap(methed ~., scales = "free_y", nrow=2) +
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
                      nrow = 2, ncol=1, plot=FALSE, rel_heights = c(1, 1))

ggsave("/home/dongbiao/word_embedding_microbiome/result/supplement_simulation_res.png", upp_plot,
       width = 30, height = 20, units = "cm")
