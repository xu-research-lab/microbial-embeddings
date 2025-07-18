library(tidyverse)
library(cowplot)
library(reshape2)
library(ggpubr)

MDS_ibd <- read.csv("../../Data/biomark/MDS_IBD_shap.csv")
MDS_crc <- read.csv("../../Data/biomark/MDS_CRC_shap.csv")

MDS_ibd$emb_type <- factor(MDS_ibd$emb_type, levels=c("SNE", "Phylo"))
MDS_crc$emb_type <- factor(MDS_crc$emb_type, levels=c("SNE", "Phylo"))
MDS_ibd$group <- factor(MDS_ibd$group, levels = c("Ctrl", "IBD"))
MDS_crc$group <- factor(MDS_crc$group, levels = c("Ctrl", "CRC"))
p1 <- ggplot(MDS_ibd, aes(x = t_SNE_1, y = t_SNE_2, color = group)) + 
    geom_point(size=3, alpha=0.5) + 
    facet_wrap(~emb_type, scales = "free") +
    labs(x = "MDS_1", y = "MDS_2") +
    theme_bw(base_size = 14) + 
    # scale_color_brewer(palette = "Set2")+
    theme(  
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14),
        text = element_text(face = "bold"),
        legend.title = element_blank(), 
        plot.margin = margin(0, 0, 0, 0),  # 调整边距
        legend.position = "top"  # 设置图例位置为左下角
    ) +  
    scale_shape_manual(values = c(16, 1)) 

p2 <- ggplot(MDS_crc, aes(x = t_SNE_1, y = t_SNE_2, color = group)) + 
    geom_point(size=3, alpha=0.5) + 
    facet_wrap(~emb_type, scales = "free") +
    labs(x = "MDS_1", y = "MDS_2") +
    theme_bw(base_size = 14) + 
    # scale_color_brewer(palette = "Set2")+
    theme(  
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14),
        text = element_text(face = "bold"),
        legend.title = element_blank(), 
        plot.margin = margin(0, 0, 0, 0),  # 调整边距
        legend.position = "top"  # 设置图例位置为左下角
    ) +  
    scale_shape_manual(values = c(16, 1)) 

p <- plot_grid(p1, p2, labels = c('IBD', 'CRC'),
               align="hv",
               nrow = 1, ncol=2, plot=FALSE, rel_heights = c(1, 1))
