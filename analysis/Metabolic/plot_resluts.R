library(dplyr)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(lemon)
library(cowplot)
library(reshape2)
library(tidyverse)
library(stringr)
library(ggbeeswarm)

library(ggdensity)
library(ggblanket)
library(ggsci)
library(reshape2)
library(gridExtra)
library(patchwork)
library(cowplot)
library(aplot)
library(ggplotify)

library(RColorBrewer)  

set1_colors <- brewer.pal(n = 9, name = "Set1")  

high_sim_res <- read.csv("Data/high_sim_res_M11.csv")
high_sim_res <- high_sim_res %>% filter(cosine > 0.9)
high_sim_res$group <- rep("cosine > 0.9", nrow(high_sim_res))
low_SNE_res <- read.csv("Data/low_SNE_res.csv", row.names = 1)
low_SNE_res$group <- rep("cosine ≈ 0", nrow(low_SNE_res))
high_sim_res <- high_sim_res[, colnames(low_SNE_res)]
plot_df <- rbind(high_sim_res, low_SNE_res)

p_mro <- t.test(mro ~ group, data = plot_df)$p.value
p_mip <- t.test(mip ~ group, data = plot_df)$p.value

p1 <- plot_df %>% ggplot(aes(x = group, y = mro, fill = group)) +
    geom_boxplot(width = 0.4, alpha = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("cosine > 0.9" = "#E41A1C", "cosine ≈ 0" = "#377EB8")) +
    theme_bw(base_size = 14) +
    labs(fill = " ", x = NULL, y = "Competition MRO") +
    scale_y_continuous(
        limits = c(0.5, 1), breaks = seq(0.5, 1, by = 0.1)) +
    theme(
        text = element_text(face = "bold"),
        legend.position = "none"
    ) +
    annotate("text", x = -Inf, y = Inf, 
             label = ifelse(!is.na(p_mro), sprintf("p=%.2g", p_mro), ""),
             color = "black", hjust = -0.05, vjust = 1.1, size = 4) +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1))
p2 <- plot_df %>% ggplot(aes(x = group, y = mip, fill = group)) +
    geom_boxplot(width = 0.4, alpha = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = c("cosine > 0.9" = "#E41A1C", "cosine ≈ 0" = "#377EB8")) +
    theme_bw(base_size = 14) +
    labs(fill = " ", x = NULL, y = "Competition MIP") +
    scale_y_continuous(
        limits = c(0, 8), breaks = seq(0, 8, by = 2)) +
    theme(
        text = element_text(face = "bold"),
        legend.position = "none"
    ) +
    annotate("text", x = -Inf, y = Inf, 
             label = ifelse(!is.na(p_mip), sprintf("p=%.2g", p_mip), ""),
             color = "black", hjust = -0.05, vjust = 1.1, size = 4) +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1))

p_all <- plot_grid(p1, p2, align="hv",
               nrow = 1, ncol=4, plot=FALSE, rel_heights = c(1, 1))


high_sim_res <- read.csv("Data/high_sim_res_M11.csv")
high_sim_res <- high_sim_res %>% filter(cosine > 0.9) %>% arrange(co_occur)
high_sim_res_low <- high_sim_res[1:100, ]
high_sim_res_low$group <- rep("Low co-occur", 100)
high_sim_res_high <- tail(high_sim_res, 100)
high_sim_res_high$group <- rep("High co-occur", 100)
high_sim_res_plot <- rbind(high_sim_res_high, high_sim_res_low)
p_mro <- t.test(mro ~ group, data = high_sim_res_plot)$p.value
p_mip <- t.test(mip ~ group, data = high_sim_res_plot)$p.value
high_sim_res_plot$mip_random <- high_sim_res_plot$mip + runif(length(high_sim_res_plot$mip), min = -0.2, max = 0.2)
p_dotplot <- ggplot(high_sim_res_plot, aes(x = mro, y = mip_random, group = group, color=group)) +  
    geom_point(alpha=1, size=1)+
    scale_color_brewer(palette = "Set2") +
    labs(x = "Competition MRO", y = "Cooperation MIP", color="") +  
    theme_bw(base_size=14) +
    scale_x_continuous(
        limits = c(0.5, 1), breaks = seq(0.5, 1, by = 0.1)) + 
    scale_y_continuous(
        limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
    theme(legend.position = "None",
          text = element_text(face = "bold"))
distribut_mro <- high_sim_res_plot %>%ggplot(aes(x = mro, y = group, color = group)) +
    geom_boxplot(width = 0.3, alpha = 0.5, outlier.shape = NA) +
    scale_color_brewer(palette = "Set2") +
    theme_bw(base_size=14) +
    labs(fill = NULL, x = NULL, y = NULL) +
    scale_x_continuous(
        limits = c(0.4, 1), breaks = seq(0.4, 1, by = 0.1)) + 
    theme(text = element_text(face = "bold"),
          legend.position = "None",
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    annotate("text", x = -Inf, y = Inf, 
             label = ifelse(!is.na(p_mro), sprintf("p=%.2g", p_mro), ""),
             color = "black", hjust = -0.05, vjust = 1.1, size = 4)
distribut_mip <- high_sim_res_plot %>% ggplot(aes(y = mip, x = group, color = group)) +
    geom_boxplot(width = 0.3, alpha = 0.5, outlier.shape = NA) +
    scale_color_brewer(palette = "Set2") +
    theme_bw(base_size=14) +
    labs(fill = NULL, x = NULL, y = NULL) +
    scale_y_continuous(
        limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
    theme(text = element_text(face = "bold"),
          legend.position = "None",
          axis.text = element_blank(),
          axis.ticks = element_blank())+
    annotate("text", x = Inf, y = Inf,  
             label = ifelse(!is.na(p_mip), sprintf("p=%.2g", p_mip), ""),
             color = "black", 
             hjust = -0.09,  
             vjust = 1, 
             angle = 270,  
             size = 4)
p1 <- p_dotplot %>% 
    insert_right(distribut_mip, width = 0.2) %>% 
    insert_top(distribut_mro, height = 0.2)
p1 <- as.ggplot(p1) + 
    labs(title = "cosine > 0.9") +  
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

low_sim_res <- read.csv("Data/low_SNE_res.csv", row.names = 1)
low_sim_res <- low_sim_res %>% arrange(co_occur)
low_sim_res_low <- low_sim_res[1:100, ]
low_sim_res_low$group <- rep("Low co-occur", 100)
low_sim_res_high <- tail(low_sim_res, 100)
low_sim_res_high$group <- rep("High co-occur", 100)
low_sim_res_plot <- rbind(low_sim_res_high, low_sim_res_low)
p_mro <- t.test(mro ~ group, data = low_sim_res_plot)$p.value
p_mip <- t.test(mip ~ group, data = low_sim_res_plot)$p.value
low_sim_res_plot$mip_random <- low_sim_res_plot$mip + runif(length(low_sim_res_plot$mip), min = -0.2, max = 0.2)
p_dotplot <- ggplot(low_sim_res_plot, aes(x = mro, y = mip_random, group = group, color=group)) +  
    geom_point(alpha=1, size=1)+
    scale_color_brewer(palette = "Set2") +
    labs(x = "Competition MRO", y = "Cooperation MIP", color="") +  
    theme_bw(base_size=14) +
    scale_x_continuous(
        limits = c(0.5, 1), breaks = seq(0.5, 1, by = 0.1)) + 
    scale_y_continuous(
        limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
    theme(# legend.position = "None",
          text = element_text(face = "bold"))
distribut_mro <- low_sim_res_plot %>%ggplot(aes(x = mro, y = group, color = group)) +
    geom_boxplot(width = 0.3, alpha = 0.5, outlier.shape = NA) +
    scale_color_brewer(palette = "Set2") +
    theme_bw(base_size=14) +
    labs(fill = NULL, x = NULL, y = NULL) +
    scale_x_continuous(
        limits = c(0.4, 1), breaks = seq(0.4, 1, by = 0.1)) + 
    theme(text = element_text(face = "bold"),
          legend.position = "None",
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    annotate("text", x = -Inf, y = Inf, 
             label = ifelse(!is.na(p_mro), sprintf("p=%.2g", p_mro), ""),
             color = "black", hjust = -0.05, vjust = 1.1, size = 4)
distribut_mip <- low_sim_res_plot %>% ggplot(aes(y = mip, x = group, color = group)) +
    geom_boxplot(width = 0.3, alpha = 0.5, outlier.shape = NA) +
    scale_color_brewer(palette = "Set2") +
    theme_bw(base_size=14) +
    labs(fill = NULL, x = NULL, y = NULL) +
    scale_y_continuous(
        limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
    theme(text = element_text(face = "bold"),
          legend.position = "None",
          axis.text = element_blank(),
          axis.ticks = element_blank())+
    annotate("text", x = Inf, y = Inf,  
             label = ifelse(!is.na(p_mip), sprintf("p=%.2g", p_mip), ""),
             color = "black", 
             hjust = -0.09,  
             vjust = 1, 
             angle = 270,  
             size = 4)
p2 <- p_dotplot %>% 
    insert_right(distribut_mip, width = 0.2) %>% 
    insert_top(distribut_mro, height = 0.2)
p2 <- as.ggplot(p2) + 
    labs(title = "cosine ≈ 0") +  
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

p_top <- plot_grid(p1, p2, align="hv",
               nrow = 1, ncol=2, plot=FALSE, rel_widths = c(0.64, 1))

p <- plot_grid(p_all, p_top,
               align="hv", labels = c('a', 'b'),
               nrow = 2, ncol=1, plot=FALSE, rel_heights = c(1, 1))
ggsave("Figures/metabolic_smetana.pdf", p, width = 20, height = 14, units = "cm")

