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

### SpectralClustering 算法评估
eigenvalues <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/cluster/eigenvalues.csv")
eigenvalues$x <- eigenvalues$x + 1
diff_eigenvalues <- diff(eigenvalues$y)  
topK <- 5  # Set your desired number of top gaps  
index_largest_gap <- order(diff_eigenvalues, decreasing = TRUE)[1:topK] 

p1 <- eigenvalues %>% ggplot(aes(x=x, y=y)) +
    geom_point(color="skyblue") +
    theme_bw() +
    labs(x="Num modules", y="Eigenvalues", title=paste0("Optimal number of modules: \n", 
                                                         index_largest_gap[1], " ", index_largest_gap[2], " ", index_largest_gap[3]))

### clusters distribution
spectral_clustering <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/cluster/spectral_clustering.csv")
co_clusters_cosine <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/cluster/co_clusters_cosine.csv")
num_modules <- 160
co_clusters_cosine <- co_clusters_cosine %>% filter(clusters_mean > 0.6)
keep_num_modules <- nrow(co_clusters_cosine)
p2 <- co_clusters_cosine %>% ggplot(aes(x = clusters_mean)) +
    geom_density(fill = "skyblue", alpha = 1) +
    theme_bw() +
    labs(x="Cosine similarity of insider module OTUs", 
         y="Density", title=paste0("Num modules: ", keep_num_modules))

spectral_clustering <- spectral_clustering %>% filter(clusters %in% co_clusters_cosine$clusters)
element_counts <- table(spectral_clustering$clusters)
clusters_num <- data.frame(clusters_id=names(element_counts), num_otus=as.vector(element_counts))
p3 <- clusters_num %>% ggplot(aes(x = num_otus)) +
    geom_density(fill = "skyblue", alpha = 1) +
    theme_bw() +
    scale_x_continuous(breaks=c(25, 50, 75, 100, 125))+
    labs(x="Num OTU of insider module", 
         y="Density", title=paste0("Num modules: ", keep_num_modules))

### module pathway sim
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
p4 <- p_kmeans_func_embed %>% insert_right(distribut_func, width=0.2) %>% 
    insert_top(distribut_embed, height = 0.2)
p4 <- as.ggplot(p4)

### Smetana
k_means_pathway_embed <- read.csv("/home/dongbiao/word_embedding_microbiome/modelseed/results/smetanta_plot.csv")
p_kmeans_func_embed <- ggplot(k_means_pathway_embed, aes(x = embed_sim, y = semtana, group = group, fill=group)) +  
    geom_hdr(show.legend = FALSE, method = "kde", n=100,
             probs = c(0.85, 0.8, 0.7, 0.6, 0.5, 0.4, 0.2, 0.25, 0.15, 0.1, 0.05)) +
    scale_fill_manual(values = c("sim" = "blue", "unsim" = "red")) +  
    labs(x = "Embedding Similarity", y = "Smetana") +  
    theme_bw(base_size=14) +
    theme(legend.position = "right", text = element_text(face = "bold"))
# Embedding and Function Similarity by Jaccard Similarity Level

distribut_embed <- k_means_pathway_embed %>% ggplot(aes(x = embed_sim)) +
    geom_density(aes(fill = group), alpha = 0.4)+
    theme_bw(base_size=14) +
    labs(fill = NULL, x = NULL, y = NULL) +
    scale_fill_manual(values = c("sim" = "blue", "unsim" = "red")) + 
    scale_color_brewer(palette = "Paired") +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          legend.position = "None",
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"), 
          axis.ticks = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_blank())


distribut_func <- k_means_pathway_embed %>% ggplot(aes(y = semtana)) +
    geom_density(aes(fill = group), alpha = 0.4, kernel="gaussian", n=512)+
    scale_fill_manual(values = c("sim" = "blue", "unsim" = "red")) +
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

upp_plot <- plot_grid(p1, p2, p3,
                      labels = c('a', 'b', 'c'),
                      align="hv",
                      nrow = 1, ncol=3, plot=FALSE, rel_widths = c(1, 1, 1))

middle_plot <- plot_grid(p4, NULL,
                         labels = c('d', ''),
                         align="hv",
                         nrow = 1, ncol=2, plot=FALSE, rel_widths = c(1, 1))

p <- plot_grid(upp_plot, middle_plot,
               align="hv",
               nrow = 2, ncol=1, plot=FALSE, rel_heights = c(0.5, 1))
ggsave("/home/dongbiao/word_embedding_microbiome/result/embed_cluster.png", p,
       width = 25, height = 20, units = "cm")

