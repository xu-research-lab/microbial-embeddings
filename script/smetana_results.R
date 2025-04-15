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
library(reshape2)
library(gridExtra)
library(patchwork)
library(cowplot)
library(aplot)
library(ggplotify)

library(RColorBrewer)  

library(ggVennDiagram)
# 获取Set1的颜色  
set1_colors <- brewer.pal(n = 9, name = "Set1")  
communities <- c("community_2", "community_4", "community_6", "community_8", "community_10")
plot_distribution <- function(plot_data, legend_bool=FALSE){
    p_dotplot <- ggplot(plot_data, aes(x = mro, y = mip, group = group, color=group)) +  
        geom_point(alpha=0.4)+
        scale_color_manual(values = c("random" = "#E41A1C", "phylogeny" = "#377EB8", 
                                      "cooccurrence"="#4DAF4A", "sparcc" = "#984EA3")) +
        labs(x = "Competition MRO", y = "Cooperation MIP") +  
        theme_bw(base_size=14) +
        theme(legend.position = "None",
              text = element_text(face = "bold"))
    if (legend_bool){
        p_dotplot <- p_dotplot + theme(legend.position = c(0.95, 0.95)) +
            labs(x = "Competition MRO", y = "Cooperation MIP", color = NULL)
    }
    
    mro_p <- plot_data %>% filter(group %in% c("random", "cooccurrence"))
    p_value <- t.test(mro ~ group, data = mro_p)$p.value
    p_value_label <- paste("p = ", format(p_value, digits = 3))
    
    distribut_mro <- plot_data %>% ggplot(aes(x = mro)) +
        geom_density(aes(fill = group), alpha = 0.5)+
        theme_bw(base_size=14) +
        labs(fill = NULL, x = NULL, y = NULL) +
        scale_fill_manual(values = c("random" = "#E41A1C", "phylogeny" = "#377EB8", 
                                     "cooccurrence"="#4DAF4A", "sparcc" = "#984EA3")) +
        theme(text = element_text(face = "bold"), # 设置所有文本为加粗
              legend.position = "None",
              plot.title = element_text(hjust = 0.5), # 居中对齐标题
              axis.title = element_text(face = "bold"), # 加粗坐标轴标题
              axis.text = element_text(face = "bold"), 
              axis.ticks = element_blank(), axis.text.x = element_blank(),
              axis.text.y = element_blank()) 
        # annotate("text", x = Inf, y = Inf, label = p_value_label,  # 在右上角添加 p 值
        #          hjust = 1.1, vjust = 1.5, size = 2.5, color = "#4DAF4A")
        # 
    mip_p <- plot_data %>% filter(group %in% c("random", "cooccurrence"))
    p_value <- t.test(mip ~ group, data = mip_p)$p.value
    p_value_label <- paste("p = ", format(p_value, digits = 3))
    
    distribut_mip <- plot_data %>% ggplot(aes(y = mip)) +
        geom_density(aes(fill = group), alpha = 0.5)+
        scale_fill_brewer(palette = "Set1") +
        theme_bw(base_size=14) +
        labs(fill = NULL, x = NULL, y = NULL) +
        scale_fill_manual(values = c("random" = "#E41A1C", "phylogeny" = "#377EB8", 
                                     "cooccurrence"="#4DAF4A", "sparcc" = "#984EA3")) +
        theme(text = element_text(face = "bold"), # 设置所有文本为加粗
              legend.position = "None",
              plot.title = element_text(hjust = 0.5), # 居中对齐标题
              axis.title = element_text(face = "bold"), # 加粗坐标轴标题
              axis.text = element_text(face = "bold"),axis.ticks = element_blank(), axis.text.x = element_blank(),
              axis.text.y = element_blank()) 
        # annotate("text", x = Inf, y = Inf, label = p_value_label,  # 在右上角添加 p 值
        #           hjust = -0.1, vjust = 1.5, size = 2.5, color = "#4DAF4A", angle = -90)
    
    p <- p_dotplot %>% insert_right(distribut_mip, width=0.2) %>% 
        insert_top(distribut_mro, height = 0.2)
    p <- as.ggplot(p)
}

smetana_res = read.csv("/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_dissimilarity/smetana_res_guild_22.csv") %>% 
    filter(group != "phylogeny")
    filter(group %in% c("random", "phylogeny", "cooccurrence", "sparcc"))
smetana_res = read.csv("/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_dissimilarity/smetana_res_guild_171.csv")
    filter(group %in% c("random", "phylogeny", "cooccurrence", "sparcc"))
plot_list <- list()
for (i in c(1:5)){
    plot_data <- smetana_res %>% filter(community == communities[i])
    if (i == 5){
        legend_bool = TRUE
    } else {
        legend_bool = FALSE
    }
    plot_list[[i]] <- plot_distribution(plot_data, legend_bool)
}
p1 <- plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]],
                labels = c('a:size 2', 'b:size 4', 'c:size 6', 'd:size 8', 'e:size 10'),
                align="hv",nrow = 1, ncol=5, plot=FALSE, rel_widths = c(0.6, 0.6, 0.6, 0.6, 1))


phylogenetic_distance <- read.csv("/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_dissimilarity/phylogenetic_distance_module_171.csv")
    # filter(group %in% c("random", "phylogeny", "cooccurrence"))
p2 <- phylogenetic_distance %>% ggplot(aes(x = phy_distance)) +
    geom_density(aes(fill = group), alpha = 0.5) +
    scale_fill_manual(values = c("random" = "#E41A1C", "sparcc" = "#984EA3",
                                 "phylogeny" = "#377EB8", "module_171"="#4DAF4A")) +
    theme_bw(base_size=14) +
    labs(fill = NULL, x = "Phylogenetic distance", y = "Density") +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          legend.position = "None",
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))

genes_table <- read.csv("/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_dissimilarity/metabolic_genes.csv") %>% 
    filter(group %in% c("random", "phylogeny", "cooccurrence"))
p3  <- genes_table %>% ggplot(aes(x = genes_num)) +
    geom_density(aes(fill = group), alpha = 0.5)+
    theme_bw(base_size=14) +
    scale_fill_manual(values = c("random" = "#E41A1C", 
                                 "phylogeny" = "#377EB8", "cooccurrence"="#4DAF4A")) +
    theme_bw(base_size=14) +
    labs(fill = NULL, x = "No. of metabolic genes", y = "Density") +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          legend.position = "None",
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))

metabolic_dissimilarity <- read.csv("/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_dissimilarity/metabolic_dissimilarity.csv")
p4 <- metabolic_dissimilarity  %>% ggplot(aes(x = jaccard_distance)) +
    geom_density(aes(fill = group), alpha = 0.5)+
    theme_bw(base_size=14) +
    scale_fill_manual(values = c("random" = "#E41A1C", 
                                 "phylogeny" = "#377EB8", "cooccurrence"="#4DAF4A")) +
    theme_bw(base_size=14) +
    labs(fill = NULL, x = "Metabolic dissimilarity", y = "Density") +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          legend.position = "None",
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))

### cooccurrence nums otu
nums_co_otus <- read.csv("/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_dissimilarity/nums_co_otus.csv")
nums_co_otus_cooccurrence <- nums_co_otus %>% filter(group == "cooccurrence") %>% 
    arrange(desc(values))
nums_co_otus_cooccurrence$module_id <- factor(nums_co_otus_cooccurrence$module_id, levels = nums_co_otus_cooccurrence$module_id)
p5_1 <- nums_co_otus_cooccurrence %>% ggplot(aes(x = module_id, y = log10(values+1))) +
    labs(x="Co-occurrence module_id", y="log10(num)") +
    geom_col(fill = "skyblue", alpha=0.5) +
    theme_bw(base_size=14) +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          legend.position = "None",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))

nums_co_otus_phylogeny <- nums_co_otus %>% filter(group == "phylogeny") %>% 
    arrange(desc(values))
nums_co_otus_phylogeny$module_id <- factor(nums_co_otus_phylogeny$module_id, levels = nums_co_otus_phylogeny$module_id)
p5_2 <- nums_co_otus_phylogeny %>% ggplot(aes(x = module_id, y = log10(values+1))) +
    labs(x="phylogeny module_id", y="log10(num)") +
    geom_col(fill = "skyblue", alpha=0.5) +
    theme_bw(base_size=14) +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          legend.position = "None",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))
p5 <- plot_grid(p5_1, p5_2, 
                nrow = 2, ncol=1, plot=FALSE)

### module 93 cooccurrence
co_module_93_mip_mro <- read.csv("/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_dissimilarity/co_module_93_mip_mro.csv")

p6 <- plot_distribution(co_module_93_mip_mro)

### module 93 produce metabolic
chi2_results <- read.csv("/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_dissimilarity/metabolic_produce_module_111.csv")
chi2_results$q_value <- p.adjust(chi2_results$p_value_res)
chi2_results <- chi2_results %>% filter(q_value < 0.05) %>% filter(rate > 1)
chi2_results <- chi2_results %>% arrange(desc(rates_module))
colnames(chi2_results) <- c("meta_id", "chi2_res", "p_value_res", "Backgroup", "Module_93", "q_value")
fact_leves <- chi2_results$meta_id
chi2_results_plot <- chi2_results[, c("meta_id", "Backgroup", "Module_93")] %>% melt(id.vars=c("meta_id"))
chi2_results_plot$meta_id <- factor(chi2_results_plot$meta_id, levels= fact_leves)
p7 <- chi2_results_plot %>% ggplot(aes(x = meta_id, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Metabolites", y = "Probability", fill = "") +
    theme_bw(base_size = 14) + 
    scale_fill_brewer(palette = "Set1") +
    # scale_x_discrete(labels = c("fe2_c"="Fe2+", "chols_c"="Choline sulfate", "bz_c"="Benzoate", 
    #                             "galman6_c"="Galactomannan", "dcyt_c"="Deoxycytidine", "gal_c"="D-Galactose",
    #                             "6pgc_c"="6-Phospho-D-gluconate", "ppt_c"="Phosphonate", "co_c"="Carbon monoxide"))+
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"),
          legend.position = c(0.9, 0.9),  # 设置图例位置，(0.9, 0.9)表示右上角  
          legend.justification = c("right", "top"),
          axis.text.x = element_text(angle = -10, hjust = 0, vjust = 1))
    
### cooccurrence with module 93 microbiome cooperative communities
cross_feed_res <- read.csv("/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_dissimilarity/cross_feed_res.csv")
p8 <- cross_feed_res %>% ggplot(aes(x=compound, y=mus, fill=group))+
    geom_boxplot(alpha=0.4) +
    scale_fill_manual(values = c("random" = "#E41A1C", 
                                 "phylogeny" = "#377EB8", "cooccurrence"="#4DAF4A")) +
    labs(x = "", y = "Mus", fill = "") +
    theme_bw(base_size = 14) + 
    scale_x_discrete(labels = c("M_fe3pyovd_kt_e"="Ferrypyoverdine", 
                                "M_mal__L_e"="L-Malate", 
                                "M_nh4_e"="Ammonium"))+
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = -10, hjust = 0, vjust = 1))

middle_plot <- plot_grid(p2, p4, p5_1, p5_2, 
                           labels = c('f', 'g', 'h', 'i'),
                           nrow = 1, ncol=4, plot=FALSE, rel_widths = c(0.7, 0.7, 1, 1))

buttom_plot <- plot_grid(p6, p7, p8, 
                           labels = c('k', 'l', 'm'),
                           nrow = 1, ncol=3, plot=FALSE)

p <- plot_grid(p1, middle_plot,buttom_plot,
               align="hv",
               nrow = 3, ncol=1, plot=FALSE, rel_heights = c(1, 1, 1))

ggsave("/home/dongbiao/word_embedding_microbiome/result/smetana.png", p,
       width = 38, height = 20, units = "cm")




### Module 33 relative with disease
smetana_res = read.csv("/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_dissimilarity/smetana_res_module_33.csv")
plot_list <- list()
for (i in c(1:5)){
    plot_data <- smetana_res %>% filter(community == communities[i])
    if (i == 5){
        legend_bool = TRUE
    } else {
        legend_bool = FALSE
    }
    plot_list[[i]] <- plot_distribution(plot_data, legend_bool)
}
p1 <- plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]],
                labels = c('a:size 2', 'b:size 4', 'c:size 6', 'd:size 8', 'e:size 10'),
                align="hv",nrow = 1, ncol=5, plot=FALSE, rel_widths = c(0.6, 0.6, 0.6, 0.6, 1))

phylogenetic_distance <- read.csv("/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_dissimilarity/phylogenetic_distance_module_33.csv") %>%
    mutate(group = ifelse(group == "module_33", "cooccurrence", group))
p2 <- phylogenetic_distance %>% ggplot(aes(x = phy_distance)) +
    geom_density(aes(fill = group), alpha = 0.5)+
    scale_fill_manual(values = c("random" = "#E41A1C", 
                                 "phylogeny" = "#377EB8", "cooccurrence"="#4DAF4A")) +
    theme_bw(base_size=14) +
    labs(fill = NULL, x = "Phylogenetic distance", y = "Density") +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          legend.position = "None",
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))


metabolic_dissimilarity <- read.csv("/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_dissimilarity/metabolic_dissimilarity_module_33.csv")
p3 <- metabolic_dissimilarity  %>% ggplot(aes(x = jaccard_distance)) +
    geom_density(aes(fill = group), alpha = 0.5)+
    theme_bw(base_size=14) +
    scale_fill_manual(values = c("random" = "#E41A1C", 
                                 "phylogeny" = "#377EB8", "cooccurrence"="#4DAF4A")) +
    theme_bw(base_size=14) +
    labs(fill = NULL, x = "Metabolic dissimilarity", y = "Density") +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          legend.position = "None",
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))

### co-occurrence
co_module_33_mip_mro <- read.csv("/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_dissimilarity/co_module_33_mip_mro.csv")

p4 <- plot_distribution(co_module_33_mip_mro)

### module 33 produce metabolic
chi2_results <- read.csv("/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_dissimilarity/metabolic_produce_module_33.csv")
chi2_results$q_value <- p.adjust(chi2_results$p_value_res)
chi2_results <- chi2_results %>% filter(q_value < 0.05) %>% filter(rate > 1)
chi2_results <- chi2_results %>% arrange(desc(rates_module))
chi2_results$meta_id <- c("Uracil", "Phosphate", "Trimethylamine", "Formaldehyde",                                    
                          "(2R,4S)-2-methyl-2,3,3,4-\ntetrahydroxytetrahydrofuran",
                          "L-Threonine", "Salmochelin-S2")
colnames(chi2_results) <- c("meta_id", "chi2_res", "p_value_res", "Backgroup", "Module_33", "q_value")
fact_leves <- chi2_results$meta_id
chi2_results_plot <- chi2_results[, c("meta_id", "Backgroup", "Module_33")] %>% melt(id.vars=c("meta_id"))
chi2_results_plot$meta_id <- factor(chi2_results_plot$meta_id, levels= fact_leves)
p5 <- chi2_results_plot %>% ggplot(aes(x = meta_id, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Metabolites", y = "Probability", fill = "") +
    theme_bw(base_size = 10) + 
    scale_fill_brewer(palette = "Set1") +
    # scale_x_discrete(labels = c("fe2_c"="Fe2+", "chols_c"="Choline sulfate", "bz_c"="Benzoate", 
    #                             "galman6_c"="Galactomannan", "dcyt_c"="Deoxycytidine", "gal_c"="D-Galactose",
    #                             "6pgc_c"="6-Phospho-D-gluconate", "ppt_c"="Phosphonate", "co_c"="Carbon monoxide"))+
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"),
          legend.position = c(0.9, 0.9),  # 设置图例位置，(0.9, 0.9)表示右上角  
          legend.justification = c("right", "top"),
          axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))

### cooccurrence with module 93 microbiome cooperative communities
cross_feed_res <- read.csv("/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_dissimilarity/cross_feed_res_module_33.csv")
p6 <- cross_feed_res %>% ggplot(aes(x=compound, y=mus, fill=group))+
    geom_boxplot(alpha=0.4) +
    scale_fill_manual(values = c("random" = "#E41A1C", 
                                 "phylogeny" = "#377EB8", "cooccurrence"="#4DAF4A")) +
    labs(x = "", y = "Mus", fill = "") +
    theme_bw(base_size = 14) + 
    scale_x_discrete(labels = c("M_ins_e"="Inosine"))+
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 0, hjust = 0, vjust = 1))

buttom_plot <- plot_grid(p2, p3, p4, p5, p6, 
                         labels = c('f', 'g', 'h', 'i', 'j'),
                         nrow = 1, ncol=5, plot=FALSE)

p <- plot_grid(p1,buttom_plot,
               align="hv",
               nrow = 2, ncol=1, plot=FALSE, rel_heights = c(1, 1))

ggsave("/home/dongbiao/word_embedding_microbiome/result/smetana_module_33.png", p,
       width = 38, height = 14, units = "cm")


### module 93 and module 94
set1_colors <- brewer.pal(n = 9, name = "Set1")  
communities <- c("community_2", "community_4", "community_6", "community_8", "community_10")
plot_distribution <- function(plot_data, legend_bool=FALSE){
    p_dotplot <- ggplot(plot_data, aes(x = mro, y = mip, group = group, color=group)) +  
        geom_point(alpha=0.4)+
        scale_color_manual(values = c("random" = "#E41A1C", "module 94"="#984EA3",
                                      "phylogeny" = "#377EB8", "module 93"="#4DAF4A", "module 33"="#FF7F00")) +
        labs(x = "Competition MRO", y = "Cooperation MIP") +  
        theme_bw(base_size=14) +
        theme(legend.position = "None",
              text = element_text(face = "bold"))
    if (legend_bool){
        p_dotplot <- p_dotplot + theme(legend.position = c(0.95, 0.95)) +
            labs(x = "Competition MRO", y = "Cooperation MIP", color = NULL)
    }
    
    distribut_mro <- plot_data %>% ggplot(aes(x = mro)) +
        geom_density(aes(fill = group), alpha = 0.5)+
        theme_bw(base_size=14) +
        labs(fill = NULL, x = NULL, y = NULL) +
        scale_fill_manual(values = c("random" = "#E41A1C", "module 94"="#984EA3",
                                     "phylogeny" = "#377EB8", "module 93"="#4DAF4A", "module 33"="#FF7F00")) +
        theme(text = element_text(face = "bold"), # 设置所有文本为加粗
              legend.position = "None",
              plot.title = element_text(hjust = 0.5), # 居中对齐标题
              axis.title = element_text(face = "bold"), # 加粗坐标轴标题
              axis.text = element_text(face = "bold"), 
              axis.ticks = element_blank(), axis.text.x = element_blank(),
              axis.text.y = element_blank()) 
    
    distribut_mip <- plot_data %>% ggplot(aes(y = mip)) +
        geom_density(aes(fill = group), alpha = 0.5)+
        scale_fill_brewer(palette = "Set1") +
        theme_bw(base_size=14) +
        labs(fill = NULL, x = NULL, y = NULL) +
        scale_fill_manual(values = c("random" = "#E41A1C", "module 94"="#984EA3",
                                     "phylogeny" = "#377EB8", "module 93"="#4DAF4A", "module 33"="#FF7F00")) +
        theme(text = element_text(face = "bold"), # 设置所有文本为加粗
              legend.position = "None",
              plot.title = element_text(hjust = 0.5), # 居中对齐标题
              axis.title = element_text(face = "bold"), # 加粗坐标轴标题
              axis.text = element_text(face = "bold"),axis.ticks = element_blank(), axis.text.x = element_blank(),
              axis.text.y = element_blank())
    
    p <- p_dotplot %>% insert_right(distribut_mip, width=0.2) %>% 
        insert_top(distribut_mro, height = 0.2)
    p <- as.ggplot(p)
}

smetana_res = read.csv("/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_dissimilarity/smetana_res_module93_module94_module33.csv")
plot_list <- list()
for (i in c(1:5)){
    plot_data <- smetana_res %>% filter(community == communities[i])
    if (i == 5){
        legend_bool = TRUE
    } else {
        legend_bool = FALSE
    }
    plot_list[[i]] <- plot_distribution(plot_data, legend_bool)
}
p1 <- plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]],
                labels = c('a:size 2', 'b:size 4', 'c:size 6', 'd:size 8', 'e:size 10'),
                align="hv",nrow = 1, ncol=5, plot=FALSE, rel_widths = c(0.6, 0.6, 0.6, 0.6, 1))

### co-occurrence
co_module_93_94_mip_mro <- read.csv("/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_dissimilarity/co_module_33_93_94_mip_mro.csv")

p2 <- plot_distribution(co_module_93_94_mip_mro)

### veen
M_93 <- co_module_93_94_mip_mro %>% filter(group == "module 93")
M_93 <- as.character(M_93$id_2)
M_94 <- co_module_93_94_mip_mro %>% filter(group == "module 94")
M_94 <- as.character(M_94$id_2)
M_33 <- co_module_93_94_mip_mro %>% filter(group == "module 33")
M_33 <- as.character(M_33$id_2)

# 将数据放入一个列表
sets <- list(A = M_93, B = M_94, C = M_33)
p3 <- ggVennDiagram(sets, category.names = c("module 93","module 94", "module 33"),
                    label = "count",
                    label_color = "black",
                    label_alpha = 0,
                    set_size = 4,
                    edge_lty = "dashed", 
                    edge_size = 1) +
    scale_fill_gradient(low="white",high = "#b9292b",name = "OTU count") +
    scale_x_continuous(expand = expansion(mult = .2))

chi2_results <- read.csv("/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_dissimilarity/metabolic_produce_module_94.csv")
chi2_results$q_value <- p.adjust(chi2_results$p_value_res)
chi2_results <- chi2_results %>% filter(q_value < 0.05) %>% filter(rate > 1)
chi2_results <- chi2_results %>% arrange(desc(rates_module))
colnames(chi2_results) <- c("meta_id", "chi2_res", "p_value_res", "Backgroup", "Module_93", "q_value")
fact_leves <- chi2_results$meta_id
chi2_results_plot <- chi2_results[, c("meta_id", "Backgroup", "Module_93")] %>% melt(id.vars=c("meta_id"))
chi2_results_plot$meta_id <- factor(chi2_results_plot$meta_id, levels= fact_leves)
p4 <- chi2_results_plot %>% ggplot(aes(x = meta_id, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Metabolites", y = "Probability", fill = "") +
    theme_bw(base_size = 14) + 
    scale_fill_brewer(palette = "Set1") +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"),
          legend.position = c(0.9, 0.9),  # 设置图例位置，(0.9, 0.9)表示右上角  
          legend.justification = c("right", "top"),
          axis.text.x = element_text(angle = -10, hjust = 0, vjust = 1))

buttom_plot <- plot_grid(p2, p3, NULL, NULL, NULL, 
                         labels = c('f', 'g'),
                         nrow = 1, ncol=5, plot=FALSE)

p <- plot_grid(p1,buttom_plot,
               align="hv",
               nrow = 2, ncol=1, plot=FALSE, rel_heights = c(1, 1))

ggsave("/home/dongbiao/word_embedding_microbiome/result/smetana_module_93_94.png", p,
       width = 38, height = 14, units = "cm")
