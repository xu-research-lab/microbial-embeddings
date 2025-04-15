library(tidyverse)
library(cowplot)
library(reshape2)

niche_otu_cos <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/niche/niche_OTU_cos.csv")
p1 <- niche_otu_cos %>% ggplot(aes(x = value)) +
    geom_density(aes(fill = group), alpha = 0.5)+
    # facet_grid(. ~ Type, scales="free") + 
    theme_bw(base_size=14) +
    labs(x = "Cosine", y = "Density", fill = "") +
    scale_color_brewer(palette = "Paired") +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))

studys <- {}
studys[["ibd"]] <- c('PRJNA422193', 'PRJNA431126', 'qiita_1629', 'qiita_2538', 'HMP', 'PRJEB13679', 'PRJNA317429')
studys[["crc"]] <- c('PRJEB36789', 'PRJNA824020', 'PRJDB11845', 'PRJNA318004', 'PRJEB6070', 'PRJNA430990', 'PRJNA290926')

types <- c("ibd", "crc")
p_replace <- list()
i <- 1
title_list <- c("IBD", "CRC")
for (type in types){
    table = read.csv(paste0("/home/dongbiao/word_embedding_microbiome/programe_test/niche/res_replace/res_", type, ".csv"))
    colnames(table) <- c("Original", "Niche", "Phylogeny", "Random")
    
    table <- table %>% melt()
    table$variable <- factor(table$variable, level=c("Original", "Niche", "Phylogeny", "Random"))
    table$study <- rep(studys[[type]], 4)
    
    p_replace[[i]] <- table %>% ggplot(aes(x = variable, y = value, group = variable)) +
        geom_boxplot() + geom_line(aes(group = study, color = study)) +
        geom_point(aes(color = study), size = 3) +
        theme_bw(base_size = 14) +
        labs(x="", y="AUC", title = paste0(title_list[i])) +
        theme(text = element_text(face = "bold"), # 设置所有文本为加粗
              plot.title = element_text(hjust = 0.5), # 居中对齐标题
              axis.title = element_text(face = "bold"), # 加粗坐标轴标题
              axis.text.x = element_text(face = "bold", angle = 15, hjust = 1, vjust = 1)) +
        scale_color_brewer(palette = "Paired")
    
    i = i + 1
}

p <- plot_grid(p1, p_replace[[1]], p_replace[[2]],
                   labels = c('a', 'b', 'c'),
                   scale = c(1, 1, 1),
                   nrow = 1, ncol=3)

ggsave("/home/dongbiao/word_embedding_microbiome/result/Niche_replace.png", p,
       width = 40, height = 12, units = "cm")
