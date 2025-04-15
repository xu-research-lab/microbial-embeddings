library(dplyr)
library(ggplot2)
library(lemon)
library(ggh4x)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(grid)
library(stringr)
library(patchwork)
library(reshape2)

# function
# read in results
get_auc_f1 <- function(file, idx){
    table <- read.table(file, sep = "\t")
    table <- table[apply(table, 1, function(x) grepl("mcc", x)), ]
    table <- gsub(",", "", table)
    table <- str_split(table, " ", simplify = TRUE)
    table <- as.data.frame(table[, idx])
}

# plot for benchmark results function
lou_plot <- function(plot_data, AUC_or_F1, vector_x, nums, titles){
    colnames(plot_data) <- c("value", "v_mean", "v_sd", "Group1", "Group2", "Study")
    plot_data$Group1 <- factor(plot_data$Group1, levels = c("Abundance-percentile", "Russell_rao_weight", "OTU"))
    p <- plot_data %>% ggplot(aes(x = Group1, y = v_mean / nums)) +
        labs(x="", y=paste0(AUC_or_F1), title=titles) +
        geom_col(fill = "skyblue", alpha=0.5) +
        geom_jitter(aes(x = Group1, y = value, color=Study), width = 0.1, 
                    height = 0, size = 2, alpha = 0.8) +
        facet_wrap(.~Group2, scales = "free_x")+
        geom_errorbar(aes(ymin = v_mean - v_sd, ymax = v_mean + v_sd), width = 0.2, color = "black", size = 0.7) +
        theme_bw(base_size = 14) +
        scale_color_brewer(palette = "Paired") +
        theme(text = element_text(face = "bold"), # 设置所有文本为加粗
              plot.title = element_text(hjust = 0.5), # 居中对齐标题
              axis.title = element_text(face = "bold"), # 加粗坐标轴标题
              axis.text = element_text(face = "bold")) + # 加粗坐标轴文本
        coord_cartesian(ylim = vector_x) + theme(legend.position = 'none') +
        theme(axis.text.x = element_text(angle = -10, hjust = 0, vjust = 1))
    
    return(p)
}

plot_banchmark <- function(lable_1, attention_res, mlp_res, rf_res, 
                           num_study, study_name, vector_x_1, vector_x_2, titles){
    plot_list <- list()
    table_1 <- get_auc_f1(attention_res, c(12, 18))
    table_2 <- get_auc_f1(mlp_res, c(12, 18))
    table_3 <- get_auc_f1(rf_res, c(4, 12))
    table <- rbind(table_1, table_2, table_3)
    table <- as.data.frame(lapply(table, as.numeric))
    colnames(table) <- c("AUC", "F1")
    table$Group1 <- c(rep(lable_1, 1, each=num_study),
                      rep(lable_1, 1, each=num_study),
                      rep(lable_1, 1, each=num_study))
    table$Group2 <- rep(c("Attention", "MLP", "RF"), 1, each= num_study * 11)
    table$Study <- rep(study_name, 33, each=1)
    table <- table %>% filter(Group1 %in% c("Russell_rao_weight", "Abundance-percentile", "OTU"))
    table_mean_sd <- table %>% group_by(Group1, Group2) %>% summarise(AUC_mean = mean(AUC), AUC_sd = sd(AUC),
                                                                      F1_mean = mean(F1), F1_sd = sd(F1))
    table <- merge(table, table_mean_sd, by = c("Group1", "Group2"))
    group <- c("Attention", "MLP", "RF")
    
    plot_list[[1]] <- lou_plot(plot_data = table[,c("AUC", "AUC_mean", "AUC_sd", "Group1", "Group2", "Study")], 
                               AUC_or_F1 = "AUC", 
                               vector_x = vector_x_1, nums = num_study, titles)
    
    plot_list[[2]] <- lou_plot(plot_data = table[,c("F1", "F1_mean", "F1_sd", "Group1", "Group2", "Study")], 
                               AUC_or_F1 = "F1", vector_x = vector_x_2, nums = num_study, titles)
    
    return(plot_list)
}

# plot for diff dataset influence model preformance function
plot_diff_dataset <- function(cross_res, study_name, num_study, titles){
    cross_diff_size <- get_auc_f1(cross_res, c(12, 18))
    colnames(cross_diff_size) <- c("AUC", "F1")
    cross_diff_size <- as.data.frame(lapply(cross_diff_size, as.numeric))
    cross_diff_size$Study <- rep(study_name, 5, each=10)
    cross_diff_size$datasize <- rep(c("5000", "10000", "20000", "40000", "80000"), 1, 
                                    each=10 * num_study)
    line_plot <- cross_diff_size %>% group_by(Study, datasize) %>% 
        summarise(AUC_mean = mean(AUC), AUC_sd = sd(AUC),
                  F1_mean = mean(F1), F1_sd = sd(F1))
    cross_diff_size <- merge(cross_diff_size, line_plot, by=c("Study", "datasize"))
    cross_diff_size$datasize <- factor(cross_diff_size$datasize, 
                                       levels=c("5000", "10000", "20000", "40000", "80000"))
    
    plot_list <- list()
    cross_diff_size$datasize <- as.numeric(cross_diff_size$datasize)
    plot_list[[1]] <- cross_diff_size %>% ggplot(aes(x = datasize, y = AUC)) + 
        geom_point(aes(color = Study), size = 2, alpha=0.8) + 
        geom_errorbar(aes(ymin = AUC_mean - AUC_sd, ymax = AUC_mean + AUC_sd), width = 0.1) + 
        geom_line(aes(group = Study, y = AUC_mean, color = Study), alpha = 0.5, linewidth = 1, linetype = "dashed") +
        
        labs(x="Datasize", y="AUC", title=titles) +
        geom_smooth(level= 0.9) +
        theme_bw(base_size = 14) +
        scale_color_brewer(palette = "Paired") +
        scale_x_continuous(breaks = c(1, 2, 3, 4, 5),  # 指定刻度位置  
                           labels = c("5000", "10000", "20000", "40000", "80000"))+
        theme(text = element_text(face = "bold"), # 设置所有文本为加粗
              plot.title = element_text(hjust = 0.5), # 居中对齐标题
              axis.title = element_text(face = "bold"), # 加粗坐标轴标题
              axis.text = element_text(face = "bold"))

    plot_list[[2]] <- cross_diff_size %>% ggplot(aes(x = datasize, y = F1)) + 
        geom_boxplot() + geom_point(aes(color = Study), size = 2, alpha=0.8) +
        geom_line(aes(group = Study, y = F1_mean, color = Study), alpha = 0.5, linewidth = 1) +
        geom_errorbar(aes(ymin = F1_mean - F1_sd, ymax = F1_mean + F1_sd), width = 0.1) + 
        labs(x="Datasize", y="F1", title=titles) +
        theme_bw(base_size = 14) +
        scale_color_brewer(palette = "Paired") +
        theme(text = element_text(face = "bold"), # 设置所有文本为加粗
              plot.title = element_text(hjust = 0.5), # 居中对齐标题
              axis.title = element_text(face = "bold"), # 加粗坐标轴标题
              axis.text = element_text(face = "bold"))
    
    return(plot_list)
}

# Banchmark
lable_1 <- c("Russell_rao_weight", "Russell_rao", "faith", "Jaccard", "Abundance-percentile", 
             "Abundance-totalsum", "Braycurtis-percentile", "Braycurtis-totalsum", 
             "Phylogeny", "PCA", "OTU")
# IBD
attention_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/filter_prevalence/all_feces/IBD_cross/attention_2.out"
mlp_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/IBD/mlp_all_feces.out"
rf_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/IBD/randownforest_all_feces.out"
IBD_p <- plot_banchmark(lable_1, attention_res, mlp_res, rf_res, 
                        num_study=7, 
                        study_name=c('PRJNA422193', 'PRJNA431126', 'qiita_1629', 'qiita_2538', 'HMP', 'PRJEB13679', 'PRJNA317429'),
                        vector_x_1=c(0.4, 1), vector_x_2=c(0.18, 0.84), titles="IBD")


# CRC
attention_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/filter_prevalence/all_feces/CRC_cross/attention_2.out"
mlp_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/CRC/mlp_all_feces.out"
rf_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/CRC/randownforest_all_feces.out"
CRC_p <- plot_banchmark(lable_1, attention_res, mlp_res, rf_res, 
                        num_study=7, 
                        study_name=c('PRJEB36789', 'PRJNA824020', 'PRJDB11845', 'PRJNA318004',
                                     'PRJEB6070', 'PRJNA430990', 'PRJNA290926'),
                        vector_x_1=c(0.32, 1), vector_x_2=c(0.14, 0.88), titles="CRC")

# fiber
# attention_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/filter_prevalence/all_feces/dietary_cross/attention_1.out"
# mlp_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/mlp_all_feces.out"
# rf_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/randownforest_all_feces.out"
# fiber_p <- plot_banchmark(lable_1, lable_2, attention_res, mlp_res, rf_res, 
#                         num_study=6, 
#                         study_name=c('qiita_13367', 'SRP067761', 'SRP120250', 'SRP128128', 'SRP345891', 'SRP219296'),
#                         vector_x_1=c(0.30, 0.9), vector_x_2=c(0.16, 0.85))
# diff data size
# IBD
cross_IBD <- plot_diff_dataset(cross_res="/home/dongbiao/word_embedding_microbiome/programe_test/filter_prevalence/sub_embedding/IBD_cross/attention.out", 
                  study_name=c('PRJNA422193', 'PRJNA431126', 'qiita_1629', 'qiita_2538', 'HMP', 'PRJEB13679', 'PRJNA317429'),
                  num_study=7, titles="IBD")

# CRC
cross_CRC <- plot_diff_dataset(cross_res="/home/dongbiao/word_embedding_microbiome/programe_test/filter_prevalence/sub_embedding/CRC_cross/attention.out", 
                               study_name=c('PRJEB36789', 'PRJNA824020', 'PRJDB11845', 'PRJNA318004',
                                            'PRJEB6070', 'PRJNA430990', 'PRJNA290926'),
                               num_study=7, titles="CRC")

### embedding remove batch effects

# load("/home/dongbiao/word_embedding_microbiome/programe_test/IBD/explain/adonise2_res.RData")
# adonis_ibd <- res
# adonis_ibd$group2 <- rep("IBD",nrow(adonis_ibd))
# load("/home/dongbiao/word_embedding_microbiome/programe_test/CRC/explain_model/adonise2_res.RData")
# adonis_crc <- res
# adonis_crc$group2 <- rep("CRC",nrow(adonis_crc))
# adonis_res <- rbind(adonis_ibd, adonis_crc)
# 
# adonis_res <- adonis_res %>% mutate(pvalue = ifelse(P < 0.05, "P < 0.05", "ns"))
# adonis_res$group2 <- factor(adonis_res$group2, levels = c("IBD", "CRC"))
# p1 <- ggplot(adonis_res, aes(x = R2, y = study, color = group, shape = pvalue)) +  
#     geom_point(size = 3) +  
#     scale_shape_manual(values = c(16, 17)) + # 圆点: ns，三角: P < 0.05  
#     theme_bw(base_size=14) +  
#     facet_wrap(.~group2, scales = "free") +
#     labs(x = "R2", y = "", color = "", shape = "Pvalue") +  
#     scale_color_brewer(palette = "Paired") +
#     theme(axis.text.y = element_text(size = 10),  
#           axis.text.x = element_text(size = 10),  
#           legend.position = "right") +
#     theme(text = element_text(face = "bold"), # 设置所有文本为加粗
#           plot.title = element_text(hjust = 0.5), # 居中对齐标题
#           axis.title = element_text(face = "bold"), # 加粗坐标轴标题
#           axis.text = element_text(face = "bold"))
adonis_ibd <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/IBD/explain/adonise2_res.csv")
adonis_crc <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/dim_reduction/CRC/adonise2_res.csv")
adonis_res <- rbind(adonis_ibd, adonis_crc)
adonis_res$group2 <- factor(adonis_res$group2, levels = c("Group", "Region", "Country", "Study"))
adonis_res$type <- factor(adonis_res$type, levels = c("IBD", "CRC"))
adonis_res %>% ggplot(aes(x=group2, y=R2, fill=group1)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(.~type, scales = "free") +
    scale_fill_brewer(palette = "Paired") +
    theme(axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          legend.position = "right") +
    labs(x = "R2", y = "", fill = "") +  
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))


### Attention weight
attention_weight_ibd <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/IBD/explain/attention_weight.csv")
attention_weight_crc <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/CRC/explain_model/attention_weight.csv")
attention_weight_ibd$group3 <- rep("IBD", n=nrow(attention_weight_ibd))
attention_weight_crc$group3 <- rep("CRC", n=nrow(attention_weight_ibd))
# attention_weight_crc$group3 <- c(rep("The ratio of biomark \n attention weight in CRC", 7), "CRC")
# attention_weight_ibd$group <- c("Model1", "Model2", "Model3", "Model4", "Model5", "Model6", "Model7", "Biomark / all taxa")
# attention_weight_crc$group <- c("Model1", "Model2", "Model3", "Model4", "Model5", "Model6", "Model7", "Biomark / all taxa")
attention_weight <- rbind(attention_weight_ibd, attention_weight_crc)
attention_weight$group3 <- factor(attention_weight$group3, c("IBD", "CRC"))
attention_weight <- attention_weight %>% group_by(group1, group2, group3) %>% summarise(attention_weight=mean(attention_weight))

p2 <- attention_weight %>% ggplot(aes(x=group2, y=attention_weight, fill=group1)) +
    geom_bar(stat='identity', position = 'dodge') + 
    facet_grid(~group3, scales = "free", space = "free_x") +
    theme_bw(base_size = 14)+
    theme(text = element_text(face = "bold"), 
          axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),
          legend.text = element_text(size = 14, color = "black", face = "bold"),
          legend.box.margin = margin(l = 0.1, unit = "in"),
          # legend.position = "top",
          legend.title.align = 0.5 
    ) +
    labs(x="", y="Average of attention weight", fill="")

niche_otu_cos <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/niche/niche_OTU_cos.csv")

p3 <- niche_otu_cos %>% ggplot(aes(x = value)) +
    geom_density(aes(fill = group), alpha = 0.5)+
    # facet_grid(. ~ Type, scales="free") + 
    theme_bw(base_size=14) +
    labs(x = "Cosine", y = "Density", fill = "") +
    scale_color_brewer(palette = "Paired") +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          legend.position = "top",
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))

studys <- {}
studys[["ibd"]] <- c('PRJNA422193', 'PRJNA431126', 'qiita_1629', 'qiita_2538', 'HMP', 'PRJEB13679', 'PRJNA317429')
studys[["crc"]] <- c('PRJEB36789', 'PRJNA824020', 'PRJDB11845', 'PRJNA318004', 'PRJEB6070', 'PRJNA430990', 'PRJNA290926')
types <- c("ibd", "crc")
p_replace <- list()
get_numbirc <- function(vectors){
    vectors <- gsub("\\[|\\]", "", vectors)
    vectors <- as.numeric(unlist(strsplit(gsub(" ", "", vectors), ",")))
    return(vectors)
}
i <- 1
title_list <- c("IBD", "CRC")
for (type in types){
    table = read.csv(paste0("/home/dongbiao/word_embedding_microbiome/programe_test/niche/res_replace/res_", type, ".csv")) %>% 
        t() %>% as.data.frame()
    Niche <- c()
    Phylogeny <- c()
    Random <- c()
    for (n in c(1:7)){
        Niche <- c(Niche, get_numbirc(table[n, 2]))
        Phylogeny <- c(Phylogeny, get_numbirc(table[n, 3]))
        Random <- c(Random, get_numbirc(table[n, 4]))
    }
    value <- c(table[,1], Niche, Phylogeny, Random)
    value <- as.numeric(value)
    Group1 <- c(rep("Original", 7), rep("Niche", 63), rep("Phylogeny", 63), rep("Random", 63))
    Group2 <- c(studys[[type]], rep(rep(studys[[type]], n=7, each=9), 3))
    table <- data.frame(value = value, Group1=Group1, Group2=Group2)
    table <- table %>% group_by(Group1, Group2) %>% summarise(AUC_mean = mean(value), AUC_sd = sd(value))
    table$Group1 <- factor(table$Group1, levels = c("Original", "Niche", "Phylogeny", "Random"))
    table$Group1 <- as.numeric(table$Group1)
    
    p_replace[[i]] <- table %>% ggplot(aes(x = Group1, y = AUC_mean)) + 
        geom_errorbar(aes(ymin = AUC_mean - AUC_sd, ymax = AUC_mean + AUC_sd), width = 0.1) + 
        geom_line(aes(group = Group2, color = Group2),  linewidth = 1, linetype = "dashed") + geom_point() +
        geom_smooth(level= 0.8, linewidth = 2) +
        scale_x_continuous(breaks = c(1, 2, 3, 4),  # 指定刻度位置  
                           labels = c("Original", "Niche", "Phylogeny", "Random"))+
        theme_bw(base_size = 14) +
        labs(x="", y="AUC", color="Study") +
        theme(text = element_text(face = "bold"), # 设置所有文本为加粗
              plot.title = element_text(hjust = 0.5), # 居中对齐标题
              axis.title = element_text(face = "bold"), # 加粗坐标轴标题
              axis.text.x = element_text(face = "bold", angle = 15, hjust = 1, vjust = 1)) +
        scale_color_brewer(palette = "Paired")
    
    i = i + 1
}

upp_plot <- plot_grid(IBD_p[[1]], cross_IBD[[1]],
                      CRC_p[[1]], cross_CRC[[1]],
               labels = c('a', 'c', 
                          'b', 'd'),
               align="hv",
               scale = c(1, 1, 
                         1, 1),
               nrow = 2, ncol=2, plot=FALSE, rel_widths = c(1, 0.8))
middle_plot <- plot_grid(p1, p2,
                         labels = c('e', 'f'),
                         align="hv",
                         scale = c(1, 1, 
                                   1, 1),
                         nrow = 1, ncol=2, plot=FALSE, rel_widths = c(1, 1))
    
bottom_plot <- plot_grid(p3, p_replace[[1]], p_replace[[2]],
                         labels = c('g', 'h', 'i'),
                         scale = c(1, 1, 1),
                         nrow = 1, ncol=3, rel_widths = c(0.8, 1, 1))

p <- plot_grid(upp_plot, middle_plot, bottom_plot,
               align="hv",
               scale = c(1, 1, 1),
               nrow = 3, ncol=1, plot=FALSE, rel_heights = c(4, 2, 2))

ggsave("/home/dongbiao/word_embedding_microbiome/result/Figure_model_performance.png", p,
       width = 50, height = 50, units = "cm")

