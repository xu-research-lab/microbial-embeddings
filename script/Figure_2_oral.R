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
    table$Group2 <- rep(c("Attention", "MLP", "RF"), 1, each= num_study * 10)
    table$Study <- rep(study_name, 10, each=1)
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
    cross_diff_size$Study <- rep(study_name, 6, each=10)
    cross_diff_size$datasize <- rep(c("5000", "10000", "20000", "40000", "80000","120000"), 1, 
                                    each=10 * num_study)
    line_plot <- cross_diff_size %>% group_by(Study, datasize) %>% 
        summarise(AUC_mean = mean(AUC), AUC_sd = sd(AUC),
                  F1_mean = mean(F1), F1_sd = sd(F1))
    cross_diff_size <- merge(cross_diff_size, line_plot, by=c("Study", "datasize"))
    cross_diff_size$datasize <- factor(cross_diff_size$datasize, 
                                       levels=c("5000", "10000", "20000", "40000", "80000","120000"))
    
    plot_list <- list()
    cross_diff_size$datasize <- as.numeric(cross_diff_size$datasize)
    plot_list[[1]] <- cross_diff_size %>% ggplot(aes(x = datasize, y = AUC)) + 
        geom_errorbar(aes(ymin = AUC_mean - AUC_sd, ymax = AUC_mean + AUC_sd), width = 0.1) + 
        geom_point(aes(color = Study), size = 2, alpha=0.8) +
        geom_line(aes(group = Study, y = AUC_mean, color = Study), alpha = 0.5, linewidth = 1, linetype = "dashed") +
        labs(x="Datasize", y="AUC", title=titles) +
        geom_smooth(level= 0.9) +
        theme_bw(base_size = 14) +
        scale_color_brewer(palette = "Paired") +
        scale_x_continuous(breaks = c(1, 2, 3, 4, 5,6),  # 指定刻度位置  
                           labels = c("5000", "10000", "20000", "40000", "80000","120000"))+
        theme(text = element_text(face = "bold"), # 设置所有文本为加粗
              plot.title = element_text(hjust = 0.5), # 居中对齐标题
              axis.title = element_text(face = "bold"), # 加粗坐标轴标题
              axis.text = element_text(face = "bold"))
    
    plot_list[[2]] <- cross_diff_size %>% ggplot(aes(x = datasize, y = F1)) +
        geom_boxplot() + geom_point(aes(color = Study), size = 2, alpha=0.8) +
        geom_errorbar(aes(ymin = F1_mean - F1_sd, ymax = F1_mean + F1_sd), width = 0.1) + 
        geom_line(aes(group = Study, y = F1_mean, color = Study), alpha = 0.5, linewidth = 1) +
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
             "Phylogeny", "OTU")
# caries
attention_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/caries/attention_oral.out"
mlp_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/caries/mlp_oral.out"
rf_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/caries/randownforest_oral.out"
Caries_p <- plot_banchmark(lable_1, attention_res, mlp_res, rf_res, 
                        num_study=6, 
                        study_name=c('PRJNA330533', 'PRJNA383868', 'PRJNA454811', 'PRJNA480252',
                                     'PRJNA555320', 'PRJNA681486'),
                        vector_x_1=c(0.4, 1), vector_x_2=c(0.18, 0.84), titles="Caries")




#OSCC
attention_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/OSCC/attention_OSCC.out"
mlp_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/OSCC/mlp_oral.out"
rf_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/OSCC/randownforest_oral.out"
OSCC_p <- plot_banchmark(lable_1, attention_res, mlp_res, rf_res, 
                           num_study=8, 
                           study_name=c('PRJNA756784', 'PRJNA751046', 'PRJEB39064', 'PRJNA700849',
                                        'OEP000837', 'PRJNA421234', 'PRJNA386665', 'PRJNA744870'),
                           vector_x_1=c(0.4, 1), vector_x_2=c(0.18, 0.84), titles="OSCC")



#periodontities
attention_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/paradentities/attention_oral.out"
mlp_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/paradentities/mlp_oral.out"
rf_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/paradentities/randownforest_oral.out"
Periodontities_p <- plot_banchmark(lable_1, attention_res, mlp_res, rf_res, 
                        num_study=6, 
                        study_name=c('PRJNA321534', 'PRJEB61123', 'PRJNA578492', 'PRJNA601054','PRJNA843376','PRJEB42371'),
                        vector_x_1=c(0.32, 1), vector_x_2=c(0.14, 0.88), titles="Periodontities")



# # Caries
cross_caries <- plot_diff_dataset(cross_res="/home/dongbiao/word_embedding_microbiome/programe_test/filter_prevalence/sub_embedding/caries_cross/attention.out",
                               study_name=c('PRJNA330533', 'PRJNA383868', 'PRJNA454811', 'PRJNA480252',
                                            'PRJNA555320', 'PRJNA681486'),
                               num_study=6, titles="Caries")

# # OSCC
cross_OSCC <- plot_diff_dataset(cross_res="/home/dongbiao/word_embedding_microbiome/programe_test/filter_prevalence/sub_embedding/oscc_cross/attention.out",
                               study_name=c('PRJNA756784', 'PRJNA751046', 'PRJEB39064', 'PRJNA700849',
                                            'OEP000837', 'PRJNA421234', 'PRJNA386665', 'PRJNA744870'),
                               num_study=8, titles="OSCC")

#
# # paradentities
cross_periodontities <- plot_diff_dataset(cross_res="/home/dongbiao/word_embedding_microbiome/programe_test/filter_prevalence/sub_embedding/paradentities_cross/attention.out",
                                study_name=c('PRJNA321534', 'PRJEB61123', 'PRJNA578492', 'PRJNA601054','PRJNA843376','PRJEB42371'),
                                num_study=6, titles="Periodontities")



# ### embedding remove batch effects
load(file = "/home/dongbiao/word_embedding_microbiome/programe_test/caries/explain/adonise2_res.RData")
anosim_caries <- plot_data
load(file = "/home/dongbiao/word_embedding_microbiome/programe_test/OSCC/explain/adonise2_res.RData")
anosim_OSCC <- plot_data
load(file = "/home/dongbiao/word_embedding_microbiome/programe_test/paradentities/explain/adonise2_res.RData")
anosim_periodontities <- plot_data
anosim_res <- rbind(anosim_caries, anosim_OSCC,anosim_periodontities)
anosim_res <- melt(anosim_res, id.vars=c("Effects", "group"))
anosim_res$group <- factor(anosim_res$group, levels = c("Caries", "OSCC","Peridontities"))
anosim_res$group <-c("Caries","Caries","OSCC","OSCC","Periodontities","Periodontities","Caries","Caries","OSCC","OSCC","Periodontities","Periodontities","Caries","Caries","OSCC","OSCC","Periodontities","Periodontities")
p1 <- anosim_res %>% ggplot(aes(x=variable, y=`value`, fill=Effects)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~group, scales = "free")+
    theme_bw() +
    labs(x="", y="R2", fill = "") +
    theme_bw(base_size = 14) +
    scale_color_brewer(palette = "Set2") +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"),
          legend.position = "top",
          axis.text.x = element_text(angle = -10, hjust = 0, vjust = 1))

# ### Attention weight
attention_weight_caries <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/caries/explain/attention_weight.csv")
attention_weight_OSCC <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/OSCC/explain/attention_weight.csv")
attention_weight_periodontities <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/paradentities/explain/attention_weight.csv")
attention_weight_caries$group3 <- rep("Caries", n=nrow(attention_weight_caries))
attention_weight_OSCC$group3 <- rep("OSCC", n=nrow(attention_weight_OSCC))
attention_weight_periodontities$group3 <- rep("Periodontities", n=nrow(attention_weight_periodontities))
attention_weight <- rbind(attention_weight_caries, attention_weight_OSCC,attention_weight_periodontities)
attention_weight$group3 <- factor(attention_weight$group3, c("Caries", "OSCC","Periodontities"))
attention_weight <- attention_weight %>% group_by(group1, group2, group3) %>% summarise(attention_weight=mean(attention_weight))

p2 <- attention_weight %>% ggplot(aes(x=group2, y=attention_weight, fill=group1)) +
    geom_bar(stat='identity', position = 'dodge') +
    facet_grid(~group3, scales = "free", space = "free_x") +
    theme_bw(base_size = 14)+
    theme(text = element_text(face = "bold"),
          axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),
          legend.text = element_text(size = 14, color = "black", face = "bold"),
          legend.box.margin = margin(l = 0.1, unit = "in"),
          legend.position = "top",
          legend.title.align = 0.5
    ) +
    labs(x="", y="Average of attention weight", fill="")

niche_otu_cos <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/niche/niche_OTU_cos_oral_gut.csv")

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
studys[["caries"]] <- c('PRJNA330533', 'PRJNA383868', 'PRJNA454811', 'PRJNA480252', 'PRJNA555320', 'PRJNA681486')
studys[["OSCC"]] <- c('PRJNA756784', 'PRJNA751046', 'PRJEB39064', 'PRJNA700849',
                      'OEP000837', 'PRJNA421234', 'PRJNA386665', 'PRJNA744870')
studys[["paradentities"]] <- c('PRJNA321534', 'PRJEB61123', 'PRJNA578492', 'PRJNA601054','PRJNA843376','PRJEB42371')
types <- c("caries", "OSCC", "paradentities")
p_replace <- list()
get_numbirc <- function(vectors){
    vectors <- gsub("\\[|\\]", "", vectors)
    vectors <- as.numeric(unlist(strsplit(gsub(" ", "", vectors), ",")))
    return(vectors)
}
i <- 1
title_list <- c("Caries", "OSCC", "Periodontities")
num_studies <- c(6, 8, 6)
for (type in types){
    m <- num_studies[i]
    table = read.csv(paste0("/home/dongbiao/word_embedding_microbiome/programe_test/", type, "/niche/res_niche.csv")) %>% t() %>% as.data.frame()
    Niche <- c()
    Phylogeny <- c()
    Random <- c()
    for (n in c(1:m)){
        Niche <- c(Niche, get_numbirc(table[n, 2]))
        Phylogeny <- c(Phylogeny, get_numbirc(table[n, 3]))
        Random <- c(Random, get_numbirc(table[n, 4]))
    }
    value <- c(table[,1], Niche, Phylogeny, Random)
    value <- as.numeric(value)
    Group1 <- c(rep("Original", m), rep("Niche", 10 * m), rep("Phylogeny", 10 * m), rep("Random", 10 * m))
    Group2 <- c(studys[[type]], rep(rep(studys[[type]], n=m, each=10), 3))
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
        labs(x="", y="AUC", color="Study", title=title_list[i]) +
        theme(text = element_text(face = "bold"), # 设置所有文本为加粗
              plot.title = element_text(hjust = 0.5), # 居中对齐标题
              axis.title = element_text(face = "bold"), # 加粗坐标轴标题
              axis.text.x = element_text(face = "bold", angle = 15, hjust = 1, vjust = 1)) +
        scale_color_brewer(palette = "Paired")
    
    i = i + 1
}

upp_plot <- plot_grid(Caries_p[[1]], cross_caries[[1]],
                      OSCC_p[[1]], cross_OSCC[[1]],
                      Periodontities_p[[1]], cross_periodontities[[1]],
                      labels = c('a', 'd', 
                                 'b', 'e',
                                 'c', 'f'),
                      align="hv",
                      nrow = 3, ncol=2, plot=FALSE, rel_widths = c(1, 0.8))

middle_plot <- plot_grid(p1, p2, p3, 
                         labels = c('g', 'h'),
                         align="hv",
                         nrow = 1, ncol=3, plot=FALSE, rel_widths = c(1, 1, 0.7))

bottom_plot <- plot_grid(p_replace[[1]], p_replace[[2]], p_replace[[3]],
                         labels = c('i', 'j', 'k'),
                         scale = c(1, 1, 1),
                         nrow = 1, ncol=3, rel_widths = c(1, 1, 1))

p <- plot_grid(upp_plot, middle_plot, bottom_plot,
               align="hv",
               scale = c(1, 1, 1),
               nrow = 3, ncol=1, plot=FALSE, rel_heights = c(4, 1.8, 2))

ggsave("/home/dongbiao/word_embedding_microbiome/result/Figure_model_performance_oral.png", p,
       width = 50, height = 50, units = "cm")

p <- plot_grid(Caries_p[[1]], OSCC_p[[1]], Periodontities_p[[1]],
               nrow = 3, ncol=1)

ggsave("/home/dongbiao/word_embedding_microbiome/result/performance_oral.png", p,
       width = 25, height = 25, units = "cm")
