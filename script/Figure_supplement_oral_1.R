library(dplyr)
library(ggplot2)
library(lemon)
library(cowplot)
library(stringr)

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
    colnames(plot_data) <- c("value", "v_mean", "v_sd", "Group1", "Study")
    plot_data <- plot_data %>% arrange(desc(v_mean))
    plot_data$Group1 <- factor(plot_data$Group1, levels = unique(plot_data$Group1))
    p <- plot_data %>% ggplot(aes(x = Group1, y = v_mean / nums)) +
        labs(x="", y=paste0(AUC_or_F1), title=titles) +
        geom_col(fill = "skyblue", alpha=0.5) +
        geom_jitter(aes(x = Group1, y = value, color=Study), width = 0.1, 
                    height = 0, size = 2, alpha = 0.8) +
        geom_errorbar(aes(ymin = v_mean - v_sd, ymax = v_mean + v_sd), width = 0.2, color = "black", size = 0.7) +
        theme_bw(base_size = 14) +
        scale_color_brewer(palette = "Paired") +
        theme(text = element_text(face = "bold"), # 设置所有文本为加粗
              plot.title = element_text(hjust = 0.5), # 居中对齐标题
              axis.title = element_text(face = "bold"), # 加粗坐标轴标题
              axis.text = element_text(face = "bold")) + # 加粗坐标轴文本
        coord_cartesian(ylim = vector_x) +
        theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1))
    
    return(p)
}

plot_banchmark <- function(lable_1, attention_res, mlp_res, rf_res, 
                           num_study, study_name, vector_x_1, vector_x_2, titles){
    plot_list <- list()
    table <- get_auc_f1(attention_res, c(12, 18))
    table <- as.data.frame(lapply(table, as.numeric))
    colnames(table) <- c("AUC", "F1")
    table$Group1 <- c(rep(lable_1, 1, each=num_study))
    table$Study <- rep(study_name, 10, each=1)
    table_mean_sd <- table %>% group_by(Group1) %>% summarise(AUC_mean = mean(AUC), AUC_sd = sd(AUC),
                                                                      F1_mean = mean(F1), F1_sd = sd(F1))
    table <- merge(table, table_mean_sd, by = c("Group1"))
    plot_list[[1]] <- lou_plot(plot_data = table[,c("AUC", "AUC_mean", "AUC_sd", "Group1", "Study")], 
                               AUC_or_F1 = "AUC", 
                               vector_x = vector_x_1, nums = num_study, titles=titles)
    
    plot_list[[2]] <- lou_plot(plot_data = table[,c("F1", "F1_mean", "F1_sd", "Group1", "Study")], 
                               AUC_or_F1 = "F1", 
                               vector_x = vector_x_2, nums = num_study, titles=titles)
    
    return(plot_list)
}

# Banchmark
lable_1 <- c("Russell_rao_weight", "Russell_rao", "faith", "Jaccard", "Abundance-percentile", 
             "Abundance-totalsum", "Braycurtis-percentile", "Braycurtis-totalsum", 
             "Phylogeny", "OTU")
# caries
attention_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/caries/attention.out"
mlp_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/caries/mlp.out"
rf_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/caries/randownforest.out"
Caries_p <- plot_banchmark(lable_1, attention_res, mlp_res, rf_res, 
                           num_study=6, 
                           study_name=c('PRJNA330533', 'PRJNA383868', 'PRJNA454811', 'PRJNA480252',
                                        'PRJNA555320', 'PRJNA681486'),
                           vector_x_1=c(0.4, 1), vector_x_2=c(0.18, 0.84), titles="Caries")




#OSCC
attention_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/OSCC/attention.out"
mlp_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/OSCC/mlp.out"
rf_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/OSCC/randownforest.out"
OSCC_p <- plot_banchmark(lable_1, attention_res, mlp_res, rf_res, 
                         num_study=8, 
                         study_name=c('PRJNA756784', 'PRJNA751046', 'PRJEB39064', 'PRJNA700849',
                                      'OEP000837', 'PRJNA421234', 'PRJNA386665', 'PRJNA744870'),
                         vector_x_1=c(0.4, 1), vector_x_2=c(0.18, 0.84), titles="OSCC")



#periodontities
attention_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/paradentities/attention.out"
mlp_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/paradentities/mlp.out"
rf_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/paradentities/randownforest.out"
Periodontities_p <- plot_banchmark(lable_1, attention_res, mlp_res, rf_res, 
                                   num_study=6, 
                                   study_name=c('PRJNA321534', 'PRJEB61123', 'PRJNA578492', 'PRJNA601054','PRJNA843376','PRJEB42371'),
                                   vector_x_1=c(0.32, 1), vector_x_2=c(0.14, 0.88), titles="Periodontities")

p <- plot_grid(Caries_p[[1]],OSCC_p[[1]],Periodontities_p[[1]],
               labels = c('a', 'b','c'),
               align="hv",
               scale = c(1, 1),
               nrow = 1, ncol=3, plot=FALSE, rel_widths = c(1, 1))

ggsave("/home/dongbiao/word_embedding_microbiome/result/supplement_figure_oral_1.png", p,
       width = 30, height = 12, units = "cm")

