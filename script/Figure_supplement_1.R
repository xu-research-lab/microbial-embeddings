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
    colnames(plot_data) <- c("value", "v_mean", "v_sd", "Group1", "Group2", "Study")
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
    table <- table %>% filter(Group2 == "Attention")
    table_mean_sd <- table %>% group_by(Group1, Group2) %>% summarise(AUC_mean = mean(AUC), AUC_sd = sd(AUC),
                                                                      F1_mean = mean(F1), F1_sd = sd(F1))
    table <- merge(table, table_mean_sd, by = c("Group1", "Group2"))
    plot_list[[1]] <- lou_plot(plot_data = table[,c("AUC", "AUC_mean", "AUC_sd", "Group1", "Group2", "Study")], 
                               AUC_or_F1 = "AUC", 
                               vector_x = vector_x_1, nums = num_study, titles=titles)
    
    plot_list[[2]] <- lou_plot(plot_data = table[,c("F1", "F1_mean", "F1_sd", "Group1", "Group2", "Study")], 
                                 AUC_or_F1 = "F1", 
                                 vector_x = vector_x_2, nums = num_study, titles=titles)
    
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


p <- plot_grid(IBD_p[[1]], CRC_p[[1]],
               labels = c('a', 'b'),
               align="hv",
               scale = c(1, 1),
               nrow = 1, ncol=2, plot=FALSE, rel_widths = c(1, 1))

ggsave("/home/dongbiao/word_embedding_microbiome/result/supplement_figure_1.png", p,
       width = 30, height = 12, units = "cm")

