library(dplyr)
library(ggplot2)
library(lemon)
library(cowplot)
library(stringr)

get_auc_f1 <- function(file, idx){
    table <- read.table(file, sep = "\t")
    table <- table[apply(table, 1, function(x) grepl("mcc", x)), ]
    table <- gsub(",", "", table)
    table <- str_split(table, " ", simplify = TRUE)
    table <- as.data.frame(table[, idx])
}
lable_1 <- c("Russell_rao_weight", "Russell_rao", "faith", "Jaccard", "Abundance-percentile", 
             "Abundance-totalsum", "Braycurtis-percentile", "Braycurtis-totalsum", 
             "Phylogeny", "PCA", "RandInit")
lable_2 <- c("Russell_rao_weight", "Russell_rao", "faith", "Jaccard", "Abundance-percentile", 
             "Abundance-totalsum", "Braycurtis-percentile", "Braycurtis-totalsum", 
             "Phylogeny", "PCA", "OTU")
train_data <- "_all_feces"
# train_data <- "" # 表示预训练在AGP数据集中进行的
group2 <- c("Attention", "MLP", "RF")

# AGP_table_1 <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/filter_prevalence/GloVe/AGP_IBD/attention_AGP.out", c(12, 18))
AGP_table_1 <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/filter_prevalence/all_feces/AGP_IBD/attention_1.out", c(12, 18))
AGP_table_2 <- get_auc_f1(paste0("/home/dongbiao/word_embedding_microbiome/programe_test/AGP_run_20220702/mlp", train_data,".out"), c(12, 18))
AGP_table_3 <- get_auc_f1(paste0("/home/dongbiao/word_embedding_microbiome/programe_test/AGP_run_20220702/randownforest", train_data,".out"), c(4, 12))
AGP_table <- rbind(AGP_table_1, AGP_table_2, AGP_table_3)
AGP_table <- as.data.frame(lapply(AGP_table, as.numeric))
colnames(AGP_table) <- c("AUC", "F1")
AGP_table$Group1 <- c(rep(lable_1, 1, each=5),
                      rep(lable_2, 1, each=5),
                      rep(lable_2, 1, each=5))
AGP_table$Group2 <- rep(c("Attention", "MLP", "RF"), 1, each=55)


# plot for five-fold data in AGP
AGP_table_mean_sd <- AGP_table %>% group_by(Group1, Group2) %>% summarise(AUC_mean = mean(AUC), AUC_sd = sd(AUC),
                                                        F1_mean = mean(F1), F1_sd = sd(F1))
AGP_table <- merge(AGP_table, AGP_table_mean_sd, by = c("Group1", "Group2"))

five_fold_plot <- function(plot_data, AUC_or_F1, classify, vector_x){
    colnames(plot_data) <- c("value", "v_mean", "v_sd", "Group1", "Group2")
    data <- plot_data %>% filter(Group2 == classify)
    data <- data %>% arrange(desc(v_mean))
    data$Group1 <- factor(data$Group1, levels = unique(data$Group1))
    p <- data %>% ggplot(aes(x = Group1, y = v_mean / 5)) +
        labs(x="", y=paste0(AUC_or_F1)) +
        geom_col(fill = "skyblue", alpha=0.5) +
        geom_jitter(aes(x = Group1, y = value), width = 0.1, height = 0, size = 2) +
        geom_errorbar(aes(ymin = v_mean - v_sd, ymax = v_mean + v_sd), width = 0.2, color = "black", size = 0.7) +
        ggtitle(classify) +
        theme_bw(base_size = 14) +
        theme(text = element_text(face = "bold"), # 设置所有文本为加粗
              plot.title = element_text(hjust = 0.5), # 居中对齐标题
              axis.title = element_text(face = "bold"), # 加粗坐标轴标题
              axis.text = element_text(face = "bold")) + # 加粗坐标轴文本
        coord_cartesian(ylim = vector_x) +
        theme(axis.text.x = element_text(angle = 65, hjust = 1, vjust = 1))
    return(p)
}

colnames(plot_data) <- c("value", "v_mean", "v_sd", "Group1", "Group2")

p <- IBD_table %>% ggplot(aes(x = Group2, y = AUC_mean / 4)) +
    geom_col(aes(fill = Group2), alpha=0.5, width = 1) +
    geom_jitter(aes(x = Group2, y = AUC, color = Study), width = 0.1, height = 0, size = 2) +
    geom_errorbar(aes(ymin = AUC_mean - AUC_sd, ymax = AUC_mean + AUC_sd), 
                  width = 0.2, color = "black", size = 0.7) +
    theme_minimal(base_size = 14) +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"),
          strip.text.x = element_text(size=12, angle=65, hjust = 1, vjust = 1)) + # 加粗坐标轴文本
    scale_x_discrete(labels = c('','','')) +
    labs(x="", y="AUC") +
    facet_wrap(vars(Group1),nrow = 1, strip.position = "bottom") + 
    scale_color_discrete(guide = FALSE)


plot_list <- list()
for (i in 1 : length(group2)){
    plot_list[[i]] <- five_fold_plot(classify = group2[i],
                        plot_data = AGP_table[,c("AUC", "AUC_mean", "AUC_sd", "Group1", "Group2")], 
                        AUC_or_F1 = "AUC", 
                        vector_x = c(0.62, 0.90))
}
AGP_res_file <- "/home/dongbiao/word_embedding_microbiome/programe_test/AGP_run_20220702/results/"
bar_plot <- plot_list %>% plot_grid(plotlist = ., align = "h", nrow = 1)
ggsave(paste0(AGP_res_file, paste0("AGP_AUC", train_data, ".png")), bar_plot,
       width = 20, height = 10, units = "cm")
plot_list <- list()

for (i in 1 : length(group2)){
    plot_list[[i]] <- five_fold_plot(classify = group2[i],
                                     plot_data = AGP_table[,c("F1", "F1_mean", "F1_sd", "Group1", "Group2")], 
                                     AUC_or_F1 = "F1", 
                                     vector_x = c(0.45, 0.82))
}
bar_plot <- plot_list %>% plot_grid(plotlist = ., align = "h", nrow = 1)
ggsave(paste0(AGP_res_file, paste0("AGP_F1", train_data, ".png")), bar_plot,
       width = 20, height = 10, units = "cm")

# get cross data results
# IBD_table_1 <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/filter_prevalence/GloVe/IBD_cross/attention_AGP.out", c(12, 18))
IBD_table_1 <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/filter_prevalence/all_feces/IBD_cross/attention_1.out", c(12, 18))
IBD_table_2 <- get_auc_f1(paste0("/home/dongbiao/word_embedding_microbiome/programe_test/IBD/mlp", train_data, ".out"), c(12, 18))
IBD_table_3 <- get_auc_f1(paste0("/home/dongbiao/word_embedding_microbiome/programe_test/IBD/randownforest", train_data, ".out"), c(4, 12))
IBD_table <- rbind(IBD_table_1, IBD_table_2, IBD_table_3)
IBD_table <- as.data.frame(lapply(IBD_table, as.numeric))
colnames(IBD_table) <- c("AUC", "F1")
IBD_table$Group1 <- c(rep(lable_1, 1, each=4),
                      rep(lable_2, 1, each=4),
                      rep(lable_2, 1, each=4))
IBD_table$Group2 <- rep(c("Attention", "MLP", "RF"), 1, each=44)
IBD_table$Study <- rep(c('qiita_1629', 'qiita_2538', 'PRJNA422193', 'PRJNA431126'), 33, each=1)
IBD_table_mean_sd <- IBD_table %>% group_by(Group1, Group2) %>% summarise(AUC_mean = mean(AUC), AUC_sd = sd(AUC),
                                                                          F1_mean = mean(F1), F1_sd = sd(F1))
IBD_table <- merge(IBD_table, IBD_table_mean_sd, by = c("Group1", "Group2"))


# CRC_table_1 <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/filter_prevalence/GloVe/CRC_cross/attention_AGP.out", c(12, 18))
CRC_table_1 <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/filter_prevalence/all_feces/CRC_cross/attention_1.out", c(12, 18))
CRC_table_2 <- get_auc_f1(paste0("/home/dongbiao/word_embedding_microbiome/programe_test/CRC/mlp", train_data, ".out"), c(12, 18))
CRC_table_3 <- get_auc_f1(paste0("/home/dongbiao/word_embedding_microbiome/programe_test/CRC/randownforest", train_data, ".out"), c(4, 12))
CRC_table <- rbind(CRC_table_1, CRC_table_2, CRC_table_3)
CRC_table <- as.data.frame(lapply(CRC_table, as.numeric))
colnames(CRC_table) <- c("AUC", "F1")
CRC_table$Group1 <- c(rep(lable_1, 1, each=7),
                      rep(lable_2, 1, each=7),
                      rep(lable_2, 1, each=7))
CRC_table$Group2 <- rep(c("Attention", "MLP", "RF"), 1, each=77)
CRC_table$Study <- rep(c('PRJEB36789', 'PRJNA824020', 'PRJDB11845', 'PRJNA318004',
                         'PRJEB6070', 'PRJNA430990', 'PRJNA290926'), 33, each=1)
# CRC_table <- CRC_table %>% filter(Study != "PRJNA763872")
CRC_table_mean_sd <- CRC_table %>% group_by(Group1, Group2) %>% summarise(AUC_mean = mean(AUC), AUC_sd = sd(AUC),
                                                                          F1_mean = mean(F1), F1_sd = sd(F1))
CRC_table <- merge(CRC_table, CRC_table_mean_sd, by = c("Group1", "Group2"))


# dietary_fber_table_1 <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/filter_prevalence/GloVe/dietary_cross/attention_AGP_1.out", c(12, 18))
dietary_fber_table_1 <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/filter_prevalence/all_feces/dietary_cross/attention_1.out", c(12, 18))
dietary_fber_table_2 <- get_auc_f1(paste0("/home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/mlp", train_data, ".out"), c(12, 18))
dietary_fber_table_3 <- get_auc_f1(paste0("/home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/randownforest", train_data, ".out"), c(4, 12))
dietary_fber_table <- rbind(dietary_fber_table_1, dietary_fber_table_2, dietary_fber_table_3)
dietary_fber_table <- as.data.frame(lapply(dietary_fber_table, as.numeric))
colnames(dietary_fber_table) <- c("AUC", "F1")
dietary_fber_table$Group1 <- c(rep(lable_1, 1, each=6),
                               rep(lable_2, 1, each=6),
                               rep(lable_2, 1, each=6))
dietary_fber_table$Group2 <- rep(c("Attention", "MLP", "RF"), 1, each=66)
dietary_fber_table$Study <- rep(c('qiita_13367', 'SRP067761', 'SRP120250', 'SRP128128', 'SRP345891', 'SRP219296'), 
                                33, each=1)
dietary_fber_table_mean_sd <- dietary_fber_table %>% group_by(Group1, Group2) %>% summarise(AUC_mean = mean(AUC), AUC_sd = sd(AUC),
                                                                          F1_mean = mean(F1), F1_sd = sd(F1))
dietary_fber_table <- merge(dietary_fber_table, dietary_fber_table_mean_sd, by = c("Group1", "Group2"))


# plot for corss data
lou_auc <- function(plot_data, AUC_or_F1, classify, vector_x, nums){
    colnames(plot_data) <- c("value", "v_mean", "v_sd", "Group1", "Group2", "Study")
    data <- plot_data %>% filter(Group2 == classify)
    data <- data %>% arrange(desc(v_mean))
    data$Group1 <- factor(data$Group1, levels = unique(data$Group1))
    p <- data %>% ggplot(aes(x = Group1, y = v_mean / nums)) +
        labs(x="", y=paste0(AUC_or_F1)) +
        geom_col(fill = "skyblue", alpha=0.5) +
        geom_jitter(aes(x = Group1, y = value, color=Study), width = 0.1, height = 0, size = 2, alpha = 0.8) +
        geom_errorbar(aes(ymin = v_mean - v_sd, ymax = v_mean + v_sd), width = 0.2, color = "black", size = 0.7) +
        ggtitle(classify) +
        theme_bw(base_size = 14) +
        scale_color_brewer(palette = "Set3") +
        theme(text = element_text(face = "bold"), # 设置所有文本为加粗
              plot.title = element_text(hjust = 0.5), # 居中对齐标题
              axis.title = element_text(face = "bold"), # 加粗坐标轴标题
              axis.text = element_text(face = "bold")) + # 加粗坐标轴文本
        coord_cartesian(ylim = vector_x) +
        theme(axis.text.x = element_text(angle = 65, hjust = 1, vjust = 1))
    
    return(p)
}

# IBD results
IDB_res_file <- "/home/dongbiao/word_embedding_microbiome/programe_test/IBD/results/"

for (i in 1 : length(group2)){
    plot_list[[i]] <- lou_auc(classify = group2[i],
                              plot_data = IBD_table[,c("AUC", "AUC_mean", "AUC_sd", "Group1", "Group2", "Study")], 
                              AUC_or_F1 = "AUC", 
                              vector_x = c(0.4, 1), nums = 4)
}
p<- grid_arrange_shared_legend(plot_list[[1]], plot_list[[2]],
                               plot_list[[3]],
                               nrow = 1, position='bottom')
ggsave(paste0(IDB_res_file, paste0("AUC_cross_data", train_data, ".png")), p,
       width = 20, height = 10, units = "cm")

for (i in 1 : length(group2)){
    plot_list[[i]] <- lou_auc(classify = group2[i],
                              plot_data = IBD_table[,c("F1", "F1_mean", "F1_sd", "Group1", "Group2", "Study")], 
                              AUC_or_F1 = "F1", 
                              vector_x = c(0.18, 0.84), nums = 4)
}
p<- grid_arrange_shared_legend(plot_list[[1]], plot_list[[2]],
                               plot_list[[3]],
                               nrow = 1, position='bottom')
ggsave(paste0(IDB_res_file, paste0("F1_cross_data", train_data, ".png")), p,
       width = 20, height = 10, units = "cm")


# CRC results
CRC_res_file <- "/home/dongbiao/word_embedding_microbiome/programe_test/CRC/results/"
for (i in 1 : length(group2)){
    plot_list[[i]] <- lou_auc(classify = group2[i],
                              plot_data = CRC_table[,c("AUC", "AUC_mean", "AUC_sd", "Group1", "Group2", "Study")], 
                              AUC_or_F1 = "AUC", 
                              vector_x = c(0.32, 0.95), nums = 7)
}
p<- grid_arrange_shared_legend(plot_list[[1]], plot_list[[2]],
                               plot_list[[3]],
                               nrow = 1, position='bottom')
ggsave(paste0(CRC_res_file, paste0("AUC_cross_data_add_study", train_data, ".png")), p,
       width = 20, height = 10, units = "cm")

for (i in 1 : length(group2)){
    plot_list[[i]] <- lou_auc(classify = group2[i],
                              plot_data = CRC_table[,c("F1", "F1_mean", "F1_sd", "Group1", "Group2", "Study")], 
                              AUC_or_F1 = "F1", 
                              vector_x = c(0.14, 0.88), nums = 7)
}
p<- grid_arrange_shared_legend(plot_list[[1]], plot_list[[2]],
                               plot_list[[3]],
                               nrow = 1, position='bottom')
ggsave(paste0(CRC_res_file, paste0("F1_cross_data_add_study", train_data, ".png")), p,
       width = 20, height = 10, units = "cm")

# dietary results
dietay_res_file <- "/home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/results/"
for (i in 1 : length(group2)){
    plot_list[[i]] <- lou_auc(classify = group2[i],
                              plot_data = dietary_fber_table[,c("AUC", "AUC_mean", "AUC_sd", "Group1", "Group2", "Study")], 
                              AUC_or_F1 = "AUC", 
                              vector_x = c(0.30, 0.9), nums = 6)
}
p<- grid_arrange_shared_legend(plot_list[[1]], plot_list[[2]],
                               plot_list[[3]],
                               nrow = 1, position='bottom')
ggsave(paste0(dietay_res_file, paste0("AUC_cross_data", train_data,".png")), p,
       width = 20, height = 10, units = "cm")

for (i in 1 : length(group2)){
    plot_list[[i]] <- lou_auc(classify = group2[i],
                              plot_data = dietary_fber_table[,c("F1", "F1_mean", "F1_sd", "Group1", "Group2", "Study")], 
                              AUC_or_F1 = "F1", 
                              vector_x = c(0.16, 0.85), nums = 6)
}
p<- grid_arrange_shared_legend(plot_list[[1]], plot_list[[2]],
                               plot_list[[3]],
                               nrow = 1, position='bottom')
ggsave(paste0(dietay_res_file, paste0("F1_cross_data", train_data, ".png")), p,
       width = 20, height = 10, units = "cm")

# compare the results of all feces with AGP
# IBD in AGP datasets
lable_3 <- c("Russell_rao_weight", "Russell_rao", "faith", "Jaccard", "Abundance-percentile", 
             "Abundance-totalsum", "Braycurtis-percentile", "Braycurtis-totalsum")
AGP_table_4 <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/AGP_run_20220702/attention_all_feces.out", c(12, 18))
AGP_table_4 <- as.data.frame(lapply(AGP_table_4[1:40, ], as.numeric))
colnames(AGP_table_4) <- c("AUC", "F1")
AGP_table_4$Group1 <- c(rep(c("Russell_rao_weight", "Russell_rao", "faith", "Jaccard", "Abundance-percentile", 
                              "Abundance-totalsum", "Braycurtis-percentile", "Braycurtis-totalsum"), 1, each=5))
AGP_table_4$Group2 <- rep(c("All feces sample"), 1, each=40)
AGP_diff_datasets <- rbind(AGP_table[AGP_table$Group1 %in% lable_3,] %>% filter(Group2 == "Attention") %>% 
    select(c("AUC", "F1", "Group1")) %>%
    mutate(Group2 = rep(c("AGP feces sample"), 1, each=40)), AGP_table_4)
p <- AGP_diff_datasets %>% ggplot(aes(x = Group2, y = AUC)) + 
    geom_boxplot() + facet_grid(.~ Group1, scales = "free_y") +
    labs(x="", y="AUC") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(paste0(AGP_res_file, "AGP_AUC_diff_datasets.png"), p,
       width = 20, height = 10, units = "cm")
p <- AGP_diff_datasets %>% ggplot(aes(x = Group2, y = F1)) + 
    geom_boxplot() + facet_grid(.~ Group1, scales = "free_y") +
    labs(x="", y="F1") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(paste0(AGP_res_file, "AGP_F1_diff_datasets.png"), p,
       width = 20, height = 10, units = "cm")

# IBD corss data
IBD_table_4 <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/IBD/attention_all_feces.out", c(12, 18))
IBD_table_4 <- as.data.frame(lapply(IBD_table_4[1:32,], as.numeric))
colnames(IBD_table_4) <- c("AUC", "F1")
IBD_table_4$Group1 <- c(rep(lable_3, 1, each=4))
IBD_table_4$Group2 <- rep(c("All feces sample"), 1, each=32)
IBD_table_4$Study <- rep(c('qiita_1629', 'qiita_2538', 'PRJNA422193', 'PRJNA431126'), 8, each=1)
IBD_diff_datasets <- rbind(IBD_table[IBD_table$Group1 %in% lable_3, ] %>% filter(Group2 == "Attention") %>% 
                               select(c("AUC", "F1", "Group1", "Study")) %>%
                               mutate(Group2 = rep(c("AGP feces sample"), 1, each=32)), IBD_table_4)
p <- IBD_diff_datasets %>% ggplot(aes(x = Group2, y = AUC)) + geom_boxplot() + geom_point(aes(color = Study), size = 3) +
    geom_line(aes(group = Study, y = AUC, color = Study), alpha = 0.5, linewidth = 2) +
    facet_grid(~ Group1, scales = "free_y") +
    labs(x="", y="AUC") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    scale_color_brewer(palette = "Set3")
ggsave(paste0(IDB_res_file, "AGP_AUC_diff_datasets.png"), p,
       width = 30, height = 10, units = "cm")
p <- IBD_diff_datasets %>% ggplot(aes(x = Group2, y = F1)) + geom_boxplot() + geom_point(aes(color = Study), size = 3) +
    geom_line(aes(group = Study, y = F1, color = Study), alpha = 0.5, linewidth = 2) +
    facet_grid(~ Group1, scales = "free_y") +
    labs(x="", y="F1") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    scale_color_brewer(palette = "Set3")
ggsave(paste0(IDB_res_file, "AGP_F1_diff_datasets.png"), p,
       width = 30, height = 10, units = "cm")

# CRC cross data
CRC_table_4 <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/CRC/attention_all_feces.out", c(12, 18))
CRC_table_4 <- as.data.frame(lapply(CRC_table_4[1:48, ], as.numeric))
colnames(CRC_table_4) <- c("AUC", "F1")
CRC_table_4$Group1 <- c(rep(lable_3, 1, each=6))
CRC_table_4$Group2 <- rep(c("All feces sample"), 1, each=48)
CRC_table_4$Study <- rep(c("PRJNA763872", "PRJEB36789", "PRJNA824020", "PRJDB11845", "PRJNA318004", "PRJEB6070"), 8, each=1)
CRC_diff_datasets <- rbind(CRC_table[CRC_table$Group1 %in% lable_3,] %>% filter(Group2 == "Attention") %>% 
                               select(c("AUC", "F1", "Group1", "Study")) %>%
                               mutate(Group2 = rep(c("AGP feces sample"), 1, each=48)), CRC_table_4)
p <- CRC_diff_datasets %>% ggplot(aes(x = Group2, y = AUC)) + geom_boxplot() + geom_point(aes(color = Study), size = 3) +
    geom_line(aes(group = Study, y = AUC, color = Study), alpha = 0.5, linewidth = 2) +
    facet_grid(~ Group1, scales = "free_y") +
    labs(x="", y="AUC") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    scale_color_brewer(palette = "Set3")
ggsave(paste0(CRC_res_file, "AGP_AUC_diff_datasets.png"), p,
       width = 30, height = 10, units = "cm")
p <- CRC_diff_datasets %>% ggplot(aes(x = Group2, y = F1)) + geom_boxplot() + geom_point(aes(color = Study), size = 3) +
    geom_line(aes(group = Study, y = F1, color = Study), alpha = 0.5, linewidth = 2) +
    facet_grid(~ Group1, scales = "free_y") +
    labs(x="", y="F1") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    scale_color_brewer(palette = "Set3")
ggsave(paste0(CRC_res_file, "AGP_F1_diff_datasets.png"), p,
       width = 30, height = 10, units = "cm")

# fiber
# dietary_fber_table_4 <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/attention_all_feces.out", c(12, 18))
dietary_fber_table_1 <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/filter_prevalence/all_feces/dietary_cross/attention.out", c(12, 18))
dietary_fber_table_4 <- as.data.frame(lapply(dietary_fber_table_4[1:40, ], as.numeric))
colnames(dietary_fber_table_4) <- c("AUC", "F1")
dietary_fber_table_4$Group1 <- c(rep(lable_3, 1, each=5))
dietary_fber_table_4$Group2 <- rep(c("All feces sample"), 1, each=40)
dietary_fber_table_4$Study <- rep(c('SRP067761', 'SRP120250', 'SRP219296', 'qiita_13367', 'SRP345891'), 8, each=1)
dietary_fber_diff_datasets <- rbind(dietary_fber_table[dietary_fber_table$Group1 %in% lable_3, ] %>% filter(Group2 == "Attention") %>% 
                               select(c("AUC", "F1", "Group1", "Study")) %>%
                               mutate(Group2 = rep(c("AGP feces sample"), 1, each=40)), dietary_fber_table_4)
p <- dietary_fber_diff_datasets %>% ggplot(aes(x = Group2, y = AUC)) + geom_boxplot() + geom_point(aes(color = Study), size = 3) +
    geom_line(aes(group = Study, y = AUC, color = Study), alpha = 0.5, linewidth = 2) +
    facet_grid(~ Group1, scales = "free_y") +
    labs(x="", y="AUC") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    scale_color_brewer(palette = "Set3")
ggsave(paste0(dietay_res_file, "AGP_AUC_diff_datasets.png"), p,
       width = 30, height = 10, units = "cm")
p <- dietary_fber_diff_datasets %>% ggplot(aes(x = Group2, y = F1)) + geom_boxplot() + geom_point(aes(color = Study), size = 3) +
    geom_line(aes(group = Study, y = F1, color = Study), alpha = 0.5, linewidth = 2) +
    facet_grid(~ Group1, scales = "free_y") +
    labs(x="", y="F1") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    scale_color_brewer(palette = "Set3")
ggsave(paste0(dietay_res_file, "AGP_F1_diff_datasets.png"), p,
       width = 35, height = 10, units = "cm")


## AGP difference dataset size
# AGP IBD
AGP_res_file <- "/home/dongbiao/word_embedding_microbiome/programe_test/AGP_run_20220702/results/"
AGP_IBD_diff_size <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/AGP_run_20220702/attention_AGP_diff_size.out", c(12, 18))
colnames(AGP_IBD_diff_size) <- c("AUC", "F1")
AGP_IBD_diff_size <- as.data.frame(lapply(AGP_IBD_diff_size, as.numeric))
AGP_IBD_diff_size$datasize <- rep(c("10", "30", "50", "70", "90"), 1, each=5)

p <- AGP_IBD_diff_size %>% ggplot(aes(x = datasize, y = AUC)) + geom_boxplot() + geom_point(size = 3) +
     labs(x="Datasize (%)", y="AUC")
ggsave(paste0(AGP_res_file, "AGP_AUC_diff_size.png"), p,
       width = 20, height = 10, units = "cm")


# cross IBD
IDB_res_file <- "/home/dongbiao/word_embedding_microbiome/programe_test/IBD/results/"
cross_IBD_diff_size <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/IBD/attention_AGP_diff_size.out", c(12, 18))
colnames(cross_IBD_diff_size) <- c("AUC", "F1")
cross_IBD_diff_size <- as.data.frame(lapply(cross_IBD_diff_size, as.numeric))
cross_IBD_diff_size$Study <- rep(c('qiita_1629', 'qiita_2538', 'PRJNA422193', 'PRJNA431126'), 5, each=1)
cross_IBD_diff_size$datasize <- rep(c("10", "30", "50", "70", "90"), 1, each=4)

p <- cross_IBD_diff_size %>% ggplot(aes(x = datasize, y = AUC)) + geom_boxplot() + geom_point(aes(color = Study), size = 3) +
     geom_line(aes(group = Study, y = AUC, color = Study), alpha = 0.5, linewidth = 2) +
     labs(x="Datasize (%)", y="AUC") +
     scale_color_brewer(palette = "Set3")
ggsave(paste0(IDB_res_file, "AGP_AUC_diff_size.png"), p,
       width = 20, height = 10, units = "cm")

# cross CRC
CRC_res_file <- "/home/dongbiao/word_embedding_microbiome/programe_test/CRC/results/"
cross_CRC_diff_size <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/CRC/attention_AGP_diff_size.out", c(12, 18))
colnames(cross_CRC_diff_size) <- c("AUC", "F1")
cross_CRC_diff_size <- as.data.frame(lapply(cross_CRC_diff_size, as.numeric))
cross_CRC_diff_size$Study <- rep(c("PRJNA763872", "PRJEB36789", "PRJNA824020", "PRJDB11845", "PRJNA318004", "PRJEB6070"), 5, each=1)
cross_CRC_diff_size$datasize <- rep(c("10", "30", "50", "70", "90"), 1, each=6)

p <- cross_CRC_diff_size %>% ggplot(aes(x = datasize, y = AUC)) + geom_boxplot() + geom_point(aes(color = Study), size = 3) +
    geom_line(aes(group = Study, y = AUC, color = Study), alpha = 0.5, linewidth = 2) +
    labs(x="Datasize (%)", y="AUC") +
    scale_color_brewer(palette = "Set3")
ggsave(paste0(CRC_res_file, "AGP_AUC_diff_size.png"), p,
       width = 20, height = 10, units = "cm")

# cross diet
dietay_res_file <- "/home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/results/"
cross_diet_diff_size <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/attention_AGP_diff_size.out", c(12, 18))
colnames(cross_diet_diff_size) <- c("AUC", "F1")
cross_diet_diff_size <- as.data.frame(lapply(cross_diet_diff_size, as.numeric))
cross_diet_diff_size$Study <- rep(c('SRP067761', 'SRP120250', 'SRP219296', 'qiita_13367', 'SRP345891'), 5, each=1)
cross_diet_diff_size$datasize <- rep(c("10", "30", "50", "70", "90"), 1, each=5)

p <- cross_diet_diff_size %>% ggplot(aes(x = datasize, y = AUC)) + geom_boxplot() + geom_point(aes(color = Study), size = 3) +
    geom_line(aes(group = Study, y = AUC, color = Study), alpha = 0.5, linewidth = 2) +
    labs(x="Datasize (%)", y="AUC") +
    scale_color_brewer(palette = "Set3")
ggsave(paste0(dietay_res_file , "AGP_AUC_diff_size.png"), p,
       width = 20, height = 10, units = "cm")

IBD_table %>% ggplot(aes(x = Group1, y = AUC_mean / 5)) +
    labs(x="", y="AUC") +
    geom_col(fill = "skyblue", alpha=0.5) +
    geom_jitter(aes(x = Group1, y = AUC_mean), width = 0.1, height = 0, size = 2) +
    theme_bw(base_size = 14) +
    facet_grid(. ~ Group2, scales = "free_x") +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold")) + # 加粗坐标轴文本
    theme(axis.text.x = element_text(angle = 65, hjust = 1, vjust = 1))

