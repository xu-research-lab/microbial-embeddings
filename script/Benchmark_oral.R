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
        theme(axis.text.x = element_text(angle = 20, hjust = 0.7, vjust = 0.8))
    
    return(p)
}


# OSCC
oscc_table_1 <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/OSCC/attention.out", c(12, 18))
oscc_table_2 <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/OSCC/mlp.out", c(12, 18))
oscc_table_3 <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/OSCC/randownforest.out", c(4, 12))
oscc_table <- rbind(oscc_table_1, oscc_table_2, oscc_table_3)
oscc_table <- as.data.frame(lapply(oscc_table, as.numeric))
colnames(oscc_table) <- c("AUC", "F1")
lable_1 <- c('PRJNA756784', 'PRJNA751046', 'PRJEB39064', 'PRJNA700849', 
             'OEP000837', 'PRJNA421234', 'PRJNA386665', 'PRJNA744870')
lable_2 <- c("Abundance-percentile", "RandInit")
lable_3 <- c("Abundance-percentile", "OTU")
oscc_table$Group1 <- c(rep(lable_2, 1, each=8), rep(lable_3, 2, each=8))
oscc_table$Group2 <- rep(c("Attention", "MLP", "RF"), 1, each=16)
oscc_table$study <- rep(lable_1, 3)
oscc_table_mean_sd <- oscc_table %>% group_by(Group1, Group2) %>% summarise(AUC_mean = mean(AUC), AUC_sd = sd(AUC),
                                                                          F1_mean = mean(F1), F1_sd = sd(F1))
oscc_table <- merge(oscc_table, oscc_table_mean_sd, by = c("Group1", "Group2"))

# caries
caries_table_1 <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/caries/attention.out", c(12, 18))
caries_table_2 <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/caries/mlp.out", c(12, 18))
caries_table_3 <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/caries/randownforest.out", c(4, 12))
caries_table <- rbind(caries_table_1, caries_table_2, caries_table_3)
caries_table <- as.data.frame(lapply(caries_table, as.numeric))
colnames(caries_table) <- c("AUC", "F1")
lable_1 <- c('PRJNA325084', 'PRJNA339212', 'PRJNA383868', 'PRJNA454811',
             'PRJNA480252', 'PRJNA495719', 'PRJNA681486')
lable_2 <- c("Abundance-percentile", "RandInit")
lable_3 <- c("Abundance-percentile", "OTU")
caries_table$Group1 <- c(rep(lable_2, 2, each=7), rep(lable_3, 1, each=7))
caries_table$Group2 <- rep(c("Attention", "MLP", "RF"), 1, each=14)
caries_table$study <- rep(lable_1, 3)
caries_table_mean_sd <- caries_table %>% group_by(Group1, Group2) %>% summarise(AUC_mean = mean(AUC), AUC_sd = sd(AUC),
                                                                            F1_mean = mean(F1), F1_sd = sd(F1))
caries_table <- merge(caries_table, caries_table_mean_sd, by = c("Group1", "Group2"))

# paradentities
paradentities_table_1 <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/paradentities/attention.out", c(12, 18))
paradentities_table_2 <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/paradentities/mlp.out", c(12, 18))
paradentities_table_3 <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/paradentities/randownforest.out", c(4, 12))
paradentities_table <- rbind(paradentities_table_1, paradentities_table_2, paradentities_table_3)
paradentities_table <- as.data.frame(lapply(paradentities_table, as.numeric))
colnames(paradentities_table) <- c("AUC", "F1")
lable_1 <- c('PRJEB42371', 'PRJEB61123', 'PRJNA321534', 'PRJNA578492',
             'PRJNA580506', 'PRJNA601054', 'PRJNA650272', 'PRJNA843376', 'SRP009299')
lable_2 <- c("Abundance-percentile", "RandInit")
lable_3 <- c("Abundance-percentile", "OTU")
paradentities_table$Group1 <- c(rep(lable_2, 1, each=9), rep(lable_3, 2, each=9))
paradentities_table$Group2 <- rep(c("Attention", "MLP", "RF"), 1, each=18)
paradentities_table$study <- rep(lable_1, 3)
paradentities_table_mean_sd <- paradentities_table %>% group_by(Group1, Group2) %>% summarise(AUC_mean = mean(AUC), AUC_sd = sd(AUC),
                                                                                F1_mean = mean(F1), F1_sd = sd(F1))
paradentities_table <- merge(paradentities_table, paradentities_table_mean_sd, by = c("Group1", "Group2"))

# plot
group2 <- c("Attention", "MLP", "RF")
plot_list <- list()
# OSCC
OSCC_res_file <- "/home/dongbiao/word_embedding_microbiome/programe_test/OSCC/results/"
for (i in 1 : length(group2)){
    plot_list[[i]] <- lou_auc(classify = group2[i],
                              plot_data = oscc_table[,c("AUC", "AUC_mean", "AUC_sd", "Group1", "Group2", "study")], 
                              AUC_or_F1 = "AUC", 
                              vector_x = c(0.4, 1),
                              nums = 8)
}
p<- grid_arrange_shared_legend(plot_list[[1]], plot_list[[2]], plot_list[[3]], 
                               nrow = 1, position='bottom')
ggsave(paste0(OSCC_res_file, "AUC_cross_data.png"), p,
       width = 18, height = 8, units = "cm")

for (i in 1 : length(group2)){
    plot_list[[i]] <- lou_auc(classify = group2[i],
                              plot_data = oscc_table[,c("F1", "F1_mean", "F1_sd", "Group1", "Group2", "study")], 
                              AUC_or_F1 = "F1", 
                              vector_x = c(0.2, 0.9),
                              nums = 8)
}
p<- grid_arrange_shared_legend(plot_list[[1]], plot_list[[2]],
                               plot_list[[3]],
                               nrow = 1, position='bottom')
ggsave(paste0(OSCC_res_file, "F1_cross_data.png"), p,
       width = 18, height = 8, units = "cm")

# caries
caries_res_file <- "/home/dongbiao/word_embedding_microbiome/programe_test/caries/results/"
for (i in 1 : length(group2)){
    plot_list[[i]] <- lou_auc(classify = group2[i],
                              plot_data = caries_table[,c("AUC", "AUC_mean", "AUC_sd", "Group1", "Group2", "study")], 
                              AUC_or_F1 = "AUC", 
                              vector_x = c(0.2, 1),
                              nums = 7)
}
p<- grid_arrange_shared_legend(plot_list[[1]], plot_list[[2]], plot_list[[3]], 
                               nrow = 1, position='bottom')
ggsave(paste0(caries_res_file, "AUC_cross_data.png"), p,
       width = 18, height = 8, units = "cm")

for (i in 1 : length(group2)){
    plot_list[[i]] <- lou_auc(classify = group2[i],
                              plot_data = caries_table[,c("F1", "F1_mean", "F1_sd", "Group1", "Group2", "study")], 
                              AUC_or_F1 = "F1", 
                              vector_x = c(0.15, 0.8),
                              nums = 7)
}
p<- grid_arrange_shared_legend(plot_list[[1]], plot_list[[2]],
                               plot_list[[3]],
                               nrow = 1, position='bottom')
ggsave(paste0(caries_res_file, "F1_cross_data.png"), p,
       width = 18, height = 8, units = "cm")

# paradentities
paradentities_res_file <- "/home/dongbiao/word_embedding_microbiome/programe_test/paradentities/results/"
for (i in 1 : length(group2)){
    plot_list[[i]] <- lou_auc(classify = group2[i],
                              plot_data = paradentities_table[,c("AUC", "AUC_mean", "AUC_sd", "Group1", "Group2", "study")], 
                              AUC_or_F1 = "AUC", 
                              vector_x = c(0.35, 1),
                              nums = 9)
}
p<- grid_arrange_shared_legend(plot_list[[1]], plot_list[[2]], plot_list[[3]], 
                               nrow = 1, position='bottom')
ggsave(paste0(paradentities_res_file, "AUC_cross_data.png"), p,
       width = 18, height = 8, units = "cm")

for (i in 1 : length(group2)){
    plot_list[[i]] <- lou_auc(classify = group2[i],
                              plot_data = paradentities_table[,c("F1", "F1_mean", "F1_sd", "Group1", "Group2", "study")], 
                              AUC_or_F1 = "F1", 
                              vector_x = c(0.2, 0.95),
                              nums = 9)
}
p<- grid_arrange_shared_legend(plot_list[[1]], plot_list[[2]],
                               plot_list[[3]],
                               nrow = 1, position='bottom')
ggsave(paste0(paradentities_res_file, "F1_cross_data.png"), p,
       width = 18, height = 8, units = "cm")

# cross OSCC
oscc_res_file <- "/home/dongbiao/word_embedding_microbiome/programe_test/filter_prevalence/diff_datasize/oscc_cross/"
cross_oscc_diff_size <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/filter_prevalence/sub_embedding/oscc_cross/attention.out", c(12, 18))
colnames(cross_oscc_diff_size) <- c("AUC", "F1")
cross_oscc_diff_size <- as.data.frame(lapply(cross_oscc_diff_size, as.numeric))
cross_oscc_diff_size$Study <- rep(c('PRJNA756784', 'PRJNA751046', 'PRJEB39064', 'PRJNA700849', 
                                    'OEP000837', 'PRJNA421234', 'PRJNA386665', 'PRJNA744870'), 4, each=3)
cross_oscc_diff_size$datasize <- rep(c("5000", "10000" ,"15000", "20000"), 1, each=24)
line_plot <- cross_oscc_diff_size %>% group_by(Study, datasize) %>% summarise_all(mean)
cross_oscc_diff_size <- merge(cross_oscc_diff_size, line_plot, by=c("Study", "datasize"))
cross_oscc_diff_size$datasize <- factor(cross_oscc_diff_size$datasize, levels=c("5000", "10000", "15000", "20000"))

p <- cross_oscc_diff_size %>% ggplot(aes(x = datasize, y = AUC.x)) + geom_boxplot() + geom_point(aes(color = Study), size = 2) +
    geom_line(aes(group = Study, y = AUC.y, color = Study), alpha = 0.5, linewidth = 1) +
    labs(x="Datasize", y="AUC") +
    theme_bw(base_size = 10) +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold")) + # 加粗坐标轴文本
    scale_color_brewer(palette = "Set3")
ggsave(paste0(oscc_res_file , "AUC_diff_size.png"), p,
       width = 12, height = 6, units = "cm")

p <- cross_oscc_diff_size %>% ggplot(aes(x = datasize, y = F1.x)) + geom_boxplot() + geom_point(aes(color = Study), size = 2) +
    geom_line(aes(group = Study, y = F1.y, color = Study), alpha = 0.5, linewidth = 1) +
    labs(x="Datasize", y="F1") +
    theme_bw(base_size = 10) +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold")) +# 加粗坐标轴文本
    scale_color_brewer(palette = "Set3")
ggsave(paste0(dietay_res_file , "F1_diff_size.png"), p,
       width = 12, height = 6, units = "cm")

# cross caries
caries_res_file <- "/home/dongbiao/word_embedding_microbiome/programe_test/filter_prevalence/diff_datasize/caries_cross/"
cross_caries_diff_size <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/filter_prevalence/sub_embedding/caries_cross/attention.out", c(12, 18))
colnames(cross_caries_diff_size) <- c("AUC", "F1")
cross_caries_diff_size <- as.data.frame(lapply(cross_caries_diff_size, as.numeric))
cross_caries_diff_size$Study <- rep(c('PRJNA325084', 'PRJNA339212', 'PRJNA383868', 'PRJNA454811',
                                      'PRJNA480252', 'PRJNA495719', 'PRJNA681486'), 4, each=3)
cross_caries_diff_size$datasize <- rep(c("5000", "10000" ,"15000", "20000"), 1, each=21)
line_plot <- cross_caries_diff_size %>% group_by(Study, datasize) %>% summarise_all(mean)
cross_caries_diff_size <- merge(cross_caries_diff_size, line_plot, by=c("Study", "datasize"))
cross_caries_diff_size$datasize <- factor(cross_caries_diff_size$datasize, levels=c("5000", "10000", "15000", "20000"))

p <- cross_caries_diff_size %>% ggplot(aes(x = datasize, y = AUC.x)) + geom_boxplot() + geom_point(aes(color = Study), size = 2) +
    geom_line(aes(group = Study, y = AUC.y, color = Study), alpha = 0.5, linewidth = 1) +
    labs(x="Datasize", y="AUC") +
    theme_bw(base_size = 10) +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold")) + # 加粗坐标轴文本
    scale_color_brewer(palette = "Set3")
ggsave(paste0(caries_res_file , "AUC_diff_size.png"), p,
       width = 12, height = 6, units = "cm")

p <- cross_caries_diff_size %>% ggplot(aes(x = datasize, y = F1.x)) + geom_boxplot() + geom_point(aes(color = Study), size = 2) +
    geom_line(aes(group = Study, y = F1.y, color = Study), alpha = 0.5, linewidth = 1) +
    labs(x="Datasize", y="F1") +
    theme_bw(base_size = 10) +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold")) +# 加粗坐标轴文本
    scale_color_brewer(palette = "Set3")
ggsave(paste0(caries_res_file , "F1_diff_size.png"), p,
       width = 12, height = 6, units = "cm")

# cross paradentities
paradentities_res_file <- "/home/dongbiao/word_embedding_microbiome/programe_test/filter_prevalence/diff_datasize/periodontitis_cross/"
cross_paradentities_diff_size <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/filter_prevalence/sub_embedding/paradentities_cross/attention.out", c(12, 18))
colnames(cross_paradentities_diff_size) <- c("AUC", "F1")
cross_paradentities_diff_size <- as.data.frame(lapply(cross_paradentities_diff_size, as.numeric))
cross_paradentities_diff_size$Study <- rep(c('PRJEB42371', 'PRJEB61123', 'PRJNA321534', 'PRJNA578492',
                                      'PRJNA580506', 'PRJNA601054', 'PRJNA650272', 'PRJNA843376'), 
                                    4, each=3)
cross_paradentities_diff_size$datasize <- rep(c("5000", "10000" ,"15000", "20000"), 1, each=24)
line_plot <- cross_paradentities_diff_size %>% group_by(Study, datasize) %>% summarise_all(mean)
cross_paradentities_diff_size <- merge(cross_paradentities_diff_size, line_plot, by=c("Study", "datasize"))
cross_paradentities_diff_size$datasize <- factor(cross_paradentities_diff_size$datasize, levels=c("5000", "10000", "15000", "20000"))

p <- cross_paradentities_diff_size %>% ggplot(aes(x = datasize, y = AUC.x)) + geom_boxplot() + geom_point(aes(color = Study), size = 2) +
    geom_line(aes(group = Study, y = AUC.y, color = Study), alpha = 0.5, linewidth = 1) +
    labs(x="Datasize", y="AUC") +
    theme_bw(base_size = 10) +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold")) + # 加粗坐标轴文本
    scale_color_brewer(palette = "Set3")
ggsave(paste0(paradentities_res_file , "AUC_diff_size.png"), p,
       width = 12, height = 6, units = "cm")

p <- cross_caries_diff_size %>% ggplot(aes(x = datasize, y = F1.x)) + geom_boxplot() + geom_point(aes(color = Study), size = 2) +
    geom_line(aes(group = Study, y = F1.y, color = Study), alpha = 0.5, linewidth = 1) +
    labs(x="Datasize", y="F1") +
    theme_bw(base_size = 10) +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold")) +# 加粗坐标轴文本
    scale_color_brewer(palette = "Set3")
ggsave(paste0(paradentities_res_file, "F1_diff_size.png"), p,
       width = 12, height = 6, units = "cm")
