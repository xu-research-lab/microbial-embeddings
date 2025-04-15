library(tidyverse)
library(cowplot)
library(reshape2)
library(ggpubr)
library(pROC) 

subsystem_predict_auc <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/res/subsystem_predict_auc.csv")

p1 <- subsystem_predict_auc %>% ggplot(aes(x=Group, y=AUC, fill=Group)) +   
    geom_boxplot()+
    scale_color_brewer(palette = "Paired") +
    labs(title = "", x = "", y = "AUC", fill="") +  
    theme_bw(base_size=14)+
    theme(text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 0),
          legend.position = "None",
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5)) + 
    stat_compare_means(method = "t.test", label.x = 1.5, label.y = 0.65, paired=TRUE)

co_subsystem_predict <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/res/co_functon_results.csv") %>%
    select(c("X", "AUC_ROC")) %>% arrange(desc(AUC_ROC))
co_subsystem_predict <- co_subsystem_predict[1:5,]
co_subsystem_predict$X <- c("Ribonuclease J family", "Na-translocating NADH-quinone\n oxidoreductase", "Na(+) H(+) antiporter",
                            "Utilization systems for glycans and\n polysaccharides (PULs, SUS)", "Lpt lipopolysaccharide transport\n system")
co_subsystem_predict$X <- factor(co_subsystem_predict$X, levels = c("Ribonuclease J family", "Na-translocating NADH-quinone\n oxidoreductase", "Na(+) H(+) antiporter",
                                                                    "Utilization systems for glycans and\n polysaccharides (PULs, SUS)", "Lpt lipopolysaccharide transport\n system"))
p2 <- co_subsystem_predict %>% ggplot(aes(x = AUC_ROC, y = X)) +
    geom_col(fill = "skyblue", alpha=0.5) +
    labs(title = "Subsystem prediction", x = "AUC", y = "") +  
    theme_bw(base_size=14)+
    coord_cartesian(xlim = c(0.6, 0.70)) +
    theme(text = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5))
im_f <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/res/feature_importance.csv", row.names = 1) %>%
    filter(Importance > 0.012) %>% arrange(desc(Importance))
im_f$Feature <- paste0("Property_", im_f$Feature)
im_f$Feature <- factor(im_f$Feature , levels = im_f$Feature)
p3 <- im_f %>% ggplot(aes(x = Importance, y = Feature)) +
    geom_col(fill = "skyblue", alpha=0.5) +
    labs(title = "Utilization systems for glycans and\n polysaccharides (PULs, SUS)", x = "Importance", y = "") +  
    theme_bw(base_size=14)+
    coord_cartesian(xlim = c(0.012, 0.02)) +
    theme(text = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 8))

rf_pro_lable <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/res/res_feature_importance.csv")
names <- unique(rf_pro_lable$group)
res_auc_plot <- list()
n <- 1
auc_res <- c()
new_names <- c("Original", "Imp. properties", "Unimp. properties")
for (i in names){
    # 计算 ROC 曲线  
    temp <- rf_pro_lable %>% filter(group == i)
    roc_result <- roc(temp$label, temp$probs)
    # 绘制 ROC 曲线  
    roc_data <- data.frame(  
        Sensitivity = rev(roc_result$sensitivities),  # 真实阳性率  
        Specificity = rev(roc_result$specificities),  # 真实阴性率  
        Thresholds = roc_result$thresholds  # 阈值  
    )  
    roc_data$group <- rep(new_names[n], n = nrow(roc_data))
    res_auc_plot[[n]] = roc_data
    auc_res <- c(auc_res, round(auc(roc(temp$label, temp$probs)), 3))
    n <- n + 1
}
res_auc_plot <- rbind(res_auc_plot[[1]], res_auc_plot[[2]], res_auc_plot[[3]])
res_auc_plot$group <- factor(res_auc_plot$group, levels = c("Original", "Imp. properties", "Unimp. properties"))
p4 <- ggplot(res_auc_plot, aes(x = 1 - Specificity, y = Sensitivity, color=group)) +  
    geom_line(size = 1) +           # 画曲线  
    scale_color_manual(values=c("Original"="#8ECFC9", "Imp. properties"="#FFBE7A", "Unimp. properties"="#82B0D2")) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + # 参考线  
    labs(title = "Utilization systems for glycans and\n polysaccharides (PULs, SUS)", color="",  
         x = "1 - Specificity",  
         y = "Sensitivity") +  
    annotate("text", x = 0.75, y = 0.4, label = paste("AUC:", auc_res[1]), size = 4, color="#8ECFC9") +
    annotate("text", x = 0.75, y = 0.32, label = paste("AUC:", auc_res[2]), size = 4, color="#FFBE7A") +
    annotate("text", x = 0.75, y = 0.24, label = paste("AUC:", auc_res[3]), size = 4, color="#82B0D2") +
    theme_bw(base_size=14)+
    theme(text = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 8))

get_auc_f1 <- function(file, idx){
    table <- read.table(file, sep = "\t")
    table <- table[apply(table, 1, function(x) grepl("mcc", x)), ]
    table <- gsub(",", "", table)
    table <- str_split(table, " ", simplify = TRUE)
    table <- as.data.frame(table[, idx])
}


Study <- c('PRJNA422193', 'PRJNA431126', 'qiita_1629', 'qiita_2538', 'HMP', 'PRJEB13679', 'PRJNA317429')
attention_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/classification/IBD/attention.out"
table <- get_auc_f1(attention_res, c(12, 18))
colnames(table) <- c("AUC", "F1")
table$Study <- rep(Study, 3, each=1)
table$AUC <- as.numeric(table$AUC)
table$Group <- rep(c("Original", "Imp. properties", "Unimp. properties"), 1, each= 7)
table$Group <- factor(table$Group, levels = c("Original", "Imp. properties", "Unimp. properties"))
p5 <- table %>% ggplot(aes(x = Group, y = AUC, group = Group)) +
    geom_boxplot() + geom_line(aes(group = Study, color = Study)) +
    geom_point(aes(color = Study), size = 3) +
    theme_bw(base_size = 14) +
    labs(x="", y="AUC", title = "IBD") +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text.x = element_text(face = "bold", angle = 15, hjust = 1, vjust = 1)) +
    scale_color_brewer(palette = "Paired")


Study <- c('PRJEB36789', 'PRJNA824020', 'PRJDB11845', 'PRJNA318004', 'PRJEB6070', 'PRJNA430990', 'PRJNA290926')
attention_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/classification/CRC/attention.out"
table <- get_auc_f1(attention_res, c(12, 18))
colnames(table) <- c("AUC", "F1")
table$Study <- rep(Study, 3, each=1)
table$AUC <- as.numeric(table$AUC)
table$Group <- rep(c("Original", "Imp. properties", "Unimp. properties"), 1, each= 7)
table$Group <- factor(table$Group, levels = c("Original", "Imp. properties", "Unimp. properties"))
p6 <- table %>% ggplot(aes(x = Group, y = AUC, group = Group)) +
    geom_boxplot() + geom_line(aes(group = Study, color = Study)) +
    geom_point(aes(color = Study), size = 3) +
    theme_bw(base_size = 14) +
    labs(x="", y="AUC", title = "CRC") +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text.x = element_text(face = "bold", angle = 15, hjust = 1, vjust = 1)) +
    scale_color_brewer(palette = "Paired")

upp_plot <- plot_grid(p1, p2,
                      labels = c('a', 'b'),
                      align="hv",
                      scale = c(1, 1),
                      nrow = 1, ncol=2, plot=FALSE, rel_widths = c(1, 1))

bottom_plot <- plot_grid(p3, p4, p5, p6,
                         labels = c('c', 'd', 'e', 'f'),
                         align="hv",
                         scale = c(1, 1, 1, 1),
                         nrow = 1, ncol=4, plot=FALSE, rel_widths = c(1, 1, 1, 1))
p <- plot_grid(p1, p2,p3, p4, p5, p6,
               labels = c('a', 'b', 'c', 'd', 'e', 'f'),
               align="hv",
               scale = c(1, 1, 1, 1, 1, 1),
               nrow = 3, ncol=2, plot=FALSE, rel_heights = c(0.8, 1, 1, 1))

ggsave("/home/dongbiao/word_embedding_microbiome/result/embedding_function_explain.png", p,
       width = 28, height = 28, units = "cm")
