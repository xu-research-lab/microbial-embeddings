library(dplyr)
library(ggplot2)
library(lemon)
library(cowplot)
library(stringr)
library(cowplot)
get_auc_f1 <- function(file, idx){
    table <- read.table(file, sep = "\t")
    table <- table[apply(table, 1, function(x) grepl("mcc", x)), ]
    table <- gsub(",", "", table)
    table <- str_split(table, " ", simplify = TRUE)
    table <- as.data.frame(table[, idx])
}

lou_auc <- function(plot_data, AUC_or_F1, vector_x){
    colnames(plot_data) <- c("value", "Study", "Group2")
    p <- plot_data %>% ggplot(aes(x = Group2, y = value, group=Group2)) +
        labs(x="", y=paste0(AUC_or_F1)) +
        geom_boxplot() +
        geom_point(aes(color=Study), size = 4, alpha = 0.8) +
        geom_line(aes(group = Study, y = value, color = Study)) +
        theme_classic(base_size = 14) +
        scale_color_brewer(palette = "Set3") +
        theme(text = element_text(face = "bold"), # 设置所有文本为加粗
              plot.title = element_text(hjust = 0.5), # 居中对齐标题
              axis.title = element_text(face = "bold"), # 加粗坐标轴标题
              axis.text = element_text(face = "bold")) + # 加粗坐标轴文本
        # coord_cartesian(ylim = vector_x) +
        theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))
    
    return(p)
}

group1 <- c('PRJNA422193', 'PRJNA431126', 'qiita_1629', 'qiita_2538', 'HMP',
            'PRJEB13679', 'PRJNA317429')
group2 <- c("origin", "cutmix_inter", "cutmix_intro", "mixup_inter", "mixup_intro", "cutmix_inter_phylogeny_co", 
            "cutmix_intro_phylogeny_co",  "mixup_inter_phylogeny_co", "mixup_intro_phylogeny_co")
IBD <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/IBD/attention_aug.out", c(12, 18))
colnames(IBD) <- c("AUC", "F1")
IBD$group1 <- rep(group1, 9, each = 1)
IBD$group2 <- rep(group2, 1, each = 7)
IBD$group2 <- factor(IBD$group2, levels = group2)
IBD <- IBD %>% filter(group2 %in% c("origin", "cutmix_intro", "cutmix_intro_phylogeny_co"))
p1<-lou_auc(IBD[,c("AUC", "group1", "group2")], "AUC", vector_x=c(0.4, 0.9))
lou_auc(IBD[,c("F1", "group1", "group2")], "F1", vector_x=c(0.4, 0.9))

group1 <- c('PRJEB36789', 'PRJNA824020', 'PRJDB11845', 'PRJNA318004',
            'PRJEB6070', 'PRJNA430990', 'PRJNA290926')
group2 <- c("origin", "cutmix_inter", "cutmix_intro", "mixup_inter", "mixup_intro", "cutmix_inter_phylogeny_co", 
            "cutmix_intro_phylogeny_co",  "mixup_inter_phylogeny_co", "mixup_intro_phylogeny_co")
CRC <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/CRC/attention_aug.out", c(12, 18))
colnames(CRC) <- c("AUC", "F1")
CRC$group1 <- rep(group1, 9, each = 1)
CRC$group2 <- rep(group2, 1, each = 7)
CRC$group2 <- factor(CRC$group2, levels = group2)
CRC <- CRC %>% filter(group2 %in% c("origin", "cutmix_intro", "cutmix_intro_phylogeny_co"))
p2<-lou_auc(CRC[,c("AUC", "group1", "group2")], "AUC", vector_x=c(0.4, 0.9))
lou_auc(CRC[,c("F1", "group1", "group2")], "F1", vector_x=c(0.4, 0.9))

group1 <- c('qiita_13367', 'SRP067761', 'SRP120250', 'SRP128128', 'SRP345891', 'SRP219296')
group2 <- c("origin", "cutmix_inter", "cutmix_intro", "mixup_inter", "mixup_intro", "cutmix_inter_phylogeny_co", 
            "cutmix_intro_phylogeny_co",  "mixup_inter_phylogeny_co", "mixup_intro_phylogeny_co")
fiber <- get_auc_f1("/home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/attention_aug.out", c(12, 18))
colnames(fiber) <- c("AUC", "F1")
fiber$group1 <- rep(group1, 9, each = 1)
fiber$group2 <- rep(group2, 1, each = 6)
fiber$group2 <- factor(fiber$group2, levels = group2)
fiber <- fiber %>% filter(group2 %in% c("origin", "cutmix_intro", "cutmix_intro_phylogeny_co"))
p3<-lou_auc(fiber[,c("AUC", "group1", "group2")], "AUC", vector_x=c(0.4, 0.9))

p <- plot_grid(p1, p2, p3,
               labels = c('a:IBD',
                          'b:CRC',
                          'c:Fiber'),
               nrow = 3, ncol=1, plot=FALSE)
ggsave("/home/dongbiao/word_embedding_microbiome/result/aug.png", p,
       width = 20, height = 25, units = "cm")
