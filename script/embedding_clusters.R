library(tidyverse)
library(ggpubr) 
library(ggsci)


pointbiserial_res <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/cluster/pointbiserial_res.csv")
pointbiserial_res$PointBiserial_Correlation <- abs(pointbiserial_res$PointBiserial_Correlation)
compire <- list(c('embed','sparcc'))
pointbiserial_res %>% filter(disease == "IBD") %>%
    ggplot(aes(x=type, y=abs(PointBiserial_Correlation))) +
    geom_boxplot() + 
    facet_wrap(.~study, scales = "free", nrow = 1) +
    theme_bw(base_size = 14) + 
    geom_signif(comparisons = compire, step_increase=0.5, map_signif_level = T, 
                test = t.test, margin_top = 0.1, vjust = 1.5) +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))


pointbiserial_res %>% filter(disease == "CRC") %>%
    ggplot(aes(x=type, y=abs(PointBiserial_Correlation))) +
    geom_boxplot() + 
    facet_wrap(.~study, scales = "free", nrow = 1) +
    theme_bw(base_size = 14) + 
    geom_signif(comparisons = compire, step_increase=0.5, map_signif_level = T, 
                test = t.test, margin_top = 0.1, vjust = 1.5) +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))

clusters_abc_labels <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/cluster/smetana/clusters_abc_labels_IBD.csv")
clusters_abc_labels$group_label <- as.factor(clusters_abc_labels$group_label)
clusters_abc_labels$group_label <- as.character(clusters_abc_labels$group_label)
clusters_abc_labels$group_label[clusters_abc_labels$group_label == 0] = "Control"
clusters_abc_labels$group_label[clusters_abc_labels$group_label == 1] = "Disease"

study_name=c('PRJNA422193', 'PRJNA431126', 'qiita_1629', 'qiita_2538', 'HMP', 'PRJEB13679', 'PRJNA317429')
clusters_abc_labels$study_group[clusters_abc_labels$study_group == "study_1"] = study_name[1]
clusters_abc_labels$study_group[clusters_abc_labels$study_group == "study_2"] = study_name[2]
clusters_abc_labels$study_group[clusters_abc_labels$study_group == "study_3"] = study_name[3]
clusters_abc_labels$study_group[clusters_abc_labels$study_group == "study_4"] = study_name[4]
clusters_abc_labels$study_group[clusters_abc_labels$study_group == "study_5"] = study_name[5]
clusters_abc_labels$study_group[clusters_abc_labels$study_group == "study_6"] = study_name[6]
clusters_abc_labels$study_group[clusters_abc_labels$study_group == "study_7"] = study_name[7]
compire <- list(c('Control','Disease'))
p1 <- clusters_abc_labels %>% ggplot(aes(x = group_label, y = clusters_abc, group = group_label)) + 
    geom_jitter() + 
    geom_boxplot() + 
    theme_bw(base_size = 14) + labs(x="", y="Clusters_93_abundance", title="IBD") +
    geom_signif(comparisons = compire, step_increase=0.5, map_signif_level = T, 
                test = wilcox.test, margin_top = 0.1, vjust = 1.5) +
    facet_wrap(.~study_group, scales = "free") +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))

clusters_abc_labels <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/cluster/smetana/clusters_abc_labels_CRC.csv")
clusters_abc_labels$group_label <- as.factor(clusters_abc_labels$group_label)
clusters_abc_labels$group_label <- as.character(clusters_abc_labels$group_label)
clusters_abc_labels$group_label[clusters_abc_labels$group_label == 0] = "Control"
clusters_abc_labels$group_label[clusters_abc_labels$group_label == 1] = "Disease"
study_name=c('PRJEB36789', 'PRJNA824020', 'PRJDB11845', 'PRJNA318004', 'PRJEB6070', 'PRJNA430990', 'PRJNA290926')
clusters_abc_labels$study_group[clusters_abc_labels$study_group == "study_1"] = study_name[1]
clusters_abc_labels$study_group[clusters_abc_labels$study_group == "study_2"] = study_name[2]
clusters_abc_labels$study_group[clusters_abc_labels$study_group == "study_3"] = study_name[3]
clusters_abc_labels$study_group[clusters_abc_labels$study_group == "study_4"] = study_name[4]
clusters_abc_labels$study_group[clusters_abc_labels$study_group == "study_5"] = study_name[5]
clusters_abc_labels$study_group[clusters_abc_labels$study_group == "study_6"] = study_name[6]
clusters_abc_labels$study_group[clusters_abc_labels$study_group == "study_7"] = study_name[7]
compire <- list(c('Control','Disease'))
p2 <- clusters_abc_labels %>% ggplot(aes(x = group_label, y = log10(clusters_abc+0.0000001), group = group_label)) + 
    geom_boxplot() + 
    geom_jitter() + 
    theme_bw(base_size = 14) + labs(x="", y="Clusters_111_abundance", title="CRC") +
    geom_signif(comparisons = compire, step_increase=0.5, map_signif_level = T, 
                test = wilcox.test, margin_top = 0.1, vjust = 1.5) +
    facet_wrap(.~study_group, scales = "free") +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"))


p <- plot_grid(p1, p2, labels = c('a', 'b'),
               align="hv",
               nrow = 1, ncol=2, plot=FALSE, rel_heights = c(1, 1))

ggsave("/home/dongbiao/word_embedding_microbiome/result/emebddeing_module_supplement.png", p,
       width = 30, height = 15, units = "cm")
