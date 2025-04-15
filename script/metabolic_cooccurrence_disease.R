library(tidyverse)

taxonomy <- read.table("/beegfs/db/greengenes/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt", sep = "\t", row.names = 1)
temp <- str_split(taxonomy$V2, ";", simplify = TRUE)
tax <- data.frame(temp)
colnames(tax) <- c("K", "P", "C", "O", "F", "G", "S")
rownames(tax) <- rownames(taxonomy)
tax["493071", ]
critical_value <- function(data){
    # 计算均值和标准差
    mean_value <- mean(data)
    std_dev <- sd(data)  # 默认使用样本标准差（ddof=1）
    # 定义显著性水平
    alpha <- 0.05
    # 计算左尾和右尾的临界值
    left_critical_value <- qnorm(alpha, mean = mean_value, sd = std_dev)
    right_critical_value <- qnorm(1 - alpha, mean = mean_value, sd = std_dev)
    return(right_critical_value)
}

cosine_sample <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/metabolic_cooccurrence/data/fid_pair_cosine_sample.csv")
right_critical_value <- round(critical_value(cosine_sample$cor), 2)
p1 <- cosine_sample %>% ggplot(aes(x = cor)) +
    geom_density(fill = "skyblue", alpha = 1) +
    geom_vline(xintercept = right_critical_value, color = "red", linetype = "dashed", size = 1) +
    # annotate("text", x = right_critical_value, y = max(df$y), label = "x = 5", vjust = -1, color = "red") + 
    theme_bw() +
    labs(x="Embedding cosine", 
         y="Density", title="Random OTU pairs")

smetana_sample <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/metabolic_cooccurrence/data/smetana_random.csv")
right_critical_value <- round(critical_value(smetana_sample$semtana), 2)
p2 <- smetana_sample %>% ggplot(aes(x = semtana)) +
    geom_density(fill = "skyblue", alpha = 1) +
    geom_vline(xintercept = right_critical_value, color = "red", linetype = "dashed", size = 1) +
    # annotate("text", x = right_critical_value, y = max(df$y), label = "x = 5", vjust = -1, color = "red") + 
    theme_bw() +
    labs(x="Smetana values", 
         y="Density", title="Random OTU pairs")

ibd_cooccurrence <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/metabolic_cooccurrence/data/ibd_cooccurrence.csv", 
                             row.names = 1)
metadata <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/IBD/data/metadata.txt", sep="\t", row.names = 1)
ibd_cooccurrence$study <- metadata[rownames(ibd_cooccurrence),]$study
tax[c("4320143", "363514"), ]
p3 <- ibd_cooccurrence %>% ggplot(aes(x=X4320143, y=X363514)) + 
    geom_point()+
    theme_bw()+
    labs()+
    labs(x="f__Ruminococcaceae_4320143", y="f__Ruminococcaceae_363514", title="IBD:Smetana: 5.26, Cosine:0.56")+
    facet_wrap(.~study, scales = "free", nrow=1)

crc_cooccurrence <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/metabolic_cooccurrence/data/crc_cooccurrence.csv", 
                             row.names = 1)
metadata <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/CRC/data/metadata.txt", sep="\t", row.names = 1)
study_name=c('PRJEB36789', 'PRJNA824020', 'PRJDB11845', 'PRJNA318004',
             'PRJEB6070', 'PRJNA430990', 'PRJNA290926')
metadata <- metadata %>% filter(BioProject %in% study_name)
crc_cooccurrence <- crc_cooccurrence[rownames(metadata),]
crc_cooccurrence$study <- metadata[rownames(crc_cooccurrence),]$BioProject
tax[c("366906", "335550"), ]
p4 <- crc_cooccurrence %>% ggplot(aes(x=X366906, y=X335550)) + 
    geom_point()+
    theme_bw()+
    labs(x="f__Lachnospiraceae_366906", y="f__Ruminococcaceae;g__Oscillospira_335550", title="CRC:Smetana: 5.10, Cosine:0.9")+
    facet_wrap(.~study, scales = "free", nrow=1)

upp_plot <- plot_grid(p1, p2,
                      labels = c('a', 'b'),
                      align="hv",
                      nrow = 1, ncol=2, plot=FALSE, rel_widths = c(1, 1, 1))

p <- plot_grid(upp_plot, p3, p4,
               align="hv",labels = c('', 'c', 'd'),
               nrow = 3, ncol=1, plot=FALSE, rel_heights = c(1, 1))

ggsave("/home/dongbiao/word_embedding_microbiome/result/metabolic_cooccurrence_disease.png", p,
       width = 40, height = 30, units = "cm")
