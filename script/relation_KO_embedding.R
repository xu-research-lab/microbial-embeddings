library(plyr)
library(ggplot2)
library(dplyr)
library(purrr)
library(stringr)
library(pheatmap)
library(cowplot)
library(RColorBrewer)

embedding_file <- "/home/dongbiao/word_embedding_microbiome/databaes/16Sdatabase/glove_1/"
embedding_type <- "all_feces"

# embedding_file <- "/home/dongbiao/word_embedding_microbiome/programe_test/AGP_run_20220702/"
# embedding_type <- "age"

methods <- c("russell_rao", "russell_rao_weight", "faith", "jaccard", "abundance-totalsum", 
             "abundance-percentile", "braycurtis-totalsum", "braycurtis-percentile", "phylogeny", "PCA")

labs<- c("Russell_rao", "Russell_rao_weight", "Faith", "Jaccard", "Abundance_totalsum", "Abundance_percentile",
         "Braycurtis_totalsum", "Braycurtis_percentile", "phylogeny_glove", "Phylogeny_PCA")

# 相关比(Correlation Ratio)
caculate_cor_ratio <- function(table, embedding_table){
    # filter pathway table
    total_col <- nrow(table)
    keep <- colSums(table) < round(total_col * 0.9, 0)
    keep2 <- colSums(table) > round(total_col * 0.1, 0)
    table <- table[, keep & keep2]
    # Match, clean, and center
    taxa_names <- intersect(rownames(table), rownames(embedding_table))
    table <- table[taxa_names, ]
    embedding_table <- embedding_table[taxa_names, ]
    colnames(embedding_table) <- paste("property_100_", seq(1, 100), sep = "")
    
    # caculate the FoldChange and P-value
    # KO table是一个由0和1组成的矩阵
    sd_all <- apply(embedding_table, 2, sd)
    mean_all <- apply(embedding_table, 2, mean)
    cor_ratio <- list()
    for (idx in 1: ncol(table)){
        KO_presence <- embedding_table[table[, idx] == 1, ]
        KO_disappear <- embedding_table[table[, idx] == 0, ]
        
        cor_ratio[[idx]] <- sqrt(((colMeans(KO_presence) - mean_all) ^2 * nrow(KO_presence) + 
                                      (colMeans(KO_disappear) - mean_all) ^2 * nrow(KO_disappear)) / (nrow(embedding_table) * sd_all^2))
        
    }
    cor_ratio <- data.frame(cor_ratio)
    colnames(cor_ratio) <- colnames(table)
    return(cor_ratio)
}

# import pathway table OTU * KEGG-pathway
KO_table <- read.table("/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/pathway/picrust/KO_predicted.tsv",
                       header = TRUE, row.names = 1)
res_dir <- "/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/res/relation/"

# caculate corelation ratio
cor_list <- list()
for (i in 1: 10){
    # import embedding vector
    file_dir <- sprintf("%s/%s/%s_100.txt", embedding_file, methods[i], methods[i])
    embedding_table <- read.table(file_dir, row.names = 1, header = FALSE)
    embedding_table <- embedding_table[rownames(embedding_table) != "<unk>", ]
    
    cor_list[[i]] <- caculate_cor_ratio(KO_table, embedding_table)
}

saveRDS(cor_list, sprintf("/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/pathway/cor_list_KO_%s.rds", embedding_type))

# 获取每个代谢通路最大的相关系数（不考虑是在哪一维度，任意维度）
max_cor <- lapply(cor_list, apply, 2, max)
max_cor <- as.data.frame(max_cor)
colnames(max_cor) <- labs

# co-embedding / PCA 
ratio_data <- max_cor[, 1:9] / max_cor$Phylogeny_PCA
p <- pheatmap(ratio_data, show_rownames = FALSE)
# write.csv(ratio_data, "/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/pathway/ratio_data.csv")

ggsave(sprintf("/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/pathway/KO_embedding_corelation_%s.png", embedding_type), p, 
       width = 40, height = 25, units = 'cm', dpi = 300)

# plot the distribution of K09155 with embedding property all feces sample


