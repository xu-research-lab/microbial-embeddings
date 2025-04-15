library(tidyverse)
library(pheatmap)

KO_predict <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/function_auc_top_IBD.csv") 
KO_predict_ratio <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/KO_ratio_IBD.csv") 
res_1 <- merge(KO_predict, KO_predict_ratio, by="subsystem") %>% filter(auc > 0.8) %>% filter(com_ratios > 2)
res_1 <- res_1[grep("spo", res_1$descript.x), ]

## Taxonomy
tax_table <- read.table("/beegfs/db/greengenes/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt", sep = "\t") %>% column_to_rownames("V1")
## KO
KO <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/pathway/picrust_filter/KO_predicted_description.tsv",
                 sep = "\t", header = 1, check.names = FALSE) %>% column_to_rownames("KO")

KO_sporulation <- KO %>% filter(description %in% res_1$descript.x)
tsne_ibd <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/t_sne_IBD.csv") %>%
    filter(emb_type == "Co_embedding")
tsne_ibd$fid <- as.character(tsne_ibd$fid)
tsne_ibd$tax <- tax_table[tsne_ibd$fid, ]
heatmap_plot_data <- KO_sporulation[, tsne_ibd$fid] %>% t() %>% as.data.frame()
row_annotation <- data.frame(group = tsne_ibd$group)
rownames(row_annotation) <- tsne_ibd$fid
p_1 <- pheatmap(heatmap_plot_data, annotation_row = row_annotation, labels_row = tsne_ibd$tax, fontsize = 8, main="IBD")
ggsave("/home/dongbiao/word_embedding_microbiome/result/supplement_spo_IBD.png", p_1,
       width = 40, height = 40, units = "cm")

tsne_crc <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/t_sne_CRC.csv") %>%
    filter(emb_type == "Co_embedding")
tsne_crc$fid <- as.character(tsne_crc$fid)
tsne_crc$tax <- tax_table[tsne_crc$fid, ]
heatmap_plot_data <- KO_sporulation[, tsne_crc$fid] %>% t() %>% as.data.frame()
row_annotation <- data.frame(group = tsne_crc$group)
rownames(row_annotation) <- tsne_crc$fid
p_2 <- pheatmap(heatmap_plot_data, annotation_row = row_annotation, labels_row = tsne_crc$tax, fontsize = 8, main="CRC")
ggsave("/home/dongbiao/word_embedding_microbiome/result/supplement_spo_CRC.png", p_2,
       width = 40, height = 40, units = "cm")
