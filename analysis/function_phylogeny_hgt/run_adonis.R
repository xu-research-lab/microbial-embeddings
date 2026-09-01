library(vegan)
library(proxy)  
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
func_id <- args[1]

co_embedding <- read.csv("../../data/social_niche_embedding_100.txt", row.names = 1, sep=" ", header = FALSE)
func_file <- "data/bac_KO_predicted.tsv"
func_table <- read.csv(func_file, sep = "\t", row.names = 1, check.names = FALSE)
                       
inter_id <- intersect(rownames(co_embedding), rownames(func_table))
co_embedding <- co_embedding[inter_id, , drop = FALSE]
func_table <- func_table[inter_id, , drop = FALSE]
func_table[func_table > 0] <- 1

metadata <- data.frame(group = func_table[inter_id, func_id])
rownames(metadata) <- inter_id
distance <- proxy::dist(co_embedding, method = "cosine")
res <- adonis2(distance ~ group, data = metadata, by = "terms")
saveRDS(res, file = paste0("data/adonis/adonis2_", func_id,".rds"))


