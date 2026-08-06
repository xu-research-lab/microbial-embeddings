library(LinDA)
library(dplyr)
library(metafor)
library(tidyverse)
library(biomformat)
input_biom <- biomformat::read_biom("data/loo_all_diseases/total_plus_IBD_final.biom")
metadata <- read.table("data/loo_all_diseases/total_plus_IBD_final_metadata.tsv", sep = "\t", header = TRUE) %>%
    column_to_rownames("sample")

otu_data <- do.call(data.frame, input_biom$data)
mat <- data.frame(do.call(rbind, input_biom$rows))
colnames(otu_data) <- mat$id

# pool data caculate DA
linda_res <- linda(t(otu_data), metadata[rownames(otu_data),], formula = '~group + (1 | study)', alpha = 0.05,
                   prev.cut = 0.001, lib.cut = 1000, winsor.quan = 0.97)

res <- linda_res$output$group
write.csv(res, "data/loo_all_diseases/linda_res.csv")