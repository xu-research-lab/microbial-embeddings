#!/bin/Rscript

library(phylolm)
library(ape)
library(tidyverse)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
func_id <- as.numeric(args[1])
pruned_tree <- read.tree("Data/pruned_tree.tre")


func_file <- "Data/bac_KO_predicted.tsv"
func_table <- read.csv(func_file, row.names = 1, sep="\t")
fid <- pruned_tree$tip.label
func_table <- func_table[fid, ]
func_table <- func_table[, colSums(func_table != 0) > 0]
non_zero_counts <- colSums(func_table != 0)
sample_threshold <- 0.8 * nrow(func_table)
func_table <- func_table[, non_zero_counts >= 100 & non_zero_counts <= sample_threshold]
func_table <- func_table %>%
    mutate(across(where(is.numeric), ~ ifelse(. > 0, 1, .)))
selected_tips <- rownames(func_table)

co_embedding <- read.csv("../../data/social_niche_embedding_100.txt",
                         row.names = 1, sep=" ", header = FALSE)
co_embedding <- co_embedding[!rownames(co_embedding) %in% "<unk>", ]
co_embedding <- co_embedding[selected_tips, ]

data <- cbind(func_table[, func_id], co_embedding)
colnames(data)[1] <- "Trait"
formula_str <- paste("Trait ~", paste("V", 2:101, sep = "", collapse = " + "))
formula <- as.formula(formula_str)
data <- data[match(pruned_tree$tip.label, rownames(data)), ]
data$Trait <- factor(data$Trait)
pruned_tree$node.label <- paste0("Node_", 1:length(pruned_tree$node.label))
pruned_tree$edge.length[pruned_tree$edge.length == 0] <- 1e-6
### null model
null_model <- phyloglm(
    Trait ~ 1,
    data = data,
    phy = pruned_tree,
    method = "logistic_MPLE",   # Max penalized likelihood; alternative: 'logistic_IG10' ,'logistic_MPLE'
    btol = 50,                  # Optimization tolerance
    start.beta = NULL,          # Starting values (optional)
    start.alpha = 1
)
saveRDS(null_model, file = paste0("Data/phylolm/co_null_model_",
                             colnames(func_table)[func_id],".rds"))

model <- phyloglm(
    formula,
    data = data,
    phy = pruned_tree,
    method = "logistic_MPLE",   # Max penalized likelihood; alternative: 'logistic_IG10'
    btol = 50,                  # Optimization tolerance
    start.beta = NULL,          # Starting values (optional)
    start.alpha = 1
)

saveRDS(model, file = paste0("Data/phylolm/co_model_",
                             colnames(func_table)[func_id],".rds"))
