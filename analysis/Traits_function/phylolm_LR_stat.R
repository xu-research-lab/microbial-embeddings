library(tidyverse)
library(cowplot)
library(reshape2)
library(ggpubr)
library(pROC) 

predict_nsti <- read.csv("Data/bac_marker_predicted_and_nsti.tsv", sep="\t")
predict_nsti <- predict_nsti %>% filter(metadata_NSTI < 1)
predict_nsti <- predict_nsti %>%
    group_by(closest_reference_genome) %>%
    filter(metadata_NSTI == min(metadata_NSTI)) %>%
    slice(1) %>%  
    ungroup()
predict_nsti <- predict_nsti %>% arrange(metadata_NSTI)
predict_nsti <- predict_nsti[1:1000, ]

# func_file <- "/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/pathway/MetaCyc_pathway.csv"
func_file <- "Data/bac_KO_predicted.tsv"
func_table <- read.csv(func_file, row.names = 1, sep="\t")
fid <- intersect(predict_nsti$sequence, rownames(func_table))
func_table <- func_table[fid, ]
func_table <- func_table[, colSums(func_table != 0) > 0]
non_zero_counts <- colSums(func_table != 0)
sample_threshold <- 0.8 * nrow(func_table)
func_table <- func_table[, non_zero_counts >= 100 & non_zero_counts <= sample_threshold]

result_path <- "Data/phylolm"

# 初始化存储向量（推荐预分配内存以提高效率）
p_values <- vector("numeric", length = ncol(func_table))
names(p_values) <- colnames(func_table)  
LR_stat <- vector("numeric", length = ncol(func_table))
names(LR_stat) <- colnames(func_table)
converged <- vector("logical", length = ncol(func_table))  
pathway <- colnames(func_table)
for (i in seq_along(func_table)) {
    if (i %% 10 == 0) message("Processing subsystem ", i, "/", ncol(func_table))
    
    model_file <- paste0(result_path, "co_model_", pathway[i], ".rds")
    null_file <- paste0(result_path, "co_null_model_", pathway[i], ".rds")
    
    if (!all(file.exists(c(model_file, null_file)))) {
        warning(paste("Model files missing for patthway", pathway[i]))
        next
    }
    
    tryCatch({
        model <- readRDS(model_file)
        null_model <- readRDS(null_file)

        if (!all(c("logLik", "convergence") %in% names(model))) {
            stop("Invalid model structure")
        }
        if (model$convergence == 0 && null_model$convergence == 0) {
            LR_stat[i] <- 2 * (model$logLik - null_model$logLik)
            p_values[i] <- pchisq(LR_stat[i], df = 100, lower.tail = FALSE)
            converged[i] <- TRUE
        } else {
            p_values[i] <- NA_real_
            LR_stat[i] <- NA_real_
            converged[i] <- FALSE
        }
    }, error = function(e) {
        warning(paste("Error in pathway", pathway[i], ":", e$message))
        p_values[i] <- NA_real_
        LR_stat[i] <- NA_real_
        converged[i] <- FALSE
    })
}
valid_p <- p_values[!is.na(p_values) & converged]
valid_LR_stat <- LR_stat[!is.na(p_values) & converged]
convergence_rate <- mean(converged) * 100

# 输出结果（可选择保存） 
results <- list(
    p_values = p_values,
    LR_stat = valid_LR_stat,
    converged = converged,
    valid_p = valid_p
)
res <- data.frame(results$valid_p, 
                   p_adj_bh = p.adjust(results$valid_p, method = "BH"),
                   LR_stat=results$LR_stat)
# res <- res %>% filter(p_adj_bh < 0.05) %>% arrange(desc(LR_stat))
ko_names <- data.frame(strsplit(rownames(res), split = "\\.", fixed = FALSE))
rownames(res) <- ko_names[2,]
write.csv(res, "Data/KO_phyloglm.csv")
res$descp <- c("K01546 kdpA; potassium-transporting \nATPase potassium-binding subunit",
               "K01547 kdpB; potassium-transporting \nATPase ATP-binding subunit",
               "K01548 kdpC; potassium-transporting \nATPase KdpC subunit",
               "K03320 amt, AMT, MEP; ammonium \ntransporter, Amt family",
               "K07118 uncharacterized protein",
               "K03827 yjaB; putative acetyltransferase",
               "K11741 sugE; quaternary ammonium \ncompound-resistance protein SugE",
               "K00138 aldB; aldehyde dehydrogenase",
               "K03297 emrE, qac, mmr, smr; small \nmultidrug resistance pump",
               "K08987 putative membrane protein",
               "K09951 cas2; CRISPR-associated \nprotein Cas2")
res$descp <- factor(res$descp, levels = res$descp)
p <- res %>% ggplot(aes(x=LR_stat, y=descp)) +
    geom_bar(stat = "identity", position = "dodge", fill = "skyblue") +
    theme(axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          legend.position = "right") +
    labs(x = "LR_stat", y = "", fill = "", title="") +
    coord_cartesian(xlim = c(150, 260)) +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"),legend.position = "None")

