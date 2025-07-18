library(dplyr)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(lemon)
library(cowplot)
library(reshape2)
library(tidyverse)
library(stringr)

library(ggdensity)
library(ggblanket)
library(ggsci)
library(reshape2)
library(gridExtra)
library(patchwork)
library(cowplot)
library(aplot)
library(ggplotify)

library(RColorBrewer)  

library(ggVennDiagram)

# Gram status, Oxygen preference, Cell shape, Spore production, Motility
### bugbase
traits <- read.csv("Data/traits_bugbase.csv", row.names = 1)
traits_name <- c("Gram_Status", "Oxygen_Preference")
traits[traits == "Anaerobic"] <- "anaerobic"
traits[traits == "Aerobic"] <- "aerobic"
traits[traits == "Facultatively_Anaerobic"] <- "facultatively"
traits[traits == "Gram_Negative"] <- "negative"
traits[traits == "Gram_Positive"] <- "positive"

p <- list()
n <- 1
for (i in traits_name){
    model <- readRDS(file = paste0("Data/bugbase/model_", i,".rds"))
    plot_data <- as.data.frame(model@scoreMN)
    plot_data$Trait <- traits[rownames(plot_data), i]
    plot_data <- plot_data %>% mutate(Trait = if_else(Trait == "Facultatively_Anaerobic", "Facultatively", Trait))
    # results <- readRDS(file = paste0("/home/dongbiao/word_embedding_microbiome/all_data/gut/pheno/plsda_auc/bugbase/model_",
    #                                  i,".rds"))
    # AUC <- round(mean(results$auc_results$AUC), 2)
    R2Y_model <- model@summaryDF[1, 2]
    Q2Y_model <- model@summaryDF[1, 3]
    pR2Y <- model@summaryDF[1, 7]
    pQ2 <- model@summaryDF[1, 8]
    p[[n]] <- plot_data %>% ggplot(aes(x=p1, y=p2, color=Trait)) +
        geom_point(size = 0.5, alpha=0.5)+
        theme_bw(base_size=10) +
        labs(color = NULL, x = "PLS_1", y = "PLS_2", title = i) +
        theme(text = element_text(face = "bold")) +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
        annotate("text", x = -Inf, y = Inf, 
                 label = sprintf("R2=%.2g", R2Y_model),
                 hjust = -0.05, vjust = 1.1, size = 3) +
        annotate("text", x = -Inf, y = Inf, 
                 label = sprintf("Q2=%.2g", Q2Y_model),
                 hjust = -0.05, vjust = 3.1, size = 3) +
        annotate("text", x = Inf, y = Inf, 
                 label = sprintf("pR2=%.2g", pR2Y),
                 hjust = 1, vjust = 1.1, size = 3) +
        annotate("text", x = Inf, y = Inf, 
                 label = sprintf("pQ2=%.2g", pQ2),
                 hjust = 1, vjust = 3.1, size = 3)
    n <- n + 1
}
p <- wrap_plots(p, nrow = 1, ncol = 2)

### traitar
traits <- read.csv("Data/trait_predcit.csv", row.names = 1)

traits[traits == 3] <- 1
traits$Oxygen_Preference <- rep(0, nrow(traits))
traits[traits$Aerobe == 1, "Oxygen_Preference"] <- "anaerobic"
traits[traits$Facultative == 1, "Oxygen_Preference"] <- "aerobic"
traits[traits$Anaerobe == 1, "Oxygen_Preference"] <- "facultatively"
traits$`Gram_Status` <- rep(0, nrow(traits))
traits[traits$`Gram.negative` == 1, "Gram_Status"] <- "negative"
traits[traits$`Gram.positive` == 1, "Gram_Status"] <- "positive"
traits$`cell_shape` <- rep(0, nrow(traits))
traits[traits$Coccus == 1, "cell_shape"] <- "Coccus"
traits[traits$`Bacillus.or.coccobacillus` == 1, "cell_shape"] <- "Bacillus or \ncoccobacillus"
traits[traits == 0] = "no"
traits[traits == 1] = "yes"
traits_name <- colnames(traits)
traits_name <- traits_name[!(traits_name %in% c("Aerobe", "Facultative", "Anaerobe",
                                                "Gram.negative", "Gram.positive"))]

traits_name <- c("Gram_Status", "Oxygen_Preference", "cell_shape", "Spore.formation", "Motile")
traits_name_list <- list()
traits_name_list[["Gram_Status"]] <- "Gram_Status"
traits_name_list[["Oxygen_Preference"]] <- "Oxygen_Preference"
traits_name_list[["cell_shape"]] <- "Cell_Shape"
traits_name_list[["Spore.formation"]] <- "Spore_Formation"
traits_name_list[["Motile"]] <- "Motility"
p <- list()
n <- 1
for (i in traits_name){
    model <- readRDS(file = paste0("Data/traitar/model_", i,".rds"))
    plot_data <- as.data.frame(model@scoreMN)
    plot_data$Trait <- as.factor(traits[rownames(plot_data), i])
    R2Y_model <- model@summaryDF[1, 2]
    Q2Y_model <- model@summaryDF[1, 3]
    pR2Y <- model@summaryDF[1, 7]
    pQ2 <- model@summaryDF[1, 8]
    if (i %in% c("Oxygen.Preference", "Gram.Status")){
        plot_data <- plot_data %>% filter(Trait != "0")
    }
    p[[n]] <- plot_data %>% ggplot(aes(x=p1, y=p2, color=Trait)) +
        geom_point(size = 1, alpha=0.5)+
        theme_bw(base_size=10) +
        labs(color = NULL, x = "PLS_1", y = "PLS_2", title = paste0(traits_name_list[[i]])) +
        theme(text = element_text(face = "bold")) +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
        annotate("text", x = -Inf, y = Inf, 
                 label = sprintf("R2=%.2g", R2Y_model),
                 hjust = -0.05, vjust = 1.1, size = 3) +
        annotate("text", x = -Inf, y = Inf, 
                 label = sprintf("Q2=%.2g", Q2Y_model),
                 hjust = -0.05, vjust = 3.1, size = 3) +
        annotate("text", x = Inf, y = Inf, 
                 label = sprintf("pR2=%.2g", pR2Y),
                 hjust = 1, vjust = 1.1, size = 3) +
        annotate("text", x = Inf, y = Inf, 
                 label = sprintf("pQ2=%.2g", pQ2),
                 hjust = 1, vjust = 3.1, size = 3)
    n <- n + 1
}
p <- wrap_plots(p, nrow = 1, ncol =5)


traits <- read.csv("Data/trait_predcit.csv", row.names = 1, check.names = FALSE)
traits_group <- read.csv("Data/traits.tsv", sep="\t", check.names = FALSE) %>%
    filter(category == "Growth: Sugar")
traits_group <- traits_group[-1, ]
traits_name <- traits_group$accession
traits <- traits[, traits_name]
traits[traits == 0] = "no"
traits[traits == 3] = "yes"
names_1 <- read.csv("Data/trait_predcit.csv", 
                    row.names = 1, check.names = FALSE) %>% colnames()
names_2 <- read.csv("Data/trait_predcit.csv", 
                    row.names = 1) %>% colnames()
names_dict <- list()
for (i in c(1: length(names_1))){
    names_dict[[names_1[i]]] <- names_2[i]
}
R2Y <- c()
Q2 <- c()
pR2Y <- c()
pQ2 <- c()
for (i in traits_name){
    i <- names_dict[[i]]
    model <- readRDS(file = paste0("Data/traitar/model_",
                                   i,".rds"))
    R2Y <- c(R2Y, model@summaryDF[1, 2])
    Q2 <- c(Q2, model@summaryDF[1, 3])
    pR2Y <- c(pR2Y, model@summaryDF[1, 7])
    pQ2 <- c(pQ2, model@summaryDF[1, 8])
}
results_traitar <- data.frame(R2Y = R2Y, Q2 = Q2, pR2Y=pR2Y, pQ2=pQ2, Traits=traits_name) %>% 
    arrange(desc(Q2)) %>% filter(Traits != "Urea hydrolysis") %>% 
    mutate(sig=ifelse(pQ2 <= 0.05, "sig.", "ns"))
results_traitar$group_1 <- rep("Traitar", nrow(results_traitar))
results_traitar$group_2 <- rep("growth: sugar")
results_traitar <- results_traitar %>% arrange(desc(Q2))
p <- list()
n <- 1
traits_name <- results_traitar$Traits
for (i in traits_name){
    title <- i
    model <- readRDS(file = paste0("Data/traitar/model_",
                                   names_dict[[i]],".rds"))
    plot_data <- as.data.frame(model@scoreMN)
    plot_data$Trait <- as.factor(traits[rownames(plot_data), i])
    R2Y_model <- model@summaryDF[1, 2]
    Q2Y_model <- model@summaryDF[1, 3]
    pR2Y <- model@summaryDF[1, 7]
    pQ2 <- model@summaryDF[1, 8]
    p[[n]] <- plot_data %>% ggplot(aes(x=p1, y=p2, color=Trait)) +
        geom_point(size = 1, alpha=0.5)+
        theme_bw(base_size=10) +
        labs(color = NULL, x = "PLS_1", y = "PLS_2", title = title) +
        theme(text = element_text(face = "bold")) +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
        annotate("text", x = -Inf, y = Inf, 
                 label = sprintf("R2=%.2g", R2Y_model),
                 hjust = -0.05, vjust = 1.1, size = 3) +
        annotate("text", x = -Inf, y = Inf, 
                 label = sprintf("Q2=%.2g", Q2Y_model),
                 hjust = -0.05, vjust = 3.1, size = 3) +
        annotate("text", x = Inf, y = Inf, 
                 label = sprintf("pR2=%.2g", pR2Y),
                 hjust = 1, vjust = 1.1, size = 3) +
        annotate("text", x = Inf, y = Inf, 
                 label = sprintf("pQ2=%.2g", pQ2),
                 hjust = 1, vjust = 3.1, size = 3)
    n <- n + 1
}
p <- wrap_plots(p, nrow = 4, ncol =5)


### bacDive
traits <- read.csv("Data/bacDive.csv")
traits <- traits[!duplicated(traits$X16s_ID),]
rownames(traits) <- traits$X16s_ID

co_embedding <- read.csv("../../data/social_niche_embedding_100.txt",
                         row.names = 1, sep=" ", header = FALSE)
accessions_num <- str_split(rownames(co_embedding), "\\.", simplify = TRUE)
accessions_num <- accessions_num[,1]
df <- data.frame(accessions=accessions_num, embed_id=rownames(co_embedding))
inter_id <- intersect(rownames(traits), df$accessions)
df <- df %>% filter(accessions %in% inter_id)
traits <- traits[df$accessions, ]
rownames(traits) <- df$embed_id

traits$`Oxygen.Preference` <- rep(NA, nrow(traits))
traits <- traits %>% mutate(Oxygen.Preference = case_when(aerobe == 1 ~ "aerobic", TRUE ~ Oxygen.Preference))
traits <- traits %>% mutate(Oxygen.Preference = case_when(`facultative.anaerobe` == 1 ~ "facultatively", TRUE ~ Oxygen.Preference))
traits <- traits %>% mutate(Oxygen.Preference = case_when(anaerobe == 1 ~ "anaerobic", TRUE ~ Oxygen.Preference))
traits[traits == 0] = "no"
traits[traits == 1] = "yes"
traits[traits == "mixed"] = NA
traits[traits == "negative;positive"] = NA
traits[traits == "negative;variable"] = NA
traits[traits == ""] = NA
traits[traits == "filament-shaped"] = NA
traits[traits == "oval-shaped"] = NA
traits[traits == "ovoid-shaped"] = NA
traits[traits == "sphere-shaped"] = NA
traits[traits == "spiral-shaped"] = NA

traits_name <- c("gram_stain", "Oxygen.Preference", "cell_shape", "spore_formation", "motility")
traits_name_list <- list()
traits_name_list[["gram_stain"]] <- "Gram_Status"
traits_name_list[["Oxygen.Preference"]] <- "Oxygen_Preference"
traits_name_list[["cell_shape"]] <- "Cell_Shape"
traits_name_list[["spore_formation"]] <- "Spore_Formation"
traits_name_list[["motility"]] <- "Motility"

p <- list()
n <- 1
for (i in traits_name){
    model <- readRDS(file = paste0("Data/bacdive/model_",
                                   i,".rds"))
    plot_data <- as.data.frame(model@scoreMN)
    plot_data$Trait <- as.factor(traits[rownames(plot_data), i])
    R2Y_model <- model@summaryDF[1, 2]
    Q2Y_model <- model@summaryDF[1, 3]
    pR2Y <- model@summaryDF[1, 7]
    pQ2 <- model@summaryDF[1, 8]
    plot_data <- plot_data[!is.na(plot_data$Trait), ]
    p[[n]] <- plot_data %>% ggplot(aes(x=p1, y=p2, color=Trait)) +
        geom_point(size = 1, alpha=0.5)+
        theme_bw(base_size=10) +
        labs(color = NULL, x = "PLS_1", y = "PLS_2", title = paste0(traits_name_list[[i]])) +
        theme(text = element_text(face = "bold")) +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
        annotate("text", x = -Inf, y = Inf, 
                 label = sprintf("R2=%.2g", R2Y_model),
                 hjust = -0.05, vjust = 1.1, size = 3) +
        annotate("text", x = -Inf, y = Inf, 
                 label = sprintf("Q2=%.2g", Q2Y_model),
                 hjust = -0.05, vjust = 3.1, size = 3) +
        annotate("text", x = Inf, y = Inf, 
                 label = sprintf("pR2=%.2g", pR2Y),
                 hjust = 1, vjust = 1.1, size = 3) +
        annotate("text", x = Inf, y = Inf, 
                 label = sprintf("pQ2=%.2g", pQ2),
                 hjust = 1, vjust = 3.1, size = 3)
    n <- n + 1
}
p <- wrap_plots(p, nrow = 1, ncol = 5)


### 
traits <- read.csv("Data/bacDive.csv")
traits <- traits[!duplicated(traits$`X16s_ID`), ]
rownames(traits) <- traits$`X16s_ID`

co_embedding <- read.csv("../../data/social_niche_embedding_100.txt",
                         row.names = 1, sep=" ", header = FALSE)
accessions_num <- str_split(rownames(co_embedding), "\\.", simplify = TRUE)
accessions_num <- accessions_num[,1]
df <- data.frame(accessions=accessions_num, embed_id=rownames(co_embedding))
inter_id <- intersect(rownames(traits), df$accessions)
df <- df %>% filter(accessions %in% inter_id)
traits <- traits[df$accessions, ]
rownames(traits) <- df$embed_id
traits[traits == "NA"] <- NA
traits[traits == "-"] <- "no"
traits[traits == "+"] <- "yes"
traits[traits == ""] = NA
traits[traits == "+;NA"] = NA
remove_id <- c("X16s_ID", "aerobe", "facultative.anaerobe", "anaerobe")
traits <- traits[, !(colnames(traits) %in% remove_id)]
traits <- traits[, colSums(is.na(traits)) < nrow(traits)]

fid <- rownames(traits)
co_embedding <- co_embedding[fid, ]

co_embedding_shuffled <- co_embedding %>% 
    mutate(across(everything(), ~ sample(.x)))

traits_name <- colnames(traits)
traits_name <- traits_name[10:length(traits_name)]
agg_bac <- read.csv("Data/agg_bac.csv") %>%
    filter(level_2 %in% c("assimilation", "builds_acid_from"))
agg_bac <- agg_bac[!duplicated(agg_bac$terms), ] %>% filter(terms %in% traits_name)

R2Y <- c()
Q2 <- c()
pR2Y <- c()
pQ2 <- c()
for (i in agg_bac$terms){
    tryCatch({
        model <- readRDS(file = paste0("Data/bacdive/model_",
                                       i,".rds"))
        if (sum(table(traits[rownames(model@scoreMN), i]) > 5) ==2){
            R2Y <- c(R2Y, model@summaryDF[1, 2])
            Q2 <- c(Q2, model@summaryDF[1, 3])
            pR2Y <- c(pR2Y, model@summaryDF[1, 7])
            pQ2 <- c(pQ2, model@summaryDF[1, 8])
        } else {
            R2Y <<- c(R2Y, 0)
            pR2Y <<- c(pR2Y, 0)
            Q2 <<- c(Q2, 0)
            pQ2 <<- c(pQ2, 0)
        }
    }, error = function(e) {
        R2Y <<- c(R2Y, 0)
        pR2Y <<- c(pR2Y, 0)
        Q2 <<- c(Q2, 0)
        pQ2 <<- c(pQ2, 0)
    })
}
results_bac_dive <- data.frame(R2Y = R2Y, Q2 = Q2, pR2Y=pR2Y, 
                               pQ2=pQ2, group_2=agg_bac$level_2,
                               Traits=str_to_title(agg_bac$level_3)) %>% 
    arrange(desc(Q2)) %>% filter(pR2Y != 0) %>%
    mutate(sig=ifelse(pQ2 <= 0.05, "sig.", "ns"))
results_bac_dive$group_1 <- rep("BacDive", nrow(results_bac_dive))
results <- rbind(results_traitar[, colnames(results_bac_dive)], results_bac_dive)
results$group_1 <- factor(results$group_1, levels = c("Traitar", "BacDive"))
p <- results %>% ggplot(aes(x=R2Y, y=Q2)) +
    geom_point(aes(color=sig, shape=group_2)) +
    theme_bw() +
    labs(x="R2", shape="", color="") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_text_repel(aes(label = Traits),
                    size = 2, # 设置标签文字大小
                    max.overlaps = Inf, # 确保所有标签都显示
                    box.padding = 0.5, # 标签周围的填充空间，增加距离
                    min.segment.length = 0 # 即使标签离点很近也画线
    ) +
    scale_color_manual(values = c("sig." = "#E41A1C", "ns" = "gray")) +
    theme(text = element_text(face = "bold")) +
    facet_wrap(.~group_1, scales = "free")


R2Y <- c()
Q2 <- c()
for (i in agg_bac$terms){
    tryCatch({
        model <- readRDS(file = paste0("Data/bacdive/model_",
                                       i,".rds"))
        if (sum(table(traits[rownames(model@scoreMN), i]) > 5) ==2){
        R2Y <- c(R2Y, model@summaryDF[1, 2])
        Q2 <- c(Q2, model@summaryDF[1, 3])
        } else {
            R2Y <<- c(R2Y, 0)
            Q2 <<- c(Q2, 0)
        }
    }, error = function(e) {
        R2Y <<- c(R2Y, 0)
        Q2 <<- c(Q2, 0)
    })
}
agg_bac$R2Y <- R2Y
agg_bac$Q2 <- Q2

agg_bac <- agg_bac %>% filter(R2Y != 0) %>% arrange(level_2, desc(Q2))
agg_bac$name <- str_to_title(agg_bac$level_3)
terms_id_dict <- list()
for (i in 1:nrow(agg_bac)){
    terms_id_dict[[agg_bac[i, 1]]] <- agg_bac[i, 9]
}

### assimilation
agg_bac_assimilation <- agg_bac %>% filter(level_2 == "assimilation") %>% 
    arrange(level_2, desc(Q2))

p <- list()
n <- 1
for (i in agg_bac_assimilation$terms){
    model <- readRDS(file = paste0("Data/bacdive/model_", i,".rds"))
    plot_data <- as.data.frame(model@scoreMN)
    plot_data$Trait <- as.factor(traits[rownames(plot_data), i])
    R2Y_model <- model@summaryDF[1, 2]
    Q2Y_model <- model@summaryDF[1, 3]
    pR2Y <- model@summaryDF[1, 7]
    pQ2 <- model@summaryDF[1, 8]
    p[[n]] <- plot_data %>% ggplot(aes(x=p1, y=p2, color=Trait)) +
        geom_point(size = 1, alpha=0.5)+
        theme_bw(base_size=10) +
        labs(color = NULL, x = "PLS_1", y = "PLS_2", title = paste0(terms_id_dict[[i]])) +
        theme(text = element_text(face = "bold")) +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
        annotate("text", x = -Inf, y = Inf, 
                 label = sprintf("R2=%.2g", R2Y_model),
                 hjust = -0.05, vjust = 1.1, size = 3) +
        annotate("text", x = -Inf, y = Inf, 
                 label = sprintf("Q2=%.2g", Q2Y_model),
                 hjust = -0.05, vjust = 3.1, size = 3) +
        annotate("text", x = Inf, y = Inf, 
                 label = sprintf("pR2=%.2g", pR2Y),
                 hjust = 1, vjust = 1.1, size = 3) +
        annotate("text", x = Inf, y = Inf, 
                 label = sprintf("pQ2=%.2g", pQ2),
                 hjust = 1, vjust = 3.1, size = 3)
    n <- n + 1
}
p <- wrap_plots(p, nrow = 1, ncol = 4)


### builds_acid_from
agg_bac_builds_acid_from <- agg_bac %>% filter(level_2 == "builds_acid_from") %>% 
    arrange(level_2, desc(Q2))

p <- list()
n <- 1
for (i in agg_bac_builds_acid_from$terms){
    model <- readRDS(file = paste0("Data/bacdive/model_", i,".rds"))
    plot_data <- as.data.frame(model@scoreMN)
    plot_data$Trait <- as.factor(traits[rownames(plot_data), i])
    R2Y_model <- model@summaryDF[1, 2]
    Q2Y_model <- model@summaryDF[1, 3]
    pR2Y <- model@summaryDF[1, 7]
    pQ2 <- model@summaryDF[1, 8]
    p[[n]] <- plot_data %>% ggplot(aes(x=p1, y=p2, color=Trait)) +
        geom_point(size = 1, alpha=0.5)+
        theme_bw(base_size=10) +
        labs(color = NULL, x = "PLS_1", y = "PLS_2", title = paste0(terms_id_dict[[i]])) +
        theme(text = element_text(face = "bold")) +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
        annotate("text", x = -Inf, y = Inf, 
                 label = sprintf("R2=%.2g", R2Y_model),
                 hjust = -0.05, vjust = 1.1, size = 3) +
        annotate("text", x = -Inf, y = Inf, 
                 label = sprintf("Q2=%.2g", Q2Y_model),
                 hjust = -0.05, vjust = 3.1, size = 3) +
        annotate("text", x = Inf, y = Inf, 
                 label = sprintf("pR2=%.2g", pR2Y),
                 hjust = 1, vjust = 1.1, size = 3) +
        annotate("text", x = Inf, y = Inf, 
                 label = sprintf("pQ2=%.2g", pQ2),
                 hjust = 1, vjust = 3.1, size = 3)
    n <- n + 1
}
p <- wrap_plots(p, nrow = 4, ncol = 5)

