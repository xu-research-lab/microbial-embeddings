library(tidyverse)
library(ggplot2)

metadata <- read.csv("/home/dongbiao/all_study/data/metadata.csv")
metadata <- metadata[!duplicated(metadata$study), ]

study_sample_num <- read.csv("/home/dongbiao/all_study/study_sample_num.csv")
study_disease <- list()
for (i in c(1: nrow(metadata))){
    study_disease[[metadata$study[i]]] <- metadata$disease.name[i]
}

study_num <- list()
for (i in c(1: nrow(study_sample_num))){
    study_num[[study_sample_num$study[i]]] <- study_sample_num$sample[i]
}
table_1 <- read.csv("/home/dongbiao/all_study/attention_model_5fold_res.csv")
table_1$group <- rep("Attention", nrow(table_1))
table_2 <- read.csv("/home/dongbiao/all_study/rf_5fold_res.csv")
table_2$group <- rep("RF", nrow(table_2))
table <- rbind(table_1, table_2)

table_mean_sd <- table %>% group_by(study, group) %>% summarise(AUC_mean = mean(auc), AUC_sd = sd(auc))
table <- merge(table, table_mean_sd, by = c("study", "group"))
table$disease <- unlist(lapply(table$study, function(x) as.vector(unlist(study_disease[x]))))
table$num <- unlist(lapply(table$study, function(x) as.vector(unlist(study_num[x]))))
table$group1 <- paste(table$disease, "\n", table$study, "_", table$num)


# disease_factor <- table %>% filter(group == "Attention") %>% arrange(desc(AUC_mean))
# disease_factor <- unique(disease_factor$disease)
# table$disease <- factor(table$disease, levels = disease_factor)

table %>% ggplot(aes(x = group, y = auc)) +
    geom_point(size=3) +
    facet_wrap(.~group1, nrow = 5, scales = "free") +
    geom_line(aes(group=fold))+
    theme_bw(base_size = 14)+
    theme(text = element_text(face = "bold"), 
          legend.text = element_text(size = 14, color = "black", face = "bold"),
          legend.box.margin = margin(l = 0.1, unit = "in"),
          legend.title.align = 0.5 
    ) +
    labs(x="", y="AUC")
    
