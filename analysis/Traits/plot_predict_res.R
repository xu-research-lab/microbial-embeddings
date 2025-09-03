library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(ggpubr)

brewer.pal(9, c("Set1"))
c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628","#F781BF", "#999999")

df_1 <- read.csv("Data/auc_res.csv") %>%
    filter(datasize == "210,000") %>% filter(group == "times_1")
df_1 <- df_1[, c("auc", "traits_type", "tax")]
colnames(df_1) <- c("AUC", "traits", "tax")
df_1$group <- "traits_1"
df_2 <- read.csv("Data/predict_metabolics_res.csv") %>%
    arrange(desc(AUC))
df_2$group <- "traits_2"

df <- rbind(df_1, df_2)
df$traits <- factor(df$traits, levels = c("Oxygen_Preference", "Gram_Status", "Motility", 
                                          "Spore_Formation", "Lactose", "Melibiose", "Glycerol", 
                                          "Maltose", "Salicin", "Sucrose", "Trehalose"))
df_summary <- df %>%
    group_by(traits, group) %>% 
    summarise(
        AUC_mean = mean(AUC),
        AUC_sd = sd(AUC),
        .groups = "drop" 
    )
p_1 <- ggplot() +
    geom_col(
        data = df_summary,
        aes(x = traits, y = AUC_mean, color = group), 
        alpha = 0, width = 0.7,
        position = position_dodge(width = 0)
    ) +
    geom_errorbar(
        data = df_summary,
        aes(x = traits, ymin = AUC_mean - AUC_sd, ymax = AUC_mean + AUC_sd, color = group),
        width = 0.2,
        position = position_dodge(width = 0.9)
    ) +
    geom_point(
        data = df, aes(x = traits, y = AUC, color = tax),
        size = 2.5,
        alpha = 0.8) +
    coord_cartesian(ylim = c(0.45, 0.95)) +
    scale_color_manual(values = c("traits_1"="#E41A1C", "traits_2"="#377EB8", 
                                  "Bacillota"="#4DAF4A", "Actinomycetota"="#984EA3", "Bacteroidota"="#A65628")) +
    labs(x = "", y = "AUC", title = "") +
    theme_bw(base_size = 10) +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1))

### pretrain data size
auc_res <- read.csv("Data/auc_res.csv")
auc_res <- auc_res %>% group_by(traits_type, tax, datasize) %>%
    summarise(AUC = mean(auc))
auc_summary <- auc_res %>% group_by(traits_type, datasize) %>%
    summarise(AUC_mean = mean(AUC))
auc_res <- merge(auc_res, auc_summary, by=c("traits_type", "datasize"))
auc_res$datasize <- factor(auc_res$datasize, levels = c("10,000", "20,000", "40,000", "80,000", "160,000", "210,000"))
auc_res <- auc_res %>% filter(datasize %in% c("80,000", "160,000", "210,000"))
auc_res$traits_type <- factor(auc_res$traits_type, levels = c("Oxygen_Preference", "Gram_Status", "Motility", "Spore_Formation"))

p_2 <- auc_res %>% ggplot(aes(x=datasize, y = AUC)) +
    geom_point(aes(color=tax), size=2) +
    geom_line(aes(group=1, y=AUC_mean), size=0.8) +
    labs(x = "", y = "AUC", color="")+
    facet_wrap(~traits_type, nrow = 1) +
    theme_bw(base_size = 10) +
    coord_cartesian(ylim = c(0.45, 0.95)) +
    scale_color_manual(values = c("Bacillota"="#4DAF4A", "Actinomycetota"="#984EA3", "Bacteroidota"="#A65628")) +
    theme(legend.position = "none") +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = -45, hjust = 0))

