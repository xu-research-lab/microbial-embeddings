library(tidyverse)
library(biomformat)
library(dplyr)
library(reshape2)
library(cowplot)
library(ggplot2)

AGP_table <- read_biom("/home/dongbiao/word_embedding_microbiome/Classification_prediction/data/AGP_test.biom")
table <- do.call(data.frame, AGP_table$data)
colnames(table) <- rownames(AGP_table)
rm(AGP_table)

sample_1 <- data.frame(table[1,]) %>% t() %>% as.data.frame()
colnames(sample_1) <- "sample1"
sample_1 <- sample_1[sample_1$sample1 > 0, ]

sample_1_totalsum <- sample_1 / sum(sample_1)

sorted_indices <- order(sample_1)
sample_1_percentile <- sorted_indices / max(sorted_indices)

### abundance percentile: min(prti, prtj) * (1 - abs(prti - prtj))
otu_1 <- c()
otu_2 <- c()
co_occurrence_1 <- c()
n <- length(sample_1_percentile)
for (i in 1 : (n - 1)){
    for (j in (i+1) : n){
        otu_1 <- c(otu_1, max(sample_1_percentile[i], sample_1_percentile[j]))
        otu_2 <- c(otu_2, min(sample_1_percentile[i], sample_1_percentile[j]))
        co = min(sample_1_percentile[i], sample_1_percentile[j]) * (1 - abs(sample_1_percentile[i] - sample_1_percentile[j]))
        co_occurrence_1 <- c(co_occurrence_1, co)
    }
}
res_percentile <- data.frame(Max_PCT = otu_1, Min_PCT = otu_2, co_occurrence = co_occurrence_1)
res_percentile <- res_percentile %>% arrange(co_occurrence_1, otu_2)
res_percentile$id <- c(1 : nrow(res_percentile))
res_percentile <- res_percentile %>% melt(id.vars = "id")


otu_1 <- c()
otu_2 <- c()
co_occurrence_2 <- c()
n <- length(sample_1_percentile)
for (i in 1 : (n - 1)){
    for (j in (i+1) : n){
        otu_1 <- c(otu_1, max(sample_1_percentile[i], sample_1_percentile[j]))
        otu_2 <- c(otu_2, min(sample_1_percentile[i], sample_1_percentile[j]))
        co = 2 * min(sample_1_percentile[i], sample_1_percentile[j]) / (sample_1_percentile[i] + sample_1_percentile[j])
        co_occurrence_2 <- c(co_occurrence_2, co)
    }
}

# bc_percentile <- data.frame(Max_PCT = otu_1, Min_PCT = otu_2, co_occurrence = co_occurrence_2)
# bc_percentile <- bc_percentile %>% arrange(co_occurrence_2, otu_2)
# bc_percentile$id <- c(1 : nrow(bc_percentile))
# bc_percentile <- bc_percentile %>% melt(id.vars = "id")

p1 <- res_percentile %>% ggplot(aes(x=id, y=value)) +
    geom_point(aes(color=variable), size=0.1) +
    theme_bw(base_size = 14) +
    labs(x="Species co-occurrence pairs", y="Value") +
    scale_color_manual(values = c("Max_PCT" = "#E84445", "Min_PCT" = "#B2DBB9", "co_occurrence"="#95BCE5")) +
    theme(legend.title = element_blank(),
          legend.position = "none",
          axis.text.x = element_blank(),
          text = element_text(face = "bold"), 
          axis.title = element_text(face = "bold"),
          axis.text = element_text(face = "bold"))

# p2 <- bc_percentile %>% ggplot(aes(x=id, y=value)) +
#     geom_point(aes(color=variable), size=0.1) +
#     theme_bw(base_size = 14) +
#     labs(x="Species co-occurrence pairs", y="Value") +
#     scale_color_manual(values = c("Max_PCT" = "#E84445", "Min_PCT" = "#B2DBB9", "co_occurrence"="#95BCE5")) +
#     theme(legend.title = element_blank(),
#           legend.position = "none",
#           axis.text.x = element_blank(),
#           text = element_text(face = "bold"), 
#           axis.title = element_text(face = "bold"),
#           axis.text = element_text(face = "bold"))

res_percentile <- data.frame(Max_PCT = otu_1, Min_PCT = otu_2, co_occurrence = co_occurrence)
res_percentile <- res_percentile %>% arrange(co_occurrence, otu_2)
res_percentile$id <- c(1 : nrow(res_percentile))
res_percentile_part <- res_percentile %>% filter(Min_PCT < 0.26) %>% filter(Min_PCT > 0.24)
res_percentile_part <- res_percentile_part %>% melt(id.vars = "id")

p2 <-  res_percentile_part %>% ggplot(aes(x=id, y=value)) +
    geom_point(aes(color=variable), size=0.1) +
    theme_bw(base_size = 14) +
    labs(x="Species co-occurrence pairs", y="Value") +
    scale_color_manual(values = c("Max_PCT" = "#E84445", "Min_PCT" = "#B2DBB9", "co_occurrence"="#95BCE5")) +
    theme(legend.title = element_blank(),
          legend.position = "right",
          axis.text.x = element_blank(),
          text = element_text(face = "bold"), 
          axis.title = element_text(face = "bold"),
          axis.text = element_text(face = "bold"))

otu_1 <- c()
otu_2 <- c()
co_occurrence <- c()
n <- length(sample_1_totalsum)
for (i in 1 : (n - 1)){
    for (j in (i+1) : n){
        otu_1 <- c(otu_1, max(sample_1_totalsum[i], sample_1_totalsum[j]))
        otu_2 <- c(otu_2, min(sample_1_totalsum[i], sample_1_totalsum[j]))
        co = min(sample_1_totalsum[i], sample_1_totalsum[j]) * (1 - abs(sample_1_totalsum[i] - sample_1_totalsum[j]))
        co_occurrence <- c(co_occurrence, co)
    }
}
res_ra <- data.frame(Max_RA = otu_1, Min_RA = otu_2, co_occurrence = co_occurrence)
res_ra <- res_ra %>% arrange(co_occurrence, Min_RA)
res_ra$id <- c(1 : nrow(res_ra))
res_ra <- res_ra %>% melt(id.vars = "id")
p3 <- res_ra %>% ggplot(aes(x=id, y=log10(value))) +
    geom_point(aes(color=variable), size=0.1) +
    theme_bw(base_size = 14) +
    labs(x="Species co-occurrence pairs", y="log10(Value)") +
    scale_color_manual(values = c("Max_RA" = "#E84445", "Min_RA" = "#B2DBB9", "co_occurrence"="#95BCE5")) +
    theme(legend.title = element_blank(),
          legend.position = "none",
          axis.text.x = element_blank(),
          text = element_text(face = "bold"), 
          axis.title = element_text(face = "bold"),
          axis.text = element_text(face = "bold"))

res_ra <- data.frame(Max_RA = otu_1, Min_RA = otu_2, co_occurrence = co_occurrence)
res_ra <- res_ra %>% arrange(co_occurrence, Min_RA)
res_ra$id <- c(1 : nrow(res_ra))
res_ra_part <- res_ra %>% filter(Min_RA < 0.0002) %>% filter(Min_RA > 0.0001)
res_ra_part <- res_ra_part %>% melt(id.vars = "id")

p4 <-  res_ra_part %>% ggplot(aes(x=id, y=log10(value))) +
    geom_point(aes(color=variable), size=0.1) +
    theme_bw(base_size = 14) +
    labs(x="Species co-occurrence pairs", y="log10(Value)") +
    scale_color_manual(values = c("Max_RA" = "#E84445", "Min_RA" = "#B2DBB9", "co_occurrence"="#95BCE5")) +
    theme(legend.title = element_blank(),
          legend.position = "right",
          axis.text.x = element_blank(),
          text = element_text(face = "bold"), 
          axis.title = element_text(face = "bold"),
          axis.text = element_text(face = "bold"))

p <- plot_grid(p1, p2, p3, p4,
               labels = c('a', 'b', 'c', 'd'),
               align="hv",
               scale = c(1, 1, 1, 1),
               nrow = 2, ncol=2, plot=FALSE, rel_widths = c(1, 1.5))

ggsave("/home/dongbiao/word_embedding_microbiome/result/methods_co_occurrence.png", p,
       width = 22, height = 15, units = "cm")
