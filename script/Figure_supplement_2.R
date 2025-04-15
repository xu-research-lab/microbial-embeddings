library(tidyverse)
library(biomformat)
library(dplyr)
library(reshape2)
library(cowplot)
library(ggplot2)
library(aplot)
library(patchwork)
library(gridExtra)
library(grid)



AGP_table <- read_biom("/home/dongbiao/word_embedding_microbiome/Classification_prediction/data/AGP_test.biom")
table <- do.call(data.frame, AGP_table$data)
colnames(table) <- rownames(AGP_table)
rm(AGP_table)

sample_1 <- data.frame(table[1,]) %>% t() %>% as.data.frame()
colnames(sample_1) <- "sample1"
sample_1 <- sample_1[sample_1$sample1 > 0, ]

sample_1_totalsum <- sample_1 / sum(sample_1)

# sorted_indices <- order(sample_1)
sample_1_percentile <- c(1:length(sample_1)) / length(sample_1)

### abundance percentile: min(prti, prtj) * (1 - abs(prti - prtj))
otu_1 <- c()
otu_2 <- c()
co_occurrence <- c()
n <- length(sample_1_percentile)
for (i in 1 : (n - 1)){
    for (j in (i+1) : n){
        otu_1 <- c(otu_1, max(sample_1_percentile[i], sample_1_percentile[j]))
        otu_2 <- c(otu_2, min(sample_1_percentile[i], sample_1_percentile[j]))
        co = min(sample_1_percentile[i], sample_1_percentile[j]) * (1 - abs(sample_1_percentile[i] - sample_1_percentile[j]))
        co_occurrence <- c(co_occurrence, co)
    }
}
coocc_range <- c(0, 1)
value_range <- c(0, 1)
scale_factor <- (coocc_range[2] - coocc_range[1]) / (value_range[2] - value_range[1])


res_percentile <- data.frame(Max_PCT = otu_1, Min_PCT = otu_2, co_occurrence = co_occurrence)
res_percentile <- res_percentile %>% arrange(co_occurrence, otu_2)
res_percentile$id <- c(1 : nrow(res_percentile))
res_percentile_all <- res_percentile %>% melt(id.vars = "id")
res_percentile_all$group <- rep("Co-occurrence in a sample", n=nrow(res_percentile_all))

res_percentile_part <- res_percentile %>% filter(Min_PCT < 0.26) %>% filter(Min_PCT > 0.24)
res_percentile_part <- res_percentile_part %>% melt(id.vars = "id")
res_percentile_part$group <- rep("The part Co-occurrence\nin a sample", n=nrow(res_percentile_part))

res_percentile <- rbind(res_percentile_all, res_percentile_part)
p1 <- res_percentile %>% 
    ggplot(aes(x=id)) +
    geom_point(aes(y=value, color=variable), size=0.1) +
    theme_bw(base_size = 14) +
    facet_wrap(.~group, scales = "free_x")+
    labs(x="OTUs co-occurrence pairs") +
    scale_color_manual(values = c("Max_PCT" = "#E84445", "Min_PCT" = "#B2DBB9", "co_occurrence"="#95BCE5")) +
    scale_y_continuous(
        name = "Abundance",
        sec.axis = sec_axis(~. * scale_factor + coocc_range[1] - value_range[1] * scale_factor, name = "Co-occurrence Value")
    ) +
    theme(
        legend.title = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold")
    )

### braycurtis percentile: 2 * min(prti, prtj) / (prti + prtj)
otu_1 <- c()
otu_2 <- c()
co_occurrence <- c()
n <- length(sample_1_percentile)
for (i in 1 : (n - 1)){
    for (j in (i+1) : n){
        otu_1 <- c(otu_1, max(sample_1_percentile[i], sample_1_percentile[j]))
        otu_2 <- c(otu_2, min(sample_1_percentile[i], sample_1_percentile[j]))
        co = 2 * min(sample_1_percentile[i], sample_1_percentile[j]) / (sample_1_percentile[i] + sample_1_percentile[j])
        co_occurrence <- c(co_occurrence, co)
    }
}
coocc_range <- c(0, 1)
value_range <- c(0, 1)
scale_factor <- (coocc_range[2] - coocc_range[1]) / (value_range[2] - value_range[1])


res_percentile <- data.frame(Max_PCT = otu_1, Min_PCT = otu_2, co_occurrence = co_occurrence)
res_percentile <- res_percentile %>% arrange(co_occurrence, otu_2)
res_percentile$id <- c(1 : nrow(res_percentile))
res_percentile_all <- res_percentile %>% melt(id.vars = "id")
res_percentile_all$group <- rep("Co-occurrence in a sample", n=nrow(res_percentile_all))

res_percentile_part <- res_percentile %>% filter(Min_PCT < 0.26) %>% filter(Min_PCT > 0.24)
res_percentile_part <- res_percentile_part %>% melt(id.vars = "id")
res_percentile_part$group <- rep("The part Co-occurrence\nin a sample", n=nrow(res_percentile_part))

res_percentile <- rbind(res_percentile_all, res_percentile_part)
p2 <- res_percentile %>% 
    ggplot(aes(x=id)) +
    geom_point(aes(y=value, color=variable), size=0.1) +
    theme_bw(base_size = 14) +
    facet_wrap(.~group, scales = "free_x")+
    labs(x="OTUs co-occurrence pairs") +
    scale_color_manual(values = c("Max_PCT" = "#E84445", "Min_PCT" = "#B2DBB9", "co_occurrence"="#95BCE5")) +
    scale_y_continuous(
        name = "Abundance",
        sec.axis = sec_axis(~. * scale_factor + coocc_range[1] - value_range[1] * scale_factor, name = "Co-occurrence Value")
    ) +
    theme(
        legend.title = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold")
    )

### Russeellrao weight: 1 / |prti - prtj|
otu_1 <- c()
otu_2 <- c()
co_occurrence <- c()
n <- length(sample_1_percentile)
sorted_indices <- c(1:n)
for (i in 1 : (n - 1)){
    for (j in (i+1) : n){
        otu_1 <- c(otu_1, max(sample_1_percentile[i], sample_1_percentile[j]))
        otu_2 <- c(otu_2, min(sample_1_percentile[i], sample_1_percentile[j]))
        co = 1 / abs(sorted_indices[i] - sorted_indices[j])
        co_occurrence <- c(co_occurrence, co)
    }
}
coocc_range <- c(0, 1)
value_range <- c(0, 1)
scale_factor <- (coocc_range[2] - coocc_range[1]) / (value_range[2] - value_range[1])


res_percentile <- data.frame(Max_PCT = otu_1, Min_PCT = otu_2, co_occurrence = co_occurrence)
res_percentile <- res_percentile %>% arrange(co_occurrence, otu_2)
res_percentile$id <- c(1 : nrow(res_percentile))
res_percentile_all <- res_percentile %>% melt(id.vars = "id")
res_percentile_all$group <- rep("Co-occurrence in a sample", n=nrow(res_percentile_all))

res_percentile_part <- res_percentile %>% filter(Min_PCT < 0.26) %>% filter(Min_PCT > 0.24)
res_percentile_part <- res_percentile_part %>% melt(id.vars = "id")
res_percentile_part$group <- rep("The part Co-occurrence\nin a sample", n=nrow(res_percentile_part))

res_percentile <- rbind(res_percentile_all, res_percentile_part)
p3 <- res_percentile %>% 
    ggplot(aes(x=id)) +
    geom_point(aes(y=value, color=variable), size=0.1) +
    theme_bw(base_size = 14) +
    facet_wrap(.~group, scales = "free_x")+
    labs(x="OTUs co-occurrence pairs") +
    scale_color_manual(values = c("Max_PCT" = "#E84445", "Min_PCT" = "#B2DBB9", "co_occurrence"="#95BCE5")) +
    scale_y_continuous(
        name = "Abundance",
        sec.axis = sec_axis(~. * scale_factor + coocc_range[1] - value_range[1] * scale_factor, name = "Co-occurrence Value")
    ) +
    theme(
        legend.title = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold")
    )

p <- plot_grid(p1, p2, p3,
               align=c("h"),
               labels = c('a', 'b', 'c'),
               scale = c(1, 1, 1),
               nrow = 3, ncol=1, plot=FALSE)

extract_legend <- function(plot) {
    g <- ggplotGrob(plot)
    legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
    return(legend)
}


legend_data <- data.frame(
    Group = c("Hight abund.", "Low abund.", "Co-occur")
)

p_color_legend <- ggplot(legend_data, aes(x = Group, y = 1, fill = Group)) +
    geom_tile() +
    scale_fill_manual(values = c("Hight abund."="#E84445", "Low abund."="#B2DBB9", "Co-occur"="#95BCE5")) +
    labs(fill = "") +
    theme_void() +
    theme(
        legend.position = "top",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)
    )
p_color_legend <- extract_legend(p_color_legend)

p_color_legend <- grid.arrange(p_color_legend, ncol = 1)

p <- p / p_color_legend
p <- p + plot_layout(
    heights = c(2, 0.1)
)

ggsave("/home/dongbiao/word_embedding_microbiome/result/supplement_figure_2.png", p,
       width = 22, height = 30, units = "cm")
