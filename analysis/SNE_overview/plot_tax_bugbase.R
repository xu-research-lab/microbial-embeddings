library(tidyverse)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)

tax <- read.delim("../Pretraining_data_profile/Data/taxmap_slv_ssu_ref_nr_138.2.txt", 
                  sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
acc <- paste(tax[[1]], tax[[2]], tax[[3]], sep = ".")
tax_split <- strsplit(tax$path, ";")
max_len <- max(sapply(tax_split, length))
tax_split <- lapply(tax_split, function(x) {
    length(x) <- max_len 
    return(x)
})

tax <- as.data.frame(do.call(rbind, tax_split))
rownames(tax) <- acc

tax <- tax[, 1:7]
colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

t_sne_1 <- read.csv("Data/t_sne_co.csv", check.names = FALSE)
t_sne_1$group <- rep("SNE", nrow(t_sne_1))
t_sne_2 <- read.csv("Data/t_sne_phylo.csv", check.names = FALSE)
t_sne_2$group <- rep("PhyloE", nrow(t_sne_2))
t_sne <- rbind(t_sne_1, t_sne_2)

tax <- tax[unique(t_sne$`0`), ]
taxa_levels <- c("Phylum", "Class", "Order", "Family", "Genus")
for (level in taxa_levels){
    tax_count <- sort(table(tax[, level]), decreasing = TRUE)
    top_14 <- names(head(tax_count, 14))
    top_14 <- top_14[top_14 != "Incertae Sedis"]
    tax_labels <- case_when(
        tax[, level] %in% top_14 ~ tax[, level],
        TRUE ~ "Others")
    
    t_sne$tax <- factor(c(tax_labels, tax_labels), levels = c(top_14, "Others"))
    
    col <- c(brewer.pal(12,'Paired'), brewer.pal(4,'Pastel2'))
    col <- setNames(col, levels(t_sne$tax))
    col["Others"] <- "grey"
    t_sne$group <- factor(t_sne$group, levels = c("SNE", "PhyloE"))
    p1 <- t_sne %>% ggplot(aes(x = `t-SNE1`, y = `t-SNE2`, color = tax))+
        facet_wrap(~group, scales = "free")+
        geom_point(alpha = 0.5, size = 1) +
        theme_bw(base_size=12) +
        scale_color_manual(values = col) +
        labs(x = "t-SNE 1", y = "t-SNE 2", color = level) +
        theme(text = element_text(face = "bold"),
             legend.title = element_text(size = 12),
             legend.text = element_text(size = 12),
             plot.title = element_text(hjust = 0.5)) +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
    
    ggsave(paste0("Figures/", level, ".png"), p1,
           width = 25, height = 10, units = "cm")
}

### Bugbase
df <- read.table("../Traits_function/Data/traits_precalculated.txt", 
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)
t_sne_1 <- read.csv("Data/t_sne_co_bugbase.csv")
t_sne_2 <- read.csv("Data/t_sne_phy_bugbase.csv")
t_sne_1$group <- rep("SNE", nrow(t_sne_1))
t_sne_2$group <- rep("PhyloE", nrow(t_sne_2))
t_sne <- rbind(t_sne_1, t_sne_2)
id <- t_sne_1[, 1]

#Aerobic_Anaerobic_Facultatively
plot <- t_sne
plot$Aerobic <- rep(df[id, ]$Aerobic, 2)
plot$Anaerobic <- rep(df[id, ]$Anaerobic, 2)
plot$`Facultatively \nAnaerobic` <- rep(df[id, ]$Facultatively_Anaerobic, 2)
plot <- plot %>% melt(id.vars = c("X0", "t.SNE1", "t.SNE2", "group" )) %>%
    filter(value > 0)

plot$group <- factor(plot$group, levels = c("SNE", "PhyloE"))
p <- ggplot(plot, aes(x = t.SNE1, y = t.SNE2, color = variable)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_brewer(palette = "Set2") +
    labs(x = "t-SNE1", y = "t-SNE2", color = "") +
    theme_bw(base_size = 12) +
    facet_wrap(~group, scales = "free")+
    theme(
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        text = element_text(face = "bold")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggsave("Figures/Aerobic_Anaerobic_Facultatively.png",
       width = 25, height = 10, units = "cm")

#Gram_negative、Gram_positive
plot <- t_sne
plot$Gram_Negative <- rep(df[id, ]$Gram_Negative, 2)
plot$Gram_Positive <- rep(df[id, ]$Gram_Positive, 2)
plot <- plot %>% melt(id.vars = c("X0", "t.SNE1", "t.SNE2", "group" )) %>%
    filter(value > 0)
plot$group <- factor(plot$group, levels = c("SNE", "PhyloE"))
p <- ggplot(plot, aes(x = t.SNE1, y = t.SNE2, color = variable)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_brewer(palette = "Set2") +
    labs(x = "t-SNE1", y = "t-SNE2", color = "") +
    theme_bw(base_size = 12) +
    facet_wrap(~group, scales = "free")+
    theme(
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        text = element_text(face = "bold")
    ) + 
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggsave("Figures/Gram.png", width = 25, height = 10, units = "cm")

