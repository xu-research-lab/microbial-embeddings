library(ape)
library(ggtree)
library(ggplot2)
library(ggtreeExtra)
library(stringr)
library(ggsci)
library(dplyr)
library(reshape2)
library(tidyverse)
library(ggnewscale)
library(RColorBrewer)
library(biomformat)
library(aplot)

tax <- read.delim("../Pretraining_data_profile/Data/taxmap_slv_ssu_ref_nr_138.2.txt",
  sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
)
acc <- paste(tax[[1]], tax[[2]], tax[[3]], sep = ".")
tax_split <- strsplit(tax$path, ";")
max_len <- max(sapply(tax_split, length))
tax_split <- lapply(tax_split, function(x) {
  length(x) <- max_len
  return(x)
})

otu_annotation <- as.data.frame(do.call(rbind, tax_split))
rownames(otu_annotation) <- acc
otu_annotation <- otu_annotation[, 1:7]
colnames(otu_annotation) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
otu_table <- read.csv("Data/OTU_prevalence_abundance.csv", stringsAsFactors = FALSE)
otu_annotation <- otu_annotation[otu_table$OTU, ]
tree <- read.tree("Data/tree.tre")

otu_annotation <- otu_annotation %>% rownames_to_column("OTU")
tippoint <- data.frame(
  OTU = tree$tip.label,
  stringsAsFactors = FALSE
) %>%
  left_join(otu_annotation %>% select(OTU, Phylum), by = "OTU") %>%
  mutate(
    OTU_factor = factor(OTU, levels = tree$tip.label),
    Taxa = factor(Phylum),
    taxon = if_else(
      Phylum %in% names(sort(table(Phylum), decreasing = TRUE)[1:12]),
      Phylum,
      "Others"
    )
  )

taxon_levels <- c(
  names(sort(table(tippoint$Phylum), decreasing = TRUE)[1:12]),
  "Others"
)
col <- c(brewer.pal(12, "Paired"), brewer.pal(4, "Pastel2"))

p <- ggtree(tree, layout = "rectangular", size = 0.5, branch.length = "none") +
  layout_dendrogram() +
  theme(
    axis.ticks.length = unit(0, "mm"),
    plot.margin = margin()
  )
p

Phylum_p <- ggplot(tippoint, aes(
  x = OTU_factor,
  y = "Phylum",
  fill = factor(taxon, levels = taxon_levels)
)) +
  geom_tile() +
  scale_fill_manual(values = setNames(col, taxon_levels)) +
  labs(
    fill = "Phylum",
    x = NULL,
    y = NULL
  ) +
  scale_y_discrete(position = "right") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.length = unit(0, "mm"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12, face = "italic"),
    plot.margin = margin()
  )

cor_colors <- colorRampPalette(c("white", "#FD763F"))(100)
tippoint$prev <- otu_table$prev[match(tippoint$OTU, otu_table$OTU)]

prev_p <- ggplot(tippoint, aes(x = OTU_factor, y = "log10(prev)", fill = log10(prev))) +
  geom_tile() +
  scale_fill_gradientn(colors = cor_colors) +
  labs(
    fill = "log10(prev)",
    x = NULL,
    y = NULL
  ) +
  scale_y_discrete(position = "right") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.length = unit(0, "mm"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.margin = margin()
  )
prev_p


# mean_abc
cor_colors <- colorRampPalette(c("white", "#23BAC5"))(100)
tippoint$mean_abc <- otu_table$mean_abc[match(tippoint$OTU, otu_table$OTU)]
mean_abc_p <- ggplot(tippoint, aes(x = OTU_factor, y = "log10(mean Abund.)", fill = log10(mean_abc))) +
  geom_tile() +
  scale_fill_gradientn(colors = cor_colors) +
  labs(
    fill = "log10(mean Abund.)",
    x = NULL,
    y = NULL
  ) +
  scale_y_discrete(position = "right") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.length = unit(0, "mm"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.margin = margin()
  )
mean_abc_p

embedding_table <- read.table("../../data/social_niche_embedding_100.txt", row.names = 1, quote = "\"", comment.char = "")
colnames(embedding_table) <- paste0("dim", 1:100)
embedding_table <- embedding_table[otu_table$OTU, ]
clust <- hclust(dist(t(embedding_table)), method = "ward.D2")
order_dim <- clust$labels[clust$order]
embedding_table <- embedding_table %>% rownames_to_column("fid")
embedding_table <- embedding_table %>% melt(id.vars = c("fid"), variable.name = "dim")
embedding_table$dim <- factor(embedding_table$dim, levels = order_dim)
col_fun <- colorRampPalette(c("#2574AA", "white", "#ED7B79"))(100)
embedding_table$fid <- factor(embedding_table$fid, levels = tree$tip.label)

embedding_p <- ggplot(embedding_table, aes(y = dim, x = fid, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = col_fun) +
  labs(
    fill = "SNE",
    x = "OTUs",
    y = NULL
  ) +
  scale_y_discrete(position = "right") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.length = unit(0, "mm"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.margin = margin()
  )
embedding_p

# 6
p_combined <- embedding_p %>%
  aplot::insert_top(mean_abc_p, height = 0.03) %>%
  aplot::insert_top(prev_p, height = 0.03) %>%
  aplot::insert_top(Phylum_p, height = 0.03) %>%
  aplot::insert_top(p, height = 0.5)
p_combined
ggsave("Figures/tree_otu_embedding_v2.png", p_combined, width = 15, height = 10, dpi = 300, units = "in", device = "png")

