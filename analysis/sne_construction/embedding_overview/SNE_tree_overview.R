library(ape)
library(aplot)
library(dplyr)
library(ggplot2)
library(ggtree)
library(RColorBrewer)
library(reshape2)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  getwd()
}
repo_root <- normalizePath(file.path(script_dir, "../../.."))
overview_data <- file.path(repo_root, "analysis/sne_construction/data/embedding_overview")
figure_dir <- file.path(script_dir, "results/figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

taxmap <- read.delim(
  file.path(repo_root, "data/taxmap_slv_ssu_ref_nr_138.2.txt"),
  sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
)
accession <- paste(taxmap[[1]], taxmap[[2]], taxmap[[3]], sep = ".")
taxa <- strsplit(taxmap$path, ";")
max_rank <- max(lengths(taxa))
taxa <- lapply(taxa, function(x) `length<-`(x, max_rank))
taxa <- as.data.frame(do.call(rbind, taxa))
rownames(taxa) <- accession
taxa <- taxa[, 1:7]
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

otu <- read.csv(file.path(overview_data, "OTU_prevalence_abundance.csv"))
tree <- read.tree(file.path(repo_root, "data/SSURefNR99_1200_slv_138_2_subset.tre"))
tip_data <- data.frame(OTU = tree$tip.label) %>%
  left_join(data.frame(OTU = rownames(taxa), Phylum = taxa$Phylum), by = "OTU") %>%
  mutate(
    OTU = factor(OTU, levels = tree$tip.label),
    Taxon = if_else(
      Phylum %in% names(sort(table(Phylum), decreasing = TRUE)[1:12]),
      Phylum,
      "Others"
    )
  )

taxon_levels <- c(names(sort(table(tip_data$Phylum), decreasing = TRUE)[1:12]), "Others")
taxon_colors <- c(brewer.pal(12, "Paired"), brewer.pal(4, "Pastel2"))

tree_plot <- ggtree(tree, layout = "rectangular", linewidth = 0.5, branch.length = "none") +
  layout_dendrogram() +
  theme(axis.ticks.length = grid::unit(0, "mm"), plot.margin = margin())

phylum_plot <- ggplot(tip_data, aes(OTU, "Phylum", fill = factor(Taxon, levels = taxon_levels))) +
  geom_tile() +
  scale_fill_manual(values = setNames(taxon_colors, taxon_levels)) +
  labs(fill = "Phylum", x = NULL, y = NULL) +
  scale_y_discrete(position = "right") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin())

tip_data$prevalence <- otu$prevalence[match(tip_data$OTU, otu$OTU)]
prevalence_plot <- ggplot(tip_data, aes(OTU, "log10(prev)", fill = log10(prevalence))) +
  geom_tile() +
  scale_fill_gradientn(colors = colorRampPalette(c("white", "#FD763F"))(100)) +
  labs(fill = "log10(prev)", x = NULL, y = NULL) +
  scale_y_discrete(position = "right") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin())

tip_data$mean_abundance <- otu$mean_abc[match(tip_data$OTU, otu$OTU)]
abundance_plot <- ggplot(tip_data, aes(OTU, "log10(mean Abund.)", fill = log10(mean_abundance))) +
  geom_tile() +
  scale_fill_gradientn(colors = colorRampPalette(c("white", "#23BAC5"))(100)) +
  labs(fill = "log10(mean Abund.)", x = NULL, y = NULL) +
  scale_y_discrete(position = "right") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin())

embedding <- read.table(
  file.path(repo_root, "data/social_niche_embedding_100.txt"),
  row.names = 1, quote = "\"", comment.char = ""
)
colnames(embedding) <- paste0("dim", seq_len(ncol(embedding)))
embedding <- embedding[otu$OTU, , drop = FALSE]
dimension_clustering <- hclust(dist(t(embedding)), method = "ward.D2")
dimension_order <- dimension_clustering$labels[dimension_clustering$order]
embedding$OTU <- rownames(embedding)
embedding <- melt(embedding, id.vars = "OTU", variable.name = "Dimension")
embedding$OTU <- factor(embedding$OTU, levels = tree$tip.label)
embedding$Dimension <- factor(embedding$Dimension, levels = dimension_order)

embedding_plot <- ggplot(embedding, aes(OTU, Dimension, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = colorRampPalette(c("#2574AA", "white", "#ED7B79"))(100)) +
  labs(fill = "SNE", x = "OTUs", y = NULL) +
  scale_y_discrete(position = "right") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), plot.margin = margin())

combined <- embedding_plot %>%
  insert_top(abundance_plot, height = 0.03) %>%
  insert_top(prevalence_plot, height = 0.03) %>%
  insert_top(phylum_plot, height = 0.03) %>%
  insert_top(tree_plot, height = 0.5)

ggsave(
  file.path(figure_dir, "tree_otu_embedding.png"),
  combined, width = 15, height = 10, dpi = 300, units = "in"
)
