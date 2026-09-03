library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  getwd()
}
repo_root <- normalizePath(file.path(script_dir, "../../.."))
table_dir <- file.path(script_dir, "results/tables")
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

read_coordinates <- function(name, group) {
  table <- read.csv(file.path(table_dir, name), check.names = FALSE)
  names(table)[1:3] <- c("OTU", "TSNE1", "TSNE2")
  table$Group <- group
  table
}

sne <- read_coordinates("t_sne_co.csv", "SNE")
phylo <- read_coordinates("t_sne_phylo.csv", "PhyloE")
coordinates <- rbind(sne, phylo)

for (rank in c("Phylum", "Class", "Order", "Family", "Genus")) {
  counts <- sort(table(taxa[sne$OTU, rank]), decreasing = TRUE)
  common_taxa <- setdiff(names(head(counts, 14)), "Incertae Sedis")
  labels <- ifelse(taxa[sne$OTU, rank] %in% common_taxa, taxa[sne$OTU, rank], "Others")
  coordinates$Taxon <- factor(rep(labels, 2), levels = c(common_taxa, "Others"))
  colors <- setNames(c(brewer.pal(12, "Paired"), brewer.pal(4, "Pastel2")), levels(coordinates$Taxon))
  colors["Others"] <- "grey"

  plot <- ggplot(coordinates, aes(TSNE1, TSNE2, color = Taxon)) +
    facet_wrap(~Group, scales = "free") +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = colors) +
    labs(x = "t-SNE 1", y = "t-SNE 2", color = rank) +
    theme_bw(base_size = 12) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  ggsave(file.path(figure_dir, paste0(rank, ".png")), plot, width = 25, height = 10, units = "cm")
}

traits <- read.delim(
  file.path(repo_root, "analysis/traits/data/traits_precalculated.txt"),
  row.names = 1, check.names = FALSE
)
trait_coordinates <- rbind(
  read_coordinates("t_sne_co_bugbase.csv", "SNE"),
  read_coordinates("t_sne_phy_bugbase.csv", "PhyloE")
)

plot_traits <- function(columns, labels, filename) {
  values <- traits[trait_coordinates$OTU, columns, drop = FALSE]
  names(values) <- labels
  plot_data <- cbind(trait_coordinates, values) |>
    melt(id.vars = c("OTU", "TSNE1", "TSNE2", "Group")) |>
    filter(value > 0)

  plot <- ggplot(plot_data, aes(TSNE1, TSNE2, color = variable)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_brewer(palette = "Set2") +
    facet_wrap(~Group, scales = "free") +
    labs(x = "t-SNE 1", y = "t-SNE 2", color = NULL) +
    theme_bw(base_size = 12) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  ggsave(file.path(figure_dir, filename), plot, width = 25, height = 10, units = "cm")
}

plot_traits(
  c("Aerobic", "Anaerobic", "Facultatively_Anaerobic"),
  c("Aerobic", "Anaerobic", "Facultatively Anaerobic"),
  "Aerobic_Anaerobic_Facultatively.png"
)
plot_traits(c("Gram_Negative", "Gram_Positive"), c("Gram Negative", "Gram Positive"), "Gram.png")
