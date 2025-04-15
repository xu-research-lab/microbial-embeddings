library(msa)
library(phangorn)
library(ape)
library(ggplot2)
library(cowplot)
library(lemon)
library(scales)
library(gridExtra)
library(phyloseq)
library(dplyr)


args <- commandArgs(trailingOnly = TRUE)
arg1 <- args[1] # the dir for embedding file
arg2 <- args[2] # the dir for results file
arg3 <- args[3] # the type of embedding file training datasets


#1. write function to build tree (hierarchical clustering) from embeddings
hierarchical_cluster <- function(mat){
    dist_mat <- dist(mat, method = 'euclidean')
    hclust_avg <- hclust(dist_mat, method = 'average')
    return(hclust_avg)
}

# tree_file <- "/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/exported-tree_1/tree.nwk"
# phy_tree = read.tree(tree_file)
#4. write function to permute rownames on matrix
# null_dists <- list()
# num_iter <- 100
# phy_tree_permute <- phy_tree
# random_embedding <- data.frame(matrix(nrow = length(phy_tree$tip.label), ncol = 100))
# rownames(random_embedding) <- phy_tree$tip.label
# for(i in seq(1, num_iter)){
#     print(i)
#     for (j in 1 : length(phy_tree$tip.label)){
#         random_embedding[j, ] <- runif(100, -1, 1)
#     }
#     hclust_random <- hierarchical_cluster(random_embedding)
#     random_tree <- as.phylo(hclust_random)
#     dist_perm <- phangorn::treedist(random_tree, phy_tree)
#     null_dists[[i]] <- dist_perm
# }
# print("calculated null distances")
# saveRDS(null_dists, sprintf("/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/null_dists.rds",))
null_dists <- readRDS("/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/null_dists.rds")
methods <- c("russell_rao_weight", "russell_rao", "faith", "jaccard", "abundance-percentile", "abundance-totalsum",
             "braycurtis-totalsum", "braycurtis-percentile", "phylogeny", "PCA")
for (method in methods){
    embedding_file = sprintf("%s/%s/%s_100.txt", arg1, method, method)
    embedding <- read.table(embedding_file, quote="\"", comment.char="", row.names = 1)
    embedding <- embedding[rownames(embedding) != "<unk>", ]

    #2. read phylogenetic tree build using fasttree and clustalo msa
    tree_file <- "/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/exported-tree_1/tree.nwk"
    phy_tree = read.tree(tree_file)
    print("read phy tree")
    embedding = embedding[rownames(embedding) %in% phy_tree$tip.label, ]

    #1. write function to build tree (hierarchical clustering) from embeddings
    hclust_true <- hierarchical_cluster(embedding)
    hclust_tree <- as.phylo(hclust_true)
    print("Tree built from embeddings")
    #保存结果
    saveRDS(hclust_tree, sprintf("%s/hclust_embed_%s_%s.rds", arg2, arg3, method))

    phy_tree <- drop.tip(phy_tree, phy_tree$tip.label[!phy_tree$tip.label %in% hclust_tree$tip.label])
    #3. write function to calculate distance between two trees
    dist_true <- phangorn::treedist(hclust_tree, phy_tree)
    #print("Calculated true distance")
    saveRDS(dist_true, sprintf("%s/dist_true_%s_%s.rds", arg2, arg3, method))
    #print("Calculated distance")
}

    #4. write function to permute rownames on matrix


    #Checks the distance between two trees in comparison to random trees
    #####################################################
    ########### load objects and check results ##########
    #####################################################
## 画图时候的缩放函数
squash_axis <- function(from, to, factor) {
    # Args:
    #   from: left end of the axis
    #   to: right end of the axis
    #   factor: the compression factor of the range [from, to]

    trans <- function(x) {
        # get indices for the relevant regions
        isq <- x > from & x < to
        ito <- x >= to

        # apply transformation
        x[isq] <- from + (x[isq] - from)/factor
        x[ito] <- from + (to - from)/factor + (x[ito] - to)

        return(x)
    }

    inv <- function(x) {
        # get indices for the relevant regions
        isq <- x > from & x < from + (to - from)/factor
        ito <- x >= from + (to - from)/factor

        # apply transformation
        x[isq] <- from + (x[isq] - from) * factor
        x[ito] <- to + (x[ito] - (from + (to - from)/factor))

        return(x)
    }

    # return the transformation
    return(trans_new("squash_axis", trans, inv))
}

matrix_name <- c("Russell_rao_weight", "Russell_rao", "Faith", "Jaccard",
                 "Abundance_percentile", "Abundance_totalsum",
                 "Braycurtis_totalsum", "Braycurtis_percentile", 
                 "Phylogeny_glove", "Phylogeny_PCA")

i = 1
dist_true <- c()
for (method in methods){
    temp <- readRDS(sprintf("%s/dist_true_%s_%s.rds", arg2, arg3, method))
    dist_true <- c(dist_true, temp[i])
}

vline_data <- data.frame(value = dist_true, group = matrix_name) %>% arrange(value)
grou_name <- vline_data$group
vline_data$group <- factor(vline_data$group, levels = grou_name)
sym_diffs_total <- c()
for(null_dist_list in null_dists){
    sym_diffs_total <- c(sym_diffs_total, null_dist_list[i])
}
color_plot <- c("Russell_rao_weight"="#B7B5A0", "Russell_rao"="#447570", "Faith"="#452A3D", "Jaccard"="#D44C3C",
                "Abundance_percentile"="#E73847", "Abundance_totalsum"="#A8DADB",
                "Braycurtis_totalsum"="#457B9D", "Braycurtis_percentile"="#1D3557", 
                "Phylogeny_glove"="#EED5B7", "Phylogeny_PCA"="#4EAB90")
p1 <- ggplot(data.frame(sym_diffs_total))+
        geom_histogram(aes(x = sym_diffs_total), fill = "lightblue", binwidth = 1, center=1)+
        geom_vline(data = vline_data, aes(xintercept = value, colour = group), linetype = 4, size = 2)+ theme_bw()+
        theme(text = element_text(size=15)) +
        scale_color_manual(values = color_plot) +
        scale_x_continuous(breaks = c(18500, 20270, 20720, 20750, 20800, 20850, 20900)) +
        coord_trans(x = squash_axis(18500, 20270, 50)) +
        labs(x="Symmetric Difference", y="", color="")


i = 2
dist_true <- c()
for (method in methods){
    temp <- readRDS(sprintf("%s/dist_true_%s_%s.rds", arg2, arg3, method))
    dist_true <- c(dist_true, temp[i])
}
vline_data <- data.frame(value = dist_true, group = matrix_name) %>% arrange(value)
vline_data$group <- factor(vline_data$group, levels = grou_name)
sym_diffs_total <- c()
for(null_dist_list in null_dists){
    sym_diffs_total <- c(sym_diffs_total, null_dist_list[i])
}

p2 <- ggplot(data.frame(sym_diffs_total))+
        geom_histogram(aes(x = sym_diffs_total), fill = "lightblue", binwidth = 0.1, center=1) +
        geom_vline(data = vline_data, aes(xintercept = value, colour = group), linetype = 4, size = 2)+ theme_bw()+
        theme(text = element_text(size=15)) +
        scale_color_manual(values = color_plot) +
        scale_x_continuous(breaks = c(5, 30, 90, 150, 250, 330, 340)) +
        coord_trans(x = squash_axis(100, 330, 20)) +
        labs(x="Branch Score Difference", y="", color="")

ggsave(sprintf("%s/sym_diffs_%s_1.png", arg2, arg3),
       p1, dpi = 300, width = 35, height = 20, units = 'cm')

ggsave(sprintf("%s/sym_diffs_%s_2.png", arg2, arg3),
       p2, dpi = 300, width = 35, height = 20, units = 'cm')

# p3 <- grid_arrange_shared_legend(p1, p2, ncol = 2,
#                                  position='bottom')
# ggsave("/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/res/dis_caculate.png",
#        p3, dpi = 300, width = 60, height = 20, units = 'cm')

