library(plyr)
library(ggplot2)
library(dplyr)
library(purrr)
library(stringr)
library(pheatmap)
library(cowplot)
library(RColorBrewer)

methods <- c("russell_rao", "russell_rao_weight", "faith", "jaccard", "abundance-totalsum", 
             "abundance-percentile", "braycurtis-totalsum", "braycurtis-percentile", "phylogeny", "PCA")

labs<- c("Russell_rao", "Russell_rao_weight", "Faith", "Jaccard", "Abundance_totalsum", "Abundance_percentile",
         "Braycurtis_totalsum", "Braycurtis_percentile", "phylogeny_glove", "Phylogeny_PCA")

# 相关比(Correlation Ratio)
caculate_cor_ratio <- function(pathway_table, embedding_table){
    # filter pathway table
    total_col <- nrow(pathway_table)
    keep <- colSums(pathway_table) < round(total_col * 0.9, 0)
    keep2 <- colSums(pathway_table) > round(total_col * 0.1, 0)
    pathway_table <- pathway_table[, keep & keep2]
    # Match, clean, and center
    taxa_names <- intersect(rownames(pathway_table), rownames(embedding_table))
    pathway_table <- pathway_table[taxa_names, ]
    embedding_table <- embedding_table[taxa_names, ]
    # embedding_table <- apply(embedding_table, 2, function(x) return((x - mean(x) ) / sd(x))) 
    colnames(embedding_table) <- paste("property_100_", seq(1, 100), sep = "")
    
    # caculate the FoldChange and P-value
    # KEGG table是一个由0和1组成的矩阵
    sd_all <- apply(embedding_table, 2, sd)
    mean_all <- apply(embedding_table, 2, mean)
    cor_ratio <- list()
    for (idx in 1: ncol(pathway_table)){
        pathway_presence <- embedding_table[pathway_table[, idx] == 1, ]
        pathway_disappear <- embedding_table[pathway_table[, idx] == 0, ]
            
        cor_ratio[[idx]] <- sqrt(((colMeans(pathway_presence) - mean_all) ^2 * nrow(pathway_presence) + 
                          (colMeans(pathway_disappear) - mean_all) ^2 * nrow(pathway_disappear)) / (nrow(embedding_table) * sd_all^2))
        
    }
    cor_ratio <- data.frame(cor_ratio)
    colnames(cor_ratio) <- colnames(pathway_table)
    return(cor_ratio)
}

# import pathway table OTU * KEGG-pathway
pathway_table <- readRDS("/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/pathway/OTU_pathway_table.rds")
res_dir <- "/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/res/relation/"

cor_list <- list()
for (i in 1: 10){
    # import embedding vector
    file_dir <- sprintf("/home/dongbiao/word_embedding_microbiome/programe_test/AGP_run_20220702/%s/%s_100.txt", 
                        methods[i], methods[i])
    embedding_table <- read.table(file_dir, row.names = 1, header = FALSE)
    embedding_table <- embedding_table[rownames(embedding_table) != "<unk>", ]
    
    cor_list[[i]] <- caculate_cor_ratio(pathway_table, embedding_table)
    row_labels <- rep(" ", nrow(cor_list[[i]]))
    breaksList = seq(0, 1, by = .01)
    p <- pheatmap(cor_list[[i]], cluster_cols = FALSE, 
                  breaks = breaksList,
                  color = colorRampPalette(colors = c("white", "red"))(100),
                  treeheight_col = 0, treeheight_row = 0,
                  main = sprintf("%s", labs[i]), 
                  border_color = NA, labels_row = row_labels)
    ggsave(paste0(res_dir, labs[i], ".png"), p, width = 25, height = 10, dpi = 300)
}

# 随机生成的embedding向量的
embedding_table <- read.table("/home/dongbiao/word_embedding_microbiome/programe_test/AGP_run_20220702/russell_rao/russell_rao_100.txt", 
                              row.names = 1, header = FALSE)
embedding_table <- embedding_table[rownames(embedding_table) != "<unk>", ]
nrows <- nrow(embedding_table)
ncols <- ncol(embedding_table)
random_embed <- matrix(rnorm(nrows * ncols), nrow = nrows, ncol = ncols)
rownames(random_embed) <- rownames(embedding_table)
cor_list[[11]] <- caculate_cor_ratio(pathway_table, random_embed)
res_dir <- "/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/res/relation/"
row_labels <- rep(" ", nrow(cor_list[[11]]))
breaksList = seq(0, 1, by = .01)
p <- pheatmap(cor_list[[11]], cluster_cols = FALSE, 
              treeheight_col = 0, treeheight_row = 0,
              main = "Random", breaks = breaksList,
              color = colorRampPalette(colors = c("white", "red"))(100),
              border_color = NA, labels_row = row_labels)
ggsave(paste0(res_dir, "Random.png"), p, width = 25, height = 10, dpi = 300)


# get max cor value for each pathway
otu_pathway_name <- readRDS("/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/pathway/OTU_pathway_name.rds")
colnames(otu_pathway_name) <- data.frame(str_split(colnames(otu_pathway_name), "\\.", simplify = TRUE, n = 2))[, 1]
max_cor_value <- list()
max_cor_id <- list()
for (i in 1: 10){
    df <- t(data.frame(cor_value = apply(cor_list[[i]], 2, function(x) max(abs(x)))))
    id <- t(data.frame(cor_value = apply(cor_list[[i]], 2, function(x) which(abs(x) == max(abs(x))))))
    max_cor_value[[i]] <- as.data.frame(df)
    max_cor_id[[i]] <- as.data.frame(id)
}
max_cor_value <- as.data.frame(t(rbind.fill(max_cor_value)))
id <- intersect(colnames(otu_pathway_name), rownames(max_cor_value))
max_cor_value <- max_cor_value[id, ]
colnames(max_cor_value) <- c("Russell_rao", "Russell_rao_weight", "Faith", "Jaccard", "Abundance_totalsum", "Abundance_percentile",
                             "Braycurtis_totalsum", "Braycurtis_percentile", "Phylogeny_glove", "Phylogeny_PCA")
rowname_table <- rownames(max_cor_value)


max_cor_id <- as.data.frame(t(rbind.fill(max_cor_id)))
id <- intersect(colnames(otu_pathway_name), rownames(max_cor_id))
max_cor_id <- max_cor_id[id, ]
colnames(max_cor_id) <- c("Russell_rao", "Russell_rao_weight", "Faith", "Jaccard", "Abundance_totalsum", "Abundance_percentile",
                          "Braycurtis_totalsum", "Braycurtis_percentile", "Phylogeny_glove", "Phylogeny_PCA")

heigh_cor_lable <- list()
for (i in 1:10){
    heigh_cor_lable[[i]] = rownames(max_cor_value)[order(max_cor_value[[i]], decreasing = TRUE)]
}

plot_cor_heatmap <- function(max_cor_value, max_cor_name, file_name, width, height){
    label_table <- t(otu_pathway_name[, rownames(max_cor_value)])
    label_1 <- label_table[max_cor_name, ] %>% as.data.frame() %>% arrange(desc(V2))
    colnames(label_1) <- c("Name", "Class")
    breaksList = seq(0, 1, by = .01)
    
    p1 <- pheatmap(max_cor_value[max_cor_name, 1:10][rownames(label_1), ],
                   treeheight_col = 0, treeheight_row = 0,
                   color = colorRampPalette(colors = c("white", "red"))(100),
                   breaks = breaksList, cluster_row = FALSE,
                   labels_row = label_1[, 1], annotation_row = label_1 %>% dplyr::select(Class))
    plot_map_OTU <- t(pathway_table[, rownames(max_cor_value)])
    p2 <- pheatmap(plot_map_OTU[max_cor_name, ][rownames(label_1), ], cluster_row = FALSE, 
                   treeheight_col = 0, treeheight_row = 0, legend = FALSE,
                   show_rownames = FALSE, show_colnames = FALSE)
    p2_1 <- plot_grid(p2[[4]], NULL, rel_heights = c(8, 1.4), ncol = 1)
    p <- plot_grid(p2_1, p1[[4]], rel_widths = c(1, 9))
    ggsave(sprintf("/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/pathway/%s.png", file_name), p, 
           width = width, height = height, units = 'cm', dpi = 300)
}
max_cor_name <- Reduce(union, list(heigh_cor_lable[[1]][1:25], heigh_cor_lable[[2]][1:25], heigh_cor_lable[[3]][1:25], 
                                   heigh_cor_lable[[4]][1:25], heigh_cor_lable[[5]][1:25], heigh_cor_lable[[6]][1:25],
                                   heigh_cor_lable[[7]][1:25], heigh_cor_lable[[8]][1:25], heigh_cor_lable[[9]][1:25],
                                   heigh_cor_lable[[10]][1:25]))
plot_cor_heatmap(max_cor_value, max_cor_name, "max_cor_pathway", width = 36, height = 25)
n <- nrow(max_cor_value)
min_cor_name <- Reduce(union, list(heigh_cor_lable[[1]][(n-24) : n], heigh_cor_lable[[2]][(n-24) : n], heigh_cor_lable[[3]][(n-24) : n], 
                                   heigh_cor_lable[[4]][(n-24) : n], heigh_cor_lable[[5]][(n-24) : n], heigh_cor_lable[[6]][(n-24) : n],
                                   heigh_cor_lable[[7]][(n-24) : n], heigh_cor_lable[[8]][(n-24) : n], heigh_cor_lable[[9]][(n-24) : n],
                                   heigh_cor_lable[[10]][n-24]))
plot_cor_heatmap(max_cor_value, 
                 min_cor_name, 
                 "lowest_cor_pathway", width = 40, height = 25)


# Get top 5 pathways for each method of 
methods <- c("russell_rao", "russell_rao_weight", "faith", "jaccard", "abundance-totalsum", 
             "abundance-percentile", "braycurtis-totalsum", "braycurtis-percentile", "phylogeny", "PCA")
labs<- c("Russell_rao", "Russell_rao_weight", "Faith", "Jaccard", "Abundance_totalsum", "Abundance_percentile",
         "Braycurtis_totalsum", "Braycurtis_percentile", "phylogeny_glove", "Phylogeny_PCA")
i=1
for (method in methods){
    glove_emb <- read.table(paste("/home/dongbiao/word_embedding_microbiome/programe_test/AGP_run_20220702/", 
                                  method, "/", method, "_100.txt", 
                                  sep = ""),
                            quote="\"", comment.char="", row.names = 1, sep = " ", header = F)
    glove_emb = glove_emb[-which(rownames(glove_emb) == '<unk>'), ]
    colnames(glove_emb) <- paste("property_100_", seq(1, 100), sep = "")
    embed_table_glove <- glove_emb
    embed_table_glove <- embed_table_glove[taxa_names, ]
    embed_table_glove <- apply(embed_table_glove, 2, function(x) return((x - mean(x) ) / sd(x)))
    
    kegg_name <- otu_pathway_name[, heigh_cor_lable[[i]][1:5]]
    table_1 <- pathway_table[, heigh_cor_lable[[i]][1:5]]
    table_2 <- embed_table_glove[,max_cor_id[heigh_cor_lable[[i]][1:5], i]]
    table_2 <- table_2[rownames(table_1), ]
    x <- c(table_1[, 1], table_1[,2], table_1[,3], table_1[, 4], table_1[, 5])
    y <- c(table_2[, 1], table_2[,2], table_2[,3], table_2[, 4], table_2[, 5])
    Group <- c(rep(paste(kegg_name[1, 1], colnames(table_2)[1], sep = "-"), length(table_1[, 1])),
               rep(paste(kegg_name[1, 2], colnames(table_2)[2], sep = "-"), length(table_1[, 1])),
               rep(paste(kegg_name[1, 3], colnames(table_2)[3], sep = "-"), length(table_1[, 1])),
               rep(paste(kegg_name[1, 4], colnames(table_2)[4], sep = "-"), length(table_1[, 1])),
               rep(paste(kegg_name[1, 5], colnames(table_2)[5], sep = "-"), length(table_1[, 1])))
    Group <- c(Group, Group)
    
    box_plot <- data.frame(x = x, y = y, Group = Group)
    box_plot$x <- as.factor(box_plot$x)
    p <- box_plot %>% ggplot(aes(x, y, group=x)) + geom_boxplot() +
        labs(x="pathway", y="property", title = paste(labs[i])) +
        facet_wrap(~ Group, ncol = 5, scales = "free")
    ggsave(paste("/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/pathway/", method, "_boxplots.png", sep = ""), p,
           width = 50, height = 25, units = 'cm', dpi = 300)
    i = i + 1

}







