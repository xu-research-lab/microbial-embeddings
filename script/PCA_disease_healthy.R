library(dplyr)
library(pvca)
library(Biobase)
library(biomformat)
library(cowplot)
library(FactoMineR)
library(factoextra)

methods <- c("russell_rao", "russell_rao_weight", "faith", "jaccard", "abundance-totalsum", 
             "abundance-percentile", "braycurtis-totalsum", "braycurtis-percentile", "phylogeny", "PCA")

labs<- c("Russell_rao", "Russell_rao_weight", "Faith", "Jaccard", "Abundance_totalsum", "Abundance_percentile",
         "Braycurtis_totalsum", "Braycurtis_percentile", "phylogeny_glove", "Phylogeny_PCA")

metadata <- read.table("/home/dongbiao/word_embedding_microbiome/programe_test/IBD/data/metadata.txt", sep="\t", 
                       header = TRUE, row.names = 1)
biom_table <- read_biom("/home/dongbiao/word_embedding_microbiome/programe_test/IBD/meta_analysis/table.biom")
table <- as.data.frame(biom_table$data)
colnames(table) <- rownames(biom_table)

embedding <- read.table("/home/dongbiao/word_embedding_microbiome/programe_test/AGP_run_20220702/russell_rao/russell_rao_100.txt", 
                        quote="\"", comment.char="", row.names = 1)
embedding <- embedding[rownames(embedding) != "<unk>", ]
keep_id <- intersect(rownames(embedding), colnames(table))
embedding <- embedding[keep_id, ]
table <- table[, keep_id]
table <- table / rowSums(table)

compute_pca_R2 <- function(table, metadata, group="group"){
    res <- list()
    res_pca <- PCA(table, scale.unit = TRUE, ncp = 2, graph = FALSE)
    res$eig_val <- get_eigenvalue(res_pca)
    res$pca_coordinates <- res_pca$ind$coord
    
    # caculate Pvca
    expr_set <- ExpressionSet(assayData = as.matrix(t(table_1)),
                              phenoData = AnnotatedDataFrame(metadata[rownames(table), ]))
    
    pct_threshold <- 0.6
    batch.factors <- c("study", "group")
    
    # raw data
    pvcaObj <- pvcaBatchAssess(expr_set, batch.factors, pct_threshold)
    res$weight_raw <- pvcaObj$dat
    
    return(res)
}


i <- 1
embedding_file <- "/home/dongbiao/word_embedding_microbiome/databaes/16Sdatabase/glove_1"
for(methd in methods){
    if (methd == "Relative_abundance"){
        res <- compute_pca_R2(table, metadata[rownames(table), ], group="group")
    }
    else {
        embedding <- read.table(sprintf("%s/%s/%s_100.txt", embedding_file, methd, methd), 
                                quote="\"", comment.char="", row.names = 1)
        embedding <- embedding[keep_id, ]
        table_1 <- as.matrix(table) %*% as.matrix(embedding) %>% as.data.frame()
        res <- compute_pca_R2(table_1, metadata[rownames(table_1), ], group="group")
    }
    plot_data <- res$pca_coordinates %>% as.data.frame()
    plot_data$group <- metadata[rownames(plot_data), "group"]
    plot_data$study <- metadata[rownames(plot_data), "study"]
    plot_data <- plot_data %>% mutate(status = if_else(group == 1, "IBD", "Control"))
    
    p1 <- plot_data %>% ggplot(aes(x = Dim.1, y = Dim.2, color = study, shape=status)) +
        geom_point(size=2) + theme_bw() +
        labs(x = paste0("PCA1", " ", round(res$eig_val[1, 2], 2), "%"),
             y = paste0("PCA2", " ", round(res$eig_val[2, 2], 2), "%"))
    
    plot_pvca <- as.matrix(res$weight_raw[1, c(2,3,4)]) %>% as.data.frame() %>% mutate(Effects = c("study", "status", "resid"))
    plot_pvca$Effects <- factor(plot_pvca$Effects, levels=c("status", "study", "resid"))
    p2 <- plot_pvca %>% ggplot(aes(x=Effects, y=V1)) + 
        geom_bar(stat = "identity", position = "dodge", width = 0.5) +
        theme_bw() +
        labs(y = "Weighted average proportion variance", fill = "Data type")
    
    p <- plot_grid(p1, p2, rel_widths = c(5, 2), labels = c('A', 'B'))
    title <- ggdraw() + draw_label(labs[i], fontface='bold')
    p <- plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
    ggsave(sprintf("/home/dongbiao/word_embedding_microbiome/programe_test/IBD/plot/%s_PCA_PVCA_all_feces.png", labs[i]),
           p, width = 8, height = 5)
    i = 1 + i
}

