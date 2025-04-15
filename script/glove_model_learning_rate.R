library(ggplot2)
library(readxl)
library(cowplot)
library(dplyr)
library(stringr)
# file_dir <- "/home/dongbiao/word_embedding_microbiome/programe_test/AGP_run_20220702"
file_dir <- "/home/dongbiao/word_embedding_microbiome/databaes/16Sdatabase/glove_1"
res_list <- list()

for (i in 1:8){
    table <- read.table(paste0(file_dir, "/", "glove_", i, ".err"), sep = "\t", skip = 15)
    res_list[[i]] <- as.numeric(str_split(table$V1, ":", simplify = TRUE)[,4])
}


#创建dataframe，收录不同的matrix对用的lr和X_max
lr_xmax = data.frame(lr = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.02, 0.05, 0.05),
                     x_max = c(0.1, 0.01, 0.2, 0.04, 0.0001, 0.05, 0.1, 0.1))
rownames(lr_xmax) <- c("Russell_rao", "Russell_rao_weight", "Jaccard", "Faith", 
                       "Abundance_totalsum", "Abundance_percentile", 
                       "Braycurtis_totalsum", "Braycurtis_percentile")
plot_list <- list()
for (i in c(1: 8)){
    table_loss <- data.frame(epochs = 1:100, loss = res_list[[i]])
    plot_list[[i]] <-  table_loss %>% ggplot(aes(x = epochs, y = loss)) +
                     geom_point(size = 0.2) + 
                     geom_line() + 
                     theme_bw() + 
                     labs(title = sprintf("%s: lr=%s, X_max=%s", 
                                          rownames(lr_xmax)[i], 
                                          lr_xmax[i, 1],
                                          lr_xmax[i, 2]), size = 24) +
                    theme(axis.text = element_text(size = 24),
                          axis.title = element_text(size = 24),
                          plot.title = element_text(size = 24))
}

p <- plot_grid(plot_list[[1]], plot_list[[2]], 
               plot_list[[3]], plot_list[[4]],
               plot_list[[5]], plot_list[[6]],
               plot_list[[7]], plot_list[[8]],
               ncol = 4)
ggsave("/home/dongbiao/word_embedding_microbiome/programe_test/glove/loss_all_feces.png", 
       p, width = 90, height = 30, units = 'cm', dpi = 300)

