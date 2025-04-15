library(tidyverse)
library(cowplot)
library(RColorBrewer)

### t-SNE taxonomy
col <- c(brewer.pal(12,'Paired'), brewer.pal(4,'Pastel2'))
df_tsne_tax_co <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/res/tax_plot/tsne_phylum_co.csv")
df_tsne_tax_co$group <- rep("Co_embedding", n=nrow(df_tsne_tax_co))
df_tsne_tax_phy <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/res/tax_plot/tsne_phylum_phy.csv")
df_tsne_tax_phy$group <- rep("Phy_embedding", n=nrow(df_tsne_tax_phy))
df_tsne_tax <- rbind(df_tsne_tax_co, df_tsne_tax_phy)
p1 <- df_tsne_tax %>% ggplot(aes(x = x, y = y, color = pick_tax))+
    facet_wrap(~group, scales = "free")+
    geom_point(alpha = 0.5, size = 1) +
    theme_bw(base_size=12) +
    scale_color_manual(values = col) +
    labs(x = "t-SNE 1", y = "t-SNE 2", color = "Phylum") +
    theme(
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))
ggsave("/home/dongbiao/word_embedding_microbiome/result/supplement_igure_7.png", p1,
       width = 24, height = 12, units = "cm")

### t-SNE-bugbase
df_bugbase_co <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/res/tax_plot/bugbase_Facultatively_Anaerobic_co.csv")
df_bugbase_co$Type[df_bugbase_co$Type == "Facultatively_Anaerobic"] <- "Facultatively \n Anaerobic"
df_bugbase_co$group <- rep("Co_embedding", n=nrow(df_bugbase_co))
df_bugbase_phy <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/res/tax_plot/bugbase_Facultatively_Anaerobic_phy.csv")
df_bugbase_phy$Type[df_bugbase_phy$Type == "Facultatively_Anaerobic"] <- "Facultatively \n Anaerobic"
df_bugbase_phy$group <- rep("Phy_embedding", n=nrow(df_bugbase_phy))
df_bugbase <-rbind(df_bugbase_co, df_bugbase_phy)
p2 <- df_bugbase %>% ggplot(aes(x = x, y = y, color = Type))+
    geom_point(alpha = 0.5, size = 1) +
    facet_wrap(~group, scales = "free")+
    theme_bw(base_size=12) +
    scale_color_manual(values = col) +
    labs(x = "t-SNE 1", y = "t-SNE 2", fill = "Type") +
    theme(
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))
ggsave("/home/dongbiao/word_embedding_microbiome/result/supplement_igure_8.png", p2,
       width = 24, height = 12, units = "cm")

### King - men + women = queen
df_cos  <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/glove/order_Facultatively_Anaerobic_Anaerobic_Aerobic.csv")

p3 <- df_cos %>% ggplot(aes(x = cosine)) +
    geom_density(aes(fill = Type), alpha = 0.5)+
    # facet_grid(. ~ Type, scales="free") + 
    theme_bw(base_size=12) +
    labs(x = "Cosine", y = "Density", fill = "") +
    theme_classic(base_size=12) +
    theme(
        legend.position = "top") 

### 拼图
p <- plot_grid(p1, p2, p3, NULL,
               labels = c('a', 'b', 'c', ''),
               align="hv", rel_heights = c(2, 1),
               rel_widths = c(1.2, 1),
               nrow = 2, ncol =2, plot=FALSE)

ggsave("/home/dongbiao/word_embedding_microbiome/result/Figure_4_embedding.png", p,
       width = 35, height = 20, units = "cm")
