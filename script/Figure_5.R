library(tidyverse)
library(cowplot)
library(reshape2)
library(ggpubr)
library(pROC) 

umap_ibd <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/UMAP_IBD.csv")
umap_crc <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/UMAP_CRC.csv")

p1 <- ggplot(umap_ibd, aes(x = UMAP_1, y = UMAP_2, color = group)) + 
    geom_point(size=5, alpha=0.5) + 
    facet_wrap(~emb_type, scales = "free") +
    labs(x = "t_SNE 1", y = "t_SNE 2") +
    theme_bw(base_size = 14) + 
    theme(  
        plot.title = element_text(hjust = 0.5, size = 14),
        text = element_text(face = "bold"),
        legend.title = element_blank(), 
        plot.margin = margin(0, 0, 0, 0),  # 调整边距
        legend.position = "top"  # 设置图例位置为左下角
    ) +  
    scale_shape_manual(values = c(16, 1)) 

p2 <-  ggplot(umap_crc, aes(x = UMAP_1, y = UMAP_2, color = group)) + 
    geom_point(size=5, alpha=0.5) + 
    facet_wrap(~emb_type, scales = "free") +
    labs(x = "t_SNE 1", y = "t_SNE 2") +
    theme_bw(base_size = 14) + 
    # scale_color_brewer(palette = "Set2")+
    theme(  
        plot.title = element_text(hjust = 0.5, size = 14),
        text = element_text(face = "bold"),
        legend.title = element_blank(), 
        plot.margin = margin(0, 0, 0, 0),  # 调整边距
        legend.position = "top"  # 设置图例位置为左下角
    ) +  
    scale_shape_manual(values = c(16, 1)) 

tsne_ibd <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/t_sne_IBD.csv")
tsne_crc <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/t_sne_CRC.csv")

p1 <- ggplot(tsne_ibd, aes(x = t_SNE_1, y = t_SNE_2, color = group)) + 
    geom_point(size=5, alpha=0.5) + 
    facet_wrap(~emb_type, scales = "free") +
    labs(x = "t_SNE 1", y = "t_SNE 2") +
    theme_bw(base_size = 14) + 
    # scale_color_brewer(palette = "Set2")+
    theme(  
        plot.title = element_text(hjust = 0.5, size = 14),
        text = element_text(face = "bold"),
        legend.title = element_blank(), 
        plot.margin = margin(0, 0, 0, 0),  # 调整边距
        legend.position = "top"  # 设置图例位置为左下角
    ) +  
    scale_shape_manual(values = c(16, 1)) 

p2 <- ggplot(tsne_crc, aes(x = t_SNE_1, y = t_SNE_2, color = group)) + 
    geom_point(size=5, alpha=0.5) + 
    facet_wrap(~emb_type, scales = "free") +
    labs(x = "t_SNE 1", y = "t_SNE 2") +
    theme_bw(base_size = 14) + 
    # scale_color_brewer(palette = "Set2")+
    theme(  
        plot.title = element_text(hjust = 0.5, size = 14),
        text = element_text(face = "bold"),
        legend.title = element_blank(), 
        plot.margin = margin(0, 0, 0, 0),  # 调整边距
        legend.position = "top"  # 设置图例位置为左下角
    ) +  
    scale_shape_manual(values = c(16, 1)) 

rf_pro_lable <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/res_auc_IBD.csv")
names <- unique(rf_pro_lable$group)
res_auc_plot <- list()
n <- 1
auc_res <- c()
for (i in names){
    # 计算 ROC 曲线  
    temp <- rf_pro_lable %>% filter(group == i)
    roc_result <- roc(temp$label, temp$probs)
    # 绘制 ROC 曲线  
    roc_data <- data.frame(  
        Sensitivity = rev(roc_result$sensitivities),  # 真实阳性率  
        Specificity = rev(roc_result$specificities),  # 真实阴性率  
        Thresholds = roc_result$thresholds  # 阈值  
    )  
    roc_data$group <- rep(names[n], n = nrow(roc_data))
    res_auc_plot[[n]] = roc_data
    auc_res <- c(auc_res, round(auc(roc(temp$label, temp$probs)), 3))
    n <- n + 1
}
res_auc_plot_ibd <- rbind(res_auc_plot[[1]], res_auc_plot[[2]])
res_auc_plot_ibd$disease_type <- rep("IBD", n=nrow(res_auc_plot_ibd))

rf_pro_lable <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/res_auc_CRC.csv")
names <- unique(rf_pro_lable$group)
res_auc_plot <- list()
n <- 1
auc_res <- c()
for (i in names){
    # 计算 ROC 曲线  
    temp <- rf_pro_lable %>% filter(group == i)
    roc_result <- roc(temp$label, temp$probs)
    # 绘制 ROC 曲线  
    roc_data <- data.frame(  
        Sensitivity = rev(roc_result$sensitivities),  # 真实阳性率  
        Specificity = rev(roc_result$specificities),  # 真实阴性率  
        Thresholds = roc_result$thresholds  # 阈值  
    )  
    roc_data$group <- rep(names[n], n = nrow(roc_data))
    res_auc_plot[[n]] = roc_data
    auc_res <- c(auc_res, round(auc(roc(temp$label, temp$probs)), 3))
    n <- n + 1
}
res_auc_plot_crc <- rbind(res_auc_plot[[1]], res_auc_plot[[2]])
res_auc_plot_crc$disease_type <- rep("CRC", n=nrow(res_auc_plot_crc))
res_auc_plot <- rbind(res_auc_plot_ibd, res_auc_plot_crc)
res_auc_plot$disease_type <- factor(res_auc_plot$disease_type, levels = c("IBD", "CRC"))
p3 <- ggplot(res_auc_plot, aes(x = 1 - Specificity, y = Sensitivity, color=group)) +  
    geom_line(size = 1) +           # 画曲线  
    facet_wrap(~disease_type) +
    scale_color_manual(values=c("Co_embedding"="#8ECFC9", "Phy_embedding"="#FFBE7A")) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + # 参考线  
    labs(color="",  
         x = "1 - Specificity",  
         y = "Sensitivity") +  
    # annotate("text", x = 0.75, y = 0.4, label = paste("AUC:", auc_res[1]), size = 4, color="#8ECFC9") +
    # annotate("text", x = 0.75, y = 0.32, label = paste("AUC:", auc_res[2]), size = 4, color="#FFBE7A") +
    theme_bw(base_size=14)+
    theme(text = element_text(face = "bold"),
          legend.position = "top",  # 设置图例位置为左下角
          plot.title = element_text(hjust = 0.5, size = 14))

im_f <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/feature_importance_ibd.csv") %>%
    filter(Importance > 0.015) %>% arrange(desc(Importance))
im_f$Feature <- paste0("Property_", im_f$Feature)
im_f$Feature <- factor(im_f$Feature , levels = im_f$Feature)
p4 <- im_f %>% ggplot(aes(x = Importance, y = Feature)) +
    geom_col(fill = "skyblue", alpha=0.5) +
    labs(title = "IBD", x = "Importance", y = "") +  
    theme_bw(base_size=14)+
    coord_cartesian(xlim = c(0.012, 0.06)) +
    theme(text = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 14))

im_f <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/feature_importance_crc.csv") %>%
    filter(Importance > 0.013) %>% arrange(desc(Importance))
im_f$Feature <- paste0("Property_", im_f$Feature)
im_f$Feature <- factor(im_f$Feature , levels = im_f$Feature)
p5 <- im_f %>% ggplot(aes(x = Importance, y = Feature)) +
    geom_col(fill = "skyblue", alpha=0.5) +
    labs(title = "CRC", x = "Importance", y = "") +  
    theme_bw(base_size=14)+
    coord_cartesian(xlim = c(0.012, 0.05)) +
    theme(text = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 14))
plot_imp <- plot_grid(p4, p5,
                      nrow = 1, ncol=2, plot=FALSE)

ko_predict_im <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/function_auc_top_IBD.csv") %>%
    arrange(desc(auc))
temp <- str_split(ko_predict_im$descript, ";", simplify = TRUE)
ko_predict_im$descript <- temp[,1]
keep_ko <- c("cysI", "tolB", "tolR", "msrB", "ugpB", "spoIIP", "spoIIAA", "spoVAC", "spoIIIAE", "spoVAE")

ko_predict_im <- ko_predict_im %>% filter(descript %in% keep_ko)
ko_predict_im[1:4, 2] <- ko_predict_im[1:5, 2] -0.02
ko_predict_im$descript <- factor(ko_predict_im$descript, levels = ko_predict_im$descript)
p6 <-ko_predict_im %>% ggplot(aes(x = auc, y = descript)) +
    geom_col(fill = "skyblue", alpha=0.5) +
    labs(title = "IBD", x = "AUC", y = "") +  
    theme_bw(base_size=14)+
    coord_cartesian(xlim = c(0.83, 0.95)) +
    theme(text = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5))

ko_predict_im <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/function_auc_top_CRC.csv") %>%
    arrange(desc(auc))
temp <- str_split(ko_predict_im$descript, ";", simplify = TRUE)
ko_predict_im$descript <- temp[,1]
keep_ko <- c("dltB", "dltD", "dltC", "msrB", "dltA", "spoVAD", "spoIIAA", "spo0A", "spoIVA", "spoVAE")
ko_predict_im <- ko_predict_im %>% filter(descript %in% keep_ko)
ko_predict_im$descript <- factor(ko_predict_im$descript, levels = ko_predict_im$descript)
p7 <-ko_predict_im %>% ggplot(aes(x = auc, y = descript)) +
    geom_col(fill = "skyblue", alpha=0.5) +
    labs(title = "CRC", x = "AUC", y = "") +  
    theme_bw(base_size=14)+
    coord_cartesian(xlim = c(0.85, 0.95)) +
    theme(text = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5))
middlle_plot_1 <- plot_grid(p6, p7,
                           nrow = 1, ncol=2, plot=FALSE)

plot_tsne_2_ibd <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/function_t_sne_IBD.csv")
plot_tsne_2_ibd$disease <- rep("IBD", n=nrow(plot_tsne_2_ibd))
plot_tsne_2_crc <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/function_t_sne_CRC.csv")
plot_tsne_2_crc$disease <- rep("CRC", n=nrow(plot_tsne_2_crc))
plot_tsne_2 <- rbind(plot_tsne_2_ibd, plot_tsne_2_crc)
plot_tsne_2$function. <- factor(plot_tsne_2$function.)
plot_tsne_2$disease <- factor(plot_tsne_2$disease, levels = c("IBD", "CRC"))
middlle_plot_2 <- ggplot(plot_tsne_2, aes(x = t_SNE_1, y = t_SNE_2, color=group, shape=function.)) + 
    geom_point(fill="black", size=5) + 
    facet_wrap(~disease, scales = "free") +
    labs(x = "t_SNE 1", y = "t_SNE 2", fill="", shape="spoIIAA") +
    # scale_color_gradient(low = "white",high = "blue")+
    theme_bw(base_size = 14) + 
    theme(  
        plot.title = element_text(hjust = 0.5, size = 14),
        text = element_text(face = "bold"),
        plot.margin = margin(0, 0, 0, 0),  # 调整边距
        legend.position = "top"  # 设置图例位置为左下角
    ) +
    scale_shape_manual(values = c(1, 16))

## plot each gene
KO <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/pathway/picrust_filter/KO_predicted_description.tsv",
                 sep = "\t", header = 1, check.names = FALSE) %>% column_to_rownames("KO")
temp <- str_split(KO$description, ";", simplify = TRUE)
KO$description <- temp[,1]
pick_ko <- KO %>% filter(description %in% c("cysI", "tolB", "tolR", "msrB", "ugpB", "spoIIP", "spoIIAA", "spoVAC", "spoIIIAE", "spoVAE")) %>% 
    rownames_to_column("KO") %>% column_to_rownames("description")
plot_tsne_2_ibd <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/function_t_sne_IBD.csv")
plot_tsne_2_ibd$disease <- rep("IBD", n=nrow(plot_tsne_2_ibd))
plot_tsne_2_ibd$function. <- factor(plot_tsne_2_ibd$function.)
p8 <- ggplot(plot_tsne_2_ibd, aes(x = t_SNE_1, y = t_SNE_2, color=group, shape=function.)) + 
    geom_point(fill="black", size=5) + 
    facet_wrap(~disease, scales = "free") +
    labs(x = "t_SNE 1", y = "t_SNE 2", fill="", shape="cysI") +
    # scale_color_gradient(low = "white",high = "blue")+
    theme_bw(base_size = 14) + 
    theme(  
        plot.title = element_text(hjust = 0.5, size = 14),
        text = element_text(face = "bold"),
        plot.margin = margin(0, 0, 0, 0),  # 调整边距
        legend.position = "top"  # 设置图例位置为左下角
    ) +
    scale_shape_manual(values = c(1, 16))

merge_table <- pick_ko[, as.character(plot_tsne_2_ibd$fid)] %>% t() %>% as.data.frame()
merge_table[merge_table > 0] <- 1
plot_tsne_2_ibd <- cbind(plot_tsne_2_ibd, merge_table)

### 使用热图展示结果
# ratios <- c()
# for (i in rownames(pick_ko)){
#     ratios <- c(ratios, sum(plot_tsne_2_ibd[,i][plot_tsne_2_ibd$group == "Enrich_control"] == 0) / sum(plot_tsne_2_ibd$group == "Enrich_control"))
#     ratios <- c(ratios, sum(plot_tsne_2_ibd[,i][plot_tsne_2_ibd$group == "Enrich_control"] == 1) / sum(plot_tsne_2_ibd$group == "Enrich_control"))
#     ratios <- c(ratios, sum(plot_tsne_2_ibd[,i][plot_tsne_2_ibd$group == "Enrich_IBD"] == 0) / sum(plot_tsne_2_ibd$group == "Enrich_IBD"))
#     ratios <- c(ratios, sum(plot_tsne_2_ibd[,i][plot_tsne_2_ibd$group == "Enrich_IBD"] == 1) / sum(plot_tsne_2_ibd$group == "Enrich_IBD"))
# }
# heat_plot <- data.frame(enrich=rep(c("Enrich_control", "Enrich_IBD"), 10, each=2), 
#                         type=rep(c("0", "1"), 20), ratios=ratios, 
#                         group=rep(rownames(pick_ko), each=4))
# heat_plot$group <- factor(heat_plot$group, levels = c("cysI", "tolB", "tolR", "msrB", "ugpB", "spoIIP", "spoIIAA", "spoVAC", "spoIIIAE", "spoVAE"))
# heat_1 <- ggplot(heat_plot, aes(x = type, y = enrich, fill = ratios)) +  
#     geom_tile() +  
#     facet_wrap(~group, nrow = 2) +
#     scale_fill_gradient(low = "white", high = "red") +  # 设置颜色渐变  
#     labs(title = "", x = "", y = "", fill = "Ratio") +  
#     theme_bw() 
# ggsave("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/results/heatmap_gene_ibd.png", heat_1,
#        width = 20, height = 10, units = "cm")

### 使用柱状图展示
counts <- c()
for (i in rownames(pick_ko)){
    counts <- c(counts, sum(plot_tsne_2_ibd[,i][plot_tsne_2_ibd$group == "Enrich_control"] == 1))
    counts <- c(counts, sum(plot_tsne_2_ibd[,i][plot_tsne_2_ibd$group == "Enrich_IBD"] == 1))
}
bar_plot <- data.frame(enrich=rep(c("Enrich_control", "Enrich_IBD"), 10), 
                       counts=counts, group=rep(rownames(pick_ko), each=2))
bar_plot$group <- factor(bar_plot$group, levels = c("cysI", "tolB", "tolR", "msrB", "ugpB", "spoIIP", "spoIIAA", "spoVAC", "spoIIIAE", "spoVAE"))
bar_1 <- ggplot(bar_plot, aes(x = enrich, y = counts)) +  
    geom_bar(stat = "identity", fill = "steelblue") +
    facet_wrap(~group, nrow = 2, scales = "free") +
    scale_fill_gradient(low = "white", high = "red") +  # 设置颜色渐变  
    labs(title = "", x = "", y = "Counts", fill = "Ratio") +  
    theme_bw() +
    theme(axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1))
ggsave("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/results/bar_gene_ibd.png", bar_1,
       width = 20, height = 10, units = "cm")

plot_tsne_2_ibd_reshape <- melt(plot_tsne_2_ibd, id.vars = c("t_SNE_1", "t_SNE_2", "group", "emb_type", "fid", "function.", "disease"),
                                variable.name = "Gene_name", value.name = "Gene")
plot_tsne_2_ibd_reshape$Gene <- factor(plot_tsne_2_ibd_reshape$Gene)
plot_tsne_2_ibd_reshape$Gene_name <- factor(plot_tsne_2_ibd_reshape$Gene_name, levels = c("cysI", "tolB", "tolR", "msrB", "ugpB", "spoIIP", "spoIIAA", "spoVAC", "spoIIIAE", "spoVAE"))
sup_1 <- ggplot(plot_tsne_2_ibd_reshape, aes(x = t_SNE_1, y = t_SNE_2, color=group, shape=Gene)) + 
    geom_point(fill="black", size=5) + 
    facet_wrap(~Gene_name, scales = "free") +
    labs(x = "t_SNE 1", y = "t_SNE 2", fill="", title = "IBD") +
    # scale_color_gradient(low = "white",high = "blue")+
    theme_bw(base_size = 14) + 
    theme(  
        plot.title = element_text(hjust = 0.5, size = 14),
        text = element_text(face = "bold"),
        plot.margin = margin(0, 0, 0, 0),  # 调整边距
        legend.position = "top"  # 设置图例位置为左下角
    ) +
    scale_shape_manual(values = c(1, 16))


pick_ko <- KO %>% filter(description %in% c("dltB", "dltD", "dltC", "msrB", "dltA", "spoVAD", "spoIIAA", "spo0A", "spoIVA", "spoVAE")) %>% 
    rownames_to_column("KO") %>% column_to_rownames("description")
plot_tsne_2_crc <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/function_t_sne_CRC.csv")
plot_tsne_2_crc$disease <- rep("CRC", n=nrow(plot_tsne_2_crc))
plot_tsne_2_crc$function. <- factor(plot_tsne_2_crc$function.)
p9 <- ggplot(plot_tsne_2_crc, aes(x = t_SNE_1, y = t_SNE_2, color=group, shape=function.)) + 
    geom_point(fill="black", size=5) + 
    facet_wrap(~disease, scales = "free") +
    labs(x = "t_SNE 1", y = "t_SNE 2", fill="", shape="dltB") +
    # scale_color_gradient(low = "white",high = "blue")+
    theme_bw(base_size = 14) + 
    theme(  
        plot.title = element_text(hjust = 0.5, size = 14),
        text = element_text(face = "bold"),
        plot.margin = margin(0, 0, 0, 0),  # 调整边距
        legend.position = "top"  # 设置图例位置为左下角
    ) +
    scale_shape_manual(values = c(1, 16))


merge_table <- pick_ko[, as.character(plot_tsne_2_crc$fid)] %>% t() %>% as.data.frame()
merge_table[merge_table > 0] <- 1
plot_tsne_2_crc <- cbind(plot_tsne_2_crc, merge_table)

### 使用热图展示结果
# ratios <- c()
# for (i in rownames(pick_ko)){
#     ratios <- c(ratios, sum(plot_tsne_2_crc[,i][plot_tsne_2_crc$group == "Enrich_control"] == 0) / sum(plot_tsne_2_crc$group == "Enrich_control"))
#     ratios <- c(ratios, sum(plot_tsne_2_crc[,i][plot_tsne_2_crc$group == "Enrich_control"] == 1) / sum(plot_tsne_2_crc$group == "Enrich_control"))
#     ratios <- c(ratios, sum(plot_tsne_2_crc[,i][plot_tsne_2_crc$group == "Enrich_CRC"] == 0) / sum(plot_tsne_2_crc$group == "Enrich_CRC"))
#     ratios <- c(ratios, sum(plot_tsne_2_crc[,i][plot_tsne_2_crc$group == "Enrich_CRC"] == 1) / sum(plot_tsne_2_crc$group == "Enrich_CRC"))
# }
# heat_plot <- data.frame(enrich=rep(c("Enrich_control", "Enrich_CRC"), 10, each=2), 
#                         type=rep(c("0", "1"), 20), ratios=ratios, 
#                         group=rep(rownames(pick_ko), each=4))
# heat_plot$group <- factor(heat_plot$group, levels = c("dltB", "dltD", "dltC", "msrB", "dltA", "spoVAD", "spoIIAA", "spo0A", "spoIVA", "spoVAE"))
# heat_2 <- ggplot(heat_plot, aes(x = type, y = enrich, fill = ratios)) +  
#     geom_tile() +  
#     facet_wrap(~group, nrow = 2) +
#     scale_fill_gradient(low = "white", high = "red") +  # 设置颜色渐变  
#     labs(title = "", x = "", y = "", fill = "Ratio") +  
#     theme_bw() 
# ggsave("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/results/heatmap_gene_crc.png", heat_2,
#        width = 20, height = 10, units = "cm")

### 使用柱状图展示
counts <- c()
for (i in rownames(pick_ko)){
    counts <- c(counts, sum(plot_tsne_2_crc[,i][plot_tsne_2_crc$group == "Enrich_control"] == 1))
    counts <- c(counts, sum(plot_tsne_2_crc[,i][plot_tsne_2_crc$group == "Enrich_CRC"] == 1))
}
bar_plot <- data.frame(enrich=rep(c("Enrich_control", "Enrich_CRC"), 10), 
                       counts=counts, group=rep(rownames(pick_ko), each=2))
bar_plot$group <- factor(bar_plot$group, levels = c("dltB", "dltD", "dltC", "msrB", "dltA", "spoVAD", "spoIIAA", "spo0A", "spoIVA", "spoVAE"))
bar_2 <- ggplot(bar_plot, aes(x = enrich, y = counts)) +  
    geom_bar(stat = "identity", fill = "steelblue") +
    facet_wrap(~group, nrow = 2, scales = "free") +
    scale_fill_gradient(low = "white", high = "red") +  # 设置颜色渐变  
    labs(title = "", x = "", y = "Counts", fill = "Ratio") +  
    theme_bw() +
    theme(axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1))
ggsave("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/results/bar_gene_crc.png", bar_2,
       width = 20, height = 10, units = "cm")

plot_tsne_2_crc_reshape <- melt(plot_tsne_2_crc, id.vars = c("t_SNE_1", "t_SNE_2", "group", "emb_type", "fid", "function.", "disease"),
                                variable.name = "Gene_name", value.name = "Gene")
plot_tsne_2_crc_reshape$Gene <- factor(plot_tsne_2_crc_reshape$Gene)
plot_tsne_2_crc_reshape$Gene_name <- factor(plot_tsne_2_crc_reshape$Gene_name, levels = c("dltB", "dltD", "dltC", "msrB", "dltA", "spoVAD", "spoIIAA", "spo0A", "spoIVA", "spoVAE"))
sup_2 <- ggplot(plot_tsne_2_crc_reshape, aes(x = t_SNE_1, y = t_SNE_2, color=group, shape=Gene)) + 
    geom_point(fill="black", size=5) + 
    facet_wrap(~Gene_name, scales = "free") +
    labs(x = "t_SNE 1", y = "t_SNE 2", fill="", title = "CRC") +
    # scale_color_gradient(low = "white",high = "blue")+
    theme_bw(base_size = 14) + 
    theme(  
        plot.title = element_text(hjust = 0.5, size = 14),
        text = element_text(face = "bold"),
        plot.margin = margin(0, 0, 0, 0),  # 调整边距
        legend.position = "top"  # 设置图例位置为左下角
    ) +
    scale_shape_manual(values = c(1, 16))

bottom_plot_1 <- plot_grid(p8, p9,
                            nrow = 1, ncol=2, plot=FALSE)
# plot_tsne_2_crc %>% group_by(group) %>% summarise(pathogenic_mean = mean(pathogenic))
# feces_imp_disease <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/feces_imp_rf.csv")
# diseases <- unique(feces_imp_disease$group)
# types <- unique(feces_imp_disease$features)
# res_auc_plot <- list()
# n <- 1
# auc_res <- c()
# for (i in diseases){
#     for ( j in types){
#         # 计算 ROC 曲线  
#         temp <- feces_imp_disease %>% filter(group == i) %>% filter(features == j)
#         roc_result <- roc(temp$label, temp$probs)
#         # 绘制 ROC 曲线  
#         roc_data <- data.frame(  
#             Sensitivity = rev(roc_result$sensitivities),  # 真实阳性率  
#             Specificity = rev(roc_result$specificities),  # 真实阴性率  
#             Thresholds = roc_result$thresholds  # 阈值  
#         )  
#         roc_data$group <- rep(i, nrow(roc_data))
#         roc_data$features <- rep(j, nrow(roc_data))
#         res_auc_plot[[n]] = roc_data
#         auc_res <- c(auc_res, round(auc(roc(temp$label, temp$probs)), 3))
#         n <- n + 1
#     }
# }
# res_auc_plot_1 <- rbind(res_auc_plot[[1]], res_auc_plot[[2]], res_auc_plot[[3]], res_auc_plot[[4]])
# 
# ko_predict <- read.csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/feces_imp_rf_ko.csv")
# types <- unique(ko_predict$features)
# res_auc_plot <- list()
# n <- 1
# auc_res <- c()
# for ( j in types){
#     # 计算 ROC 曲线  
#     temp <- ko_predict %>% filter(features == j)
#     roc_result <- roc(temp$label, temp$probs)
#     # 绘制 ROC 曲线  
#     roc_data <- data.frame(  
#         Sensitivity = rev(roc_result$sensitivities),  # 真实阳性率  
#         Specificity = rev(roc_result$specificities),  # 真实阴性率  
#         Thresholds = roc_result$thresholds  # 阈值  
#     )  
#     roc_data$features <- rep(j, nrow(roc_data))
#     res_auc_plot[[n]] = roc_data
#     auc_res <- c(auc_res, round(auc(roc(temp$label, temp$probs)), 3))
#     n <- n + 1
# }
# res_auc_plot_2 <- rbind(res_auc_plot[[1]], res_auc_plot[[2]])
# res_auc_plot_2$group <- rep("msrB", nrow(res_auc_plot_2)) 
# res_auc_plot <- rbind(res_auc_plot_1, res_auc_plot_2)
# 
# res_auc_plot$group <- factor(res_auc_plot$group, levels = c("IBD", "CRC", "msrB"))
# bottom_plot_1 <- ggplot(res_auc_plot, aes(x = 1 - Specificity, y = Sensitivity, color=features)) +  
#     geom_line(size = 1) +           # 画曲线  
#     facet_wrap(~group) +
#     scale_color_manual(values=c("Imp."="#8ECFC9", "Unimp."="#FFBE7A")) +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + # 参考线  
#     labs(color="",  
#          x = "1 - Specificity",  
#          y = "Sensitivity") +  
#     # annotate("text", x = 0.75, y = 0.4, label = paste("AUC:", auc_res[1]), size = 4, color="#8ECFC9") +
#     # annotate("text", x = 0.75, y = 0.32, label = paste("AUC:", auc_res[2]), size = 4, color="#FFBE7A") +
#     theme_bw(base_size=14)+
#     theme(text = element_text(face = "bold"),
#           legend.position = "top",  # 设置图例位置为左下角
#           plot.title = element_text(hjust = 0.5, size = 14))
# 
# get_auc_f1 <- function(file, idx){
#     table <- read.table(file, sep = "\t")
#     table <- table[apply(table, 1, function(x) grepl("mcc", x)), ]
#     table <- gsub(",", "", table)
#     table <- str_split(table, " ", simplify = TRUE)
#     table <- as.data.frame(table[, idx])
# }
# 
# 
# Study <- c('PRJNA422193', 'PRJNA431126', 'qiita_1629', 'qiita_2538', 'HMP', 'PRJEB13679', 'PRJNA317429')
# attention_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/classification/IBD/attention.out"
# table <- get_auc_f1(attention_res, c(12, 18))
# colnames(table) <- c("AUC", "F1")
# table$Study <- rep(Study, 3, each=1)
# table$AUC <- as.numeric(table$AUC)
# table$Group <- rep(c("Original", "Imp. properties", "Unimp. properties"), 1, each= 7)
# table$Group <- factor(table$Group, levels = c("Original", "Imp. properties", "Unimp. properties"))
# p9 <- table %>% ggplot(aes(x = Group, y = AUC, group = Group)) +
#     geom_boxplot() + geom_line(aes(group = Study, color = Study)) +
#     geom_point(aes(color = Study), size = 3) +
#     theme_bw(base_size = 14) +
#     labs(x="", y="AUC", title = "IBD") +
#     theme(text = element_text(face = "bold"), # 设置所有文本为加粗
#           plot.title = element_text(hjust = 0.5), # 居中对齐标题
#           axis.title = element_text(face = "bold"), # 加粗坐标轴标题
#           axis.text.x = element_text(face = "bold", angle = 15, hjust = 1, vjust = 1)) +
#     scale_color_brewer(palette = "Paired")
# 
# 
# Study <- c('PRJEB36789', 'PRJNA824020', 'PRJDB11845', 'PRJNA318004', 'PRJEB6070', 'PRJNA430990', 'PRJNA290926')
# attention_res <- "/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/enrich_analysis/classification/CRC/attention.out"
# table <- get_auc_f1(attention_res, c(12, 18))
# colnames(table) <- c("AUC", "F1")
# table$Study <- rep(Study, 3, each=1)
# table$AUC <- as.numeric(table$AUC)
# table$Group <- rep(c("Original", "Imp. properties", "Unimp. properties"), 1, each= 7)
# table$Group <- factor(table$Group, levels = c("Original", "Imp. properties", "Unimp. properties"))
# p10 <- table %>% ggplot(aes(x = Group, y = AUC, group = Group)) +
#     geom_boxplot() + geom_line(aes(group = Study, color = Study)) +
#     geom_point(aes(color = Study), size = 3) +
#     theme_bw(base_size = 14) +
#     labs(x="", y="AUC", title = "CRC") +
#     theme(text = element_text(face = "bold"), # 设置所有文本为加粗
#           plot.title = element_text(hjust = 0.5), # 居中对齐标题
#           axis.title = element_text(face = "bold"), # 加粗坐标轴标题
#           axis.text.x = element_text(face = "bold", angle = 15, hjust = 1, vjust = 1)) +
#     scale_color_brewer(palette = "Paired")
# 
# bottom_plot_2 <- plot_grid(p9, p10,
#                            nrow = 1, ncol=2, plot=FALSE)
# p <- plot_grid(p1, p2, p3, plot_imp, middlle_plot_1, middlle_plot_2, bottom_plot_1, bottom_plot_2,
#                labels = c('a:IBD', 'b:CRC', 'c', 'd', 'e', 'f', "g", "h"),
#                align="hv",
#                nrow = 4, ncol=2, plot=FALSE, rel_heights = c(1, 1, 1))

p <- plot_grid(p1, p2, p3, plot_imp, middlle_plot_1, middlle_plot_2, p8, p9,
               labels = c('a:IBD', 'b:CRC', 'c', 'd', 'e', 'f', "g", "h"),
               align="hv",
               nrow = 4, ncol=2, plot=FALSE, rel_heights = c(1, 1, 1))

ggsave("/home/dongbiao/word_embedding_microbiome/result/embedding_explain_enrich_healthy_disease.png", p,
       width = 40, height = 40, units = "cm")
