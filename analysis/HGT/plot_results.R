library(tidyverse)
library(ggpubr)
library(readr)
library(biomformat)
library(cowplot)
library(pROC)  
library(scales)
library(ggrepel)
library(patchwork)
library(ggpubr) 
library(ggsci)


table_group <- read.csv("Data/tax_group.csv")
table_group$tax <- factor(table_group$tax, levels = c("Phylum", "Class", "Order", "Family", "Genus"), ordered = TRUE)


table_group <- table_group %>%
    group_by(tax) %>%
    mutate(
        r = cor(embed_cosine_mean, phy_dist_mean), # Pearson 
        slope = lm(phy_dist_mean ~ embed_cosine_mean)$coefficients[2], 
        intercept = lm(phy_dist_mean ~ embed_cosine_mean)$coefficients[1] 
    )
    
p1 <- ggplot(table_group, aes(x = embed_cosine_mean, y = phy_dist_mean)) +

    geom_errorbar(aes(ymin = phy_dist_mean - phy_dist_sd, ymax = phy_dist_mean + phy_dist_sd),
                  width = 0.02, color = "gray", linewidth = 0.3
    ) +
    geom_errorbarh(aes(xmin = embed_cosine_mean - embed_cosine_sd, xmax = embed_cosine_mean + embed_cosine_sd),
                   height = 0.02, color = "gray", linewidth = 0.3
    ) +
    geom_point(size = 2, alpha = 0.4, color = "blue") +
    geom_abline(aes(slope = slope, intercept = intercept), color = "red", linetype = "dashed", linewidth = 0.8) +
    geom_text(aes(x = 0.2, y = 1.5, label = paste0("r = ", round(r, 2))), size = 4, color = "black") +

    facet_wrap(~tax, nrow = 2) +

    labs(
        x = "SNE cosine",
        y = "Phylo distance",
        title = ""
    ) +

    theme_bw() +
    theme(
        legend.position = "top",
        plot.title = element_text(hjust = 0.5),
    )


### SNE Phy
hgt_plot_res <- read.csv("/home/dongbiao/word_embedding_microbiome/HGT/results/hgt_plot_res.csv")
hgt_plot_res$distance[hgt_plot_res$group == "SNE"] <- 1 - hgt_plot_res$distance[hgt_plot_res$group == "SNE"]

df <- hgt_plot_res %>% filter(group == "Phylo")
summary(lm(hgt_rate ~ distance, df))

phylogeny_data <- hgt_plot_res %>% filter(group != "SNE")
x_primary_range <- range(phylogeny_data$distance)

coembed_data <- hgt_plot_res %>% filter(group == "SNE") %>% 
    mutate(scaled_x = rescale(distance, to = x_primary_range))

coembed_data$scaled_x <- rescale(coembed_data$distance, 
                                 to = range(phylogeny_data$distance))
primary_range <- range(phylogeny_data$distance)
secondary_range <- range(1 - coembed_data$distance) 
phylogeny_data$group <- factor(phylogeny_data$group, levels = c("Phylo", "Phylo SNEsim > 0.6", "Phylo SNEsim < 0.6"))
Phylo_paired <- phylogeny_data %>%
    dplyr::filter(group != "Phylo")

group1_data <- Phylo_paired$hgt_rate[Phylo_paired$group == "Phylo SNEsim > 0.6"]
group2_data <- Phylo_paired$hgt_rate[Phylo_paired$group == "Phylo SNEsim < 0.6"]

paired_test_res <- t.test(group2_data, group1_data, paired = TRUE)


p2 <- ggplot() +
    # Phylogeny数据层（三个子组）
    geom_point(
        data = phylogeny_data,
        aes(x = distance, y = hgt_rate, color = group, shape = group),
        size = 2
    ) +
    geom_smooth(
        data = phylogeny_data,
        aes(x = distance, y = hgt_rate, color = group, linetype = group),
        method = "loess", se = TRUE, size = 0.6, level = 0.85
    ) +
    # 修正Phylo组的R²标签
    ggpubr::stat_cor(
        data = phylogeny_data[phylogeny_data$group == "Phylo", ],
        aes(x = distance, y = hgt_rate, label = ..rr.label..),
        label.x.npc = 0.2,  
        label.y.npc = 0.55,    # 标签y轴位置
        color = "#E41A1C",      # 强制指定颜色与Phylo组一致
        size = 3,
        show.legend = FALSE
    ) +
    # Co-embedding数据层
    geom_point(
        data = coembed_data,
        aes(x = scaled_x, y = hgt_rate, color = "SNE"),
        shape = 16, size = 2  
    ) +
    geom_smooth(
        data = coembed_data,
        aes(x = scaled_x, y = hgt_rate, color = "SNE"),
        method = "loess", se = TRUE, size = 0.6
    ) +
    # 修正SNE组的R²标签
    ggpubr::stat_cor(
        data = coembed_data,
        aes(x = scaled_x, y = hgt_rate, label = ..rr.label..),
        label.x.npc = 0.2,
        label.y.npc = 0.3,
        color = "#377EB8",       # 强制指定颜色与SNE组一致
        size = 3,
        show.legend = FALSE
    ) +
    
    # --- 双轴系统配置 ---
    scale_x_continuous(
        name = "Phylo distance",
        limits = range(phylogeny_data$distance),
        sec.axis = sec_axis(
            ~ rescale(., to = range(1 - coembed_data$distance)),
            name = "1 - SNESim",
            breaks = pretty_breaks(n = 5)
        )
    ) +
    
    # --- 视觉映射系统 ---
    scale_color_manual(
        name = "Group",
        values = c(
            "Phylo" = "#E41A1C",
            "Phylo SNEsim > 0.6" = "#E41A1C",
            "Phylo SNEsim < 0.6" = "#E41A1C",
            "SNE" = "#377EB8"
        ),
    ) +
    scale_shape_manual(
        name = "Group",
        values = c(
            "Phylo" = 16,      # 圆形
            "Phylo SNEsim > 0.6" = 17, # 三角形
            "Phylo SNEsim < 0.6" = 15,  # 正方形
            "SNE" = 16    # 菱形
        ),
        guide = "none"  # 隐藏独立线型图例
    ) +
    scale_linetype_manual(
        name = "Group",
        values = c(
            "Phylo" = "solid",
            "Phylo SNEsim > 0.6" = "dashed",
            "Phylo SNEsim < 0.6" = "dotted",
            "SNE" = "solid"
        ),
        guide = "none"  # 隐藏独立线型图例
    ) +
    
    # --- 图例合并系统 ---
    guides(
        color = guide_legend(
            override.aes = list(
                linetype = c("solid", "dashed", "dotted", "solid"),
                shape = c(16, 17, 15, 16)
            )
        )
    ) +
    geom_text(
        data = data.frame(
            x = 0.8, 
            y = 20, 
            label = paste0("p-value: ", 
                           format(paired_test_res$p.value, format = "e", digits = 3, scientific=TRUE))
        ),
        aes(x = x, y = y, label = label),
        size = 3.5, 
        color = "black"
    )

# --- 主题优化 ---
p2 <- p2 + theme_bw(base_size = 12) +
    theme(
        legend.title = element_blank(),
        axis.title.x.bottom = element_text(color = "#E41A1C", face = "bold"),
        axis.text.x.bottom = element_text(color = "#E41A1C"),
        axis.title.x.top = element_text(color = "#377EB8", face = "bold"),
        axis.text.x.top = element_text(color = "#377EB8")
    ) +
    labs(y = "HGT Per 100 Comparisons")

### HGT predict
hgt_predict_res <- read.csv("/home/dongbiao/word_embedding_microbiome/HGT/results/hgt_predict_res_all.csv")
auc_sne <- c()
auc_phyloe <- c()
auc_sne_phyloe <- c()
for (i in c(0:4)){
    temp <- hgt_predict_res %>% filter(test == paste0("test_", i)) %>% filter(group == "SNE")
    roc_obj <- roc(temp$labels, temp$proba)
    auc_sne <- c(auc_sne, auc(roc_obj))
    temp <- hgt_predict_res %>% filter(test == paste0("test_", i)) %>% filter(group == "PhyloE")
    roc_obj <- roc(temp$labels, temp$proba)
    auc_phyloe <- c(auc_phyloe, auc(roc_obj))
    temp <- hgt_predict_res %>% filter(test == paste0("test_", i)) %>% filter(group == "SNE+PhyloE")
    roc_obj <- roc(temp$labels, temp$proba)
    auc_sne_phyloe <- c(auc_sne_phyloe, auc(roc_obj))
}

auc_res <- data.frame(fold = rep(c("test_1", "test_2", "test_3", "test_4", "test_5"), 3),
                      auc = c(auc_sne, auc_phyloe, auc_sne_phyloe),
                      group = c(rep("SNE", 5), rep("PhyloE", 5), rep("SNE_PhyloE", 5)))
auc_res$group <- factor(auc_res$group, levels = c("SNE", "PhyloE", "SNE_PhyloE"))
line_plot <- auc_res %>% group_by(group) %>% 
    summarise(auc_mean = mean(auc), auc_sd = sd(auc))
auc_res <- merge(auc_res, line_plot, by=c("group"))
p1 <- auc_res %>% ggplot(aes(x=group, y=auc)) +
    geom_point() +
    geom_errorbar(aes(ymin = auc_mean - auc_sd, ymax = auc_mean + auc_sd), width = 0.1) + 
    geom_line(aes(group=fold, y = auc_mean), alpha = 0.5, linewidth = 1, linetype = "solid") +
    theme_bw(base_size = 14) + 
    labs(x="", y="AUC") +
    theme(text = element_text(face = "bold"))
imp_features <- read.csv("/home/dongbiao/word_embedding_microbiome/HGT/results/hgt_importance_dim.csv")
temp <- imp_features$importance[1:400] + imp_features$importance[401:800] + imp_features$importance[801:1200] + imp_features$importance[1201:1600] + imp_features$importance[1601:2000]
temp <- data.frame(importance = temp)
temp$group <- c(rep("SNE", 100), rep("PhyloE", 100), rep("SNE", 100), rep("PhyloE", 100))


mean_values <- temp %>%
    group_by(group) %>%
    summarise(mean_importance = mean(importance))

# 2. 统计检验（t检验或Wilcoxon检验）
p_value = t.test(importance ~ group, data = temp)$p.value
p_text <- paste("p =", format.pval(p_value, digits = 2))

# 3. 绘图
p2 <- ggplot(temp, aes(x = importance, fill = group)) +
    # 密度曲线
    geom_density(aes(y=after_stat(density)), alpha = 0.5) +
    # 添加均值垂直线
    geom_vline(
        data = mean_values,
        aes(xintercept = mean_importance, color = group),
        linetype = "dashed", linewidth = 0.8, show.legend = FALSE
    ) +
    # 添加p值标注
    annotate(
        "text", x = Inf, y = Inf, label = paste0(p_text),
        hjust = 1.1, vjust = 1.5, size = 3, fontface = "bold"
    ) +
    # 颜色和主题设置
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Set2") +
    labs(x = "Importance", y = "Density", fill = "Group") +
    theme_bw(base_size = 10) +
    theme(
        text = element_text(face = "bold"),
        legend.position = "top",
        legend.title = element_blank()
    )

combined_academic <- p1 +
    inset_element(p2,
                  left = 0.55, bottom = 0.01,
                  right = 0.95, top = 0.6) +
    theme(plot.tag = element_text(size = 16, face = "bold"),
          plot.tag.position = c(0.02, 0.98))



### Phyloglm KO with SNE
res <- read.csv("Data/KO_phyloglm.csv")
res <- res[1:10, ]
res$descp <- c("K01547 kdpB; potassium-transporting ATPase \nATP-binding subunit",
               "K01548 kdpC; potassium-transporting \nATPase KdpC subunit",
               "K01546 kdpA; potassium-transporting \nATPase potassium-binding subunit",
               "K07646 kdpD; two-component system, OmpR family, \nsensor histidine kinase KdpD",
               "K03320 amt, AMT, MEP; ammonium transporter, \nAmt family",
               "K05846 opuBD; osmoprotectant transport \nsystem permease protein",
               "K14698 irtA, ybtP; ATP-binding cassette, \nsubfamily B, bacterial IrtA/YbtP",
               "K05845 opuC; osmoprotectant transport \nsystem substrate-binding protein",
               "K14699 irtB, ybtQ; ATP-binding cassette, \nsubfamily B, bacterial IrtB/YbtQ",
               "K06400 spoIVCA; site-specific DNA recombinase")
res$descp <- factor(res$descp, levels = res$descp)
p2 <- res %>% ggplot(aes(x=LR_stat, y=descp)) +
    geom_bar(stat = "identity", position = "dodge", fill = "skyblue") +
    theme(axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          legend.position = "right") +
    theme_bw() + 
    labs(x = "LR_stat", y = "", fill = "", title="") +
    coord_cartesian(xlim = c(200, 400)) +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold"),legend.position = "None")

