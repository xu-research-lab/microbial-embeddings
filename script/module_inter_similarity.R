library(dplyr)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(lemon)
library(cowplot)
library(reshape2)
library(tidyverse)
library(stringr)

library(ggdensity)
library(ggblanket)
library(ggsci)
library(reshape2)
library(gridExtra)
library(patchwork)
library(cowplot)
library(aplot)
library(ggplotify)

library(RColorBrewer)  

library(ggVennDiagram)

get_cohen_d <- function(data1, data2) {
    mean_diff <- mean(data1) - mean(data2)
    pooled_var <- (var(data1) + var(data2)) / 2
    pooled_sd <- sqrt(pooled_var)
    cohen_d <- mean_diff / pooled_sd
    
    return(cohen_d)
}

plot_res <- function(cohen_d, title){
    mean_value <- mean(cohen_d)
    ggplot(data.frame(x = cohen_d), aes(x = x)) +
        # 1. 绘制频率分布直方图
        geom_histogram(
            aes(y = ..density..),  # 使用密度而非计数
            bins = 30,            # 设置分箱数
            fill = "skyblue", 
            color = "white",
            alpha = 0.7
        ) +
        # 2. 添加正态分布拟合曲线
        stat_function(
            fun = dnorm, 
            args = list(mean = mean(cohen_d), sd = sd(cohen_d)),
            color = "red", 
            linewidth = 1.2
        ) +
        geom_vline(xintercept = mean_value, color = "red", linetype = "dashed", size = 1) +
        # 4. 添加图例和标签
        labs(
            title = title,
            x = "Cohen_d",
            y = "Density"
        ) +
        
        # 5. 美化主题
        theme_bw(base_size=14) +
        theme(text = element_text(face = "bold"), # 设置所有文本为加粗
              legend.position = "None",
              plot.title = element_text(hjust = 0.5), # 居中对齐标题
              axis.title = element_text(face = "bold"), # 加粗坐标轴标题
              axis.text = element_text(face = "bold"))
}

p_plot <- list()
n <- 1
title_names <- c("Pathway", "dbcan", "SCFA")
for (i in c("pathway", "dbcan", "SCFA", "vfdb", "resistance")){
    table <- read.csv(paste0("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/cluster/module_function/res_", i, ".csv"))

    module_id <- unique(table$module_id)
    cohen_d <- c()
    for (i in module_id){
        temp <- table %>% filter(module_id == i)
        data1 <- temp %>% filter(group == "Module") %>% select(sim_pathway)
        data2 <- temp %>% filter(group == "Random") %>% select(sim_pathway)
        cohen_d <- c(cohen_d, get_cohen_d(data1$sim_pathway, data2$sim_pathway))
    
    }
    
    p_plot[[n]] <- plot_res(cohen_d, title_names[n])
    n <- n + 1
}

p1 <- plot_grid(p_plot[[1]], p_plot[[2]], p_plot[[3]],
               align="hv", labels = c('a', 'c', 'd'),
               nrow = 1, ncol=3, plot=FALSE, rel_heights = c(1, 1))

p_plot <- list()
i <- "pathway"
table <- read.csv(paste0("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/cluster/module_function/res_", i, ".csv"))
module_id <- unique(table$module_id)
cohen_d <- c()
for (i in module_id){
    temp <- table %>% filter(module_id == i)
    data1 <- temp %>% filter(group == "Module") %>% select(sim_pathway)
    data2 <- temp %>% filter(group == "Random") %>% select(sim_pathway)
    cohen_d <- c(cohen_d, get_cohen_d(data1$sim_pathway, data2$sim_pathway))
    
}
module_cohen <- data.frame(module_id=module_id, cohen_d) %>% 
    arrange(desc(cohen_d))
max_id <- module_cohen[1,1]
temp <- table %>% filter(module_id == max_id)

# 计算均值和t检验结果
mean_values <- aggregate(sim_pathway ~ group, temp, mean)
t_test <- t.test(sim_pathway ~ group, data = temp)
p_value <- format.pval(t_test$p.value, digits = 3)

# 创建基础图形
p_plot[[1]] <- ggplot(temp, aes(x = sim_pathway, fill = group)) +
    
    # 1. 绘制直方图（密度比例）
    geom_histogram(aes(y = ..density..),
                   position = "identity",
                   alpha = 0.5,
                   bins = 30) +
    
    # 2. 添加核密度曲线
    geom_density(alpha = 0.3, color = NA) +
    
    # 3. 添加均值垂直线
    geom_vline(data = mean_values,
               aes(xintercept = sim_pathway, color = group),
               linetype = "dashed",
               linewidth = 0.8,
               show.legend = FALSE) +
    
    # 4. 添加p值标注
    annotate("text",
             x = Inf, y = Inf,
             label = paste("p =", p_value),
             hjust = 1.1, vjust = 1.5,
             size = 3,
             color = "black",
             fontface = "italic") +
    
    # 5. 设置颜色和标签
    scale_fill_manual(values = c("Module" = "blue", "Random" = "orange"), 
                      labels=c("Module"=paste0("Guild", "_", max_id),"Random"="Random")) +
    scale_color_manual(values = c("Module" = "blue", "Random" = "orange"),
                       labels=c("Module"=paste0("Guild", "_", max_id),"Random"="Random")) +
    labs(x = "Pathway similarity", 
         y = "Density",
         fill = "Group") +
    
    # 6. 调整主题
    theme_bw(base_size = 14) +
    theme(text = element_text(face = "bold"),
        legend.position = c(0.15, 0.85),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        panel.grid.minor = element_blank()
    )


i <- "dbcan"
table <- read.csv(paste0("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/cluster/module_function/res_", i, ".csv"))
module_id <- unique(table$module_id)
cohen_d <- c()
for (i in module_id){
    temp <- table %>% filter(module_id == i)
    data1 <- temp %>% filter(group == "Module") %>% select(sim_pathway)
    data2 <- temp %>% filter(group == "Random") %>% select(sim_pathway)
    cohen_d <- c(cohen_d, get_cohen_d(data1$sim_pathway, data2$sim_pathway))
    
}
module_cohen <- data.frame(module_id=module_id, cohen_d) %>% 
    arrange(desc(cohen_d))
max_id <- module_cohen[1,1]
temp <- table %>% filter(module_id == max_id)

# 计算均值和t检验结果
mean_values <- aggregate(sim_pathway ~ group, temp, mean)
t_test <- t.test(sim_pathway ~ group, data = temp)
p_value <- format.pval(t_test$p.value, digits = 3)

# 创建基础图形
p_plot[[2]] <- ggplot(temp, aes(x = sim_pathway, fill = group)) +
    
    # 1. 绘制直方图（密度比例）
    geom_histogram(aes(y = ..density..),
                   position = "identity",
                   alpha = 0.5,
                   bins = 30) +
    
    # 2. 添加核密度曲线
    geom_density(alpha = 0.3, color = NA) +
    
    # 3. 添加均值垂直线
    geom_vline(data = mean_values,
               aes(xintercept = sim_pathway, color = group),
               linetype = "dashed",
               linewidth = 0.8,
               show.legend = FALSE) +
    
    # 4. 添加p值标注
    annotate("text",
             x = Inf, y = Inf,
             label = paste("p =", p_value),
             hjust = 1.1, vjust = 1.5,
             size = 3,
             color = "black",
             fontface = "italic") +
    
    # 5. 设置颜色和标签
    scale_fill_manual(values = c("Module" = "blue", "Random" = "orange"), 
                      labels=c("Module"=paste0("Guild", "_", max_id),"Random"="Random")) +
    scale_color_manual(values = c("Module" = "blue", "Random" = "orange"),
                       labels=c("Module"=paste0("Guild", "_", max_id),"Random"="Random")) +
    labs(x = "dbcan similarity", 
         y = "Density",
         fill = "Group") +
    
    # 6. 调整主题
    theme_bw(base_size = 14) +
    theme(text = element_text(face = "bold"),
          legend.position = c(0.15, 0.85),
          legend.background = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 8),
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 8),
          panel.grid.minor = element_blank()
    )


i <- "SCFA"
table <- read.csv(paste0("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/cluster/module_function/res_", i, ".csv"))
module_id <- unique(table$module_id)
cohen_d <- c()
for (i in module_id){
    temp <- table %>% filter(module_id == i)
    data1 <- temp %>% filter(group == "Module") %>% select(sim_pathway)
    data2 <- temp %>% filter(group == "Random") %>% select(sim_pathway)
    cohen_d <- c(cohen_d, get_cohen_d(data1$sim_pathway, data2$sim_pathway))
    
}
module_cohen <- data.frame(module_id=module_id, cohen_d) %>% 
    arrange(desc(cohen_d))
max_id <- module_cohen[1,1]
temp <- table %>% filter(module_id == max_id)

# 计算均值和t检验结果
mean_values <- aggregate(sim_pathway ~ group, temp, mean)
t_test <- t.test(sim_pathway ~ group, data = temp)
p_value <- format.pval(t_test$p.value, digits = 3)

# 创建基础图形
p_plot[[3]] <- ggplot(temp, aes(x = sim_pathway, fill = group)) +
    
    # 1. 绘制直方图（密度比例）
    geom_histogram(aes(y = ..density..),
                   position = "identity",
                   alpha = 0.5,
                   bins = 30) +
    
    # 2. 添加核密度曲线
    # geom_density(alpha = 0.3, color = NA) +
    
    # 3. 添加均值垂直线
    # geom_vline(data = mean_values,
    #            aes(xintercept = sim_pathway, color = group),
    #            linetype = "dashed",
    #            linewidth = 0.8,
    #            show.legend = FALSE) +
    
    # 4. 添加p值标注
    annotate("text",
             x = Inf, y = Inf,
             label = paste("p =", p_value),
             hjust = 1.1, vjust = 1.5,
             size = 3,
             color = "black",
             fontface = "italic") +
    
    # 5. 设置颜色和标签
    scale_fill_manual(values = c("Module" = "blue", "Random" = "orange"), 
                      labels=c("Module"=paste0("Guild", "_", max_id),"Random"="Random")) +
    scale_color_manual(values = c("Module" = "blue", "Random" = "orange"),
                       labels=c("Module"=paste0("Guild", "_", max_id),"Random"="Random")) +
    labs(x = "SCFA similarity", 
         y = "Density",
         fill = "Group") +
    
    # 6. 调整主题
    theme_bw(base_size = 14) +
    theme(text = element_text(face = "bold"),
          legend.position = c(0.15, 0.85),
          legend.background = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 8),
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 8),
          panel.grid.minor = element_blank()
    )


p2 <- plot_grid(p_plot[[1]], p_plot[[2]], p_plot[[3]],
                align="hv", labels = c('b', 'd', 'e'), 
                nrow = 1, ncol=3, plot=FALSE, rel_heights = c(1, 1))


p <- plot_grid(p1, p2,
               align="hv",
               nrow = 2, ncol=1, plot=FALSE, rel_heights = c(1, 1, 1))

ggsave("/home/dongbiao/word_embedding_microbiome/result/module_function_similarity.png", p,
       width = 30, height = 15, units = "cm")
