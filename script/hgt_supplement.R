library(tidyverse)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggbreak)
library(scales)  
module_hgt <- read.csv("/home/dongbiao/word_embedding_microbiome/HGT/results/moduel_hgt.csv")

# 创建一个自定义的变换函数  
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

p <- module_hgt %>% ggplot(aes(x = Value, fill = Category)) +
    
    # 1. 绘制直方图（密度比例）
    geom_histogram(aes(y = ..density..),
                   position = "identity",
                   alpha = 0.5,
                   bins = 30) +
    
    # 2. 添加核密度曲线
    geom_density(alpha = 0.3, color = NA) +
    # 4. 添加p值标注
    labs(x = "log10(HGT rate % + 1)", 
         y = "Density",
         fill = "") +
    theme_bw(base_size = 14) +
    annotate("text",
             x = Inf, y = Inf,
             label = paste("p =", p_value),
             hjust = 1.1, vjust = 1.5,
             size = 3,
             color = "black",
             fontface = "italic") + 
    scale_fill_manual(values=c("between module"="orange","within module"="blue"),
                      labels=c("between module"="between guild","within module"="within guild")) +
    theme(text = element_text(face = "bold"),
        legend.position = c(0.8, 0.8),
        legend.background = element_blank(),
        legend.text = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        panel.grid.minor = element_blank()
    ) + coord_trans(y = squash_axis(1, 12, 20))
ggsave("/home/dongbiao/word_embedding_microbiome/result/hgt_guild.png", p,
       width = 8, height = 6, units = "cm")

# 绘制小提琴图
# p <- ggplot(module_hgt, aes(x = Category, y = Value, fill = Category)) +
#     geom_violin(trim = FALSE, alpha = 0.7) +  # 小提琴图
#     theme_bw(base_size=14) +  # 使用简洁主题
#     labs(title = "",
#          x = "",
#          y = "log10(HGT rate % + 1)", fill="") +
#     stat_compare_means(
#         method = "t.test",  # 使用t检验
#         comparisons = list(c("between module", "within module")),
#         label = "p.signif",  # 显示p值
#         label.y = c(2.5, 2.6, 2.7)  # 标注位置
#     ) +
#     theme(text = element_text(face = "bold"), # 设置所有文本为加粗
#           legend.position = "top",
#           plot.title = element_text(hjust = 0.5), # 居中对齐标题
#           axis.title = element_text(face = "bold"), # 加粗坐标轴标题
#           axis.text = element_text(face = "bold"))
