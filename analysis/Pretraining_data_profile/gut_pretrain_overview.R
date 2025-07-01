# load libraries
library(ggplot2)
library(purrr)
library(dplyr)
library(maps)
library(scales)
library(treemapify)
library(ggsci)
library(scales)
library(jsonlite)  
library(tidyverse)
library(patchwork)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(Polychrome)

# load data
metadata <- read.csv("Data/meta_data_include_shannon.csv", header = TRUE)

studies <- metadata %>%
    group_by(project) %>%
    summarise(n = n()) %>%
    arrange(desc(n))
# studies <- studies %>% filter(n > 9)
# 定义分组区间
breaks <- c(0, 50, 100, 150, 200, 300, 500, 1000, 2500, 5000, 7500, 10000, 25000)

# 定义分组标签（确保与 breaks 匹配）
labels <- c(
    "0-50", "51-100", "101-150", "151-200", "201-300",
    "301-500", "501-1000", "1001-2500", "2501-5000", "5001-7500", "7500-10000", ">20000"
)
# 使用 cut() 分组（推荐使用 right = TRUE 和 include.lowest = TRUE）
studies$class <- cut(
    studies$n,
    breaks = breaks,
    include.lowest = TRUE,  # 包含最小值（0）
    right = TRUE,           # 区间左开右闭 (a, b]
    labels = labels
)

p1 <- ggplot(studies, aes(area = n, fill = class, label = project)) +
    geom_treemap() +
    scale_fill_brewer(palette = "Set3")  +
    # geom_treemap_text(colour = "white", place = "centre") +
    labs(
        fill = "No. of samples",
        title = ""
    ) + 
    theme_bw() +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold")) 

### samples reads
df <- read.csv("Data/read_depth_samples.csv", header = FALSE)
colnames(df) <- "reads"

median_reads <- median(df$reads)
min_reads <- min(df$reads)
max_reads <- max(df$reads)

p2 <- ggplot(df, aes(x = reads)) +
    geom_histogram(
        breaks = 10^seq(log10(min_reads), log10(max_reads), length.out = 40),
        fill = "black", color = NA, alpha = 0.9) +
    geom_vline(
        xintercept = median_reads, color = "gray70",
        linetype = "solid", linewidth = 1.2
    ) +
    annotate(
        "text", x = median_reads,y = Inf, 
        label = paste("Median:", comma(median_reads)),
        vjust = 1.5, hjust = -0.1, color = "black", size = 4) +

    scale_x_log10(
        labels = trans_format("log10", math_format(10^.x)),
        breaks = trans_breaks("log10", function(x) 10^x),
        expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
        expand = expansion(mult = c(0, 0.1)), # 顶部留 10% 空间
        labels = comma
    ) +
    
    labs(
        x = "No. of reads",
        y = "No. of samples",
        title = ""
    ) +
    theme_bw(base_size = 14) +
    theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white"),
        
        panel.grid.major = element_line(color = "gray70",
            linewidth = 0.3, linetype = "dashed"
        ),
        panel.grid.minor = element_line(
            color = "gray90", linewidth = 0.2, linetype = "dotted"
        ),
        axis.line = element_line(color = "black", linewidth = 0.6),
        axis.ticks = element_line(color = "black", linewidth = 0.6),
        
        # 文本设置
        text = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black"))

### prevalence distribution
prevalence <- read.csv("Data/prevalence.csv")
min_prevalence <- min(prevalence$prevalence)
max_prevalence <- max(prevalence$prevalence)
median_prevalence <- median(prevalence$prevalence)

# 假设 prevalence 是一个数值向量  
p3 <- ggplot(prevalence, aes(x = prevalence)) +
    # 添加直方图层
    geom_histogram(bins = 30,   
                   fill = "black") +    
    
    # *** 应用对数刻度 (base 10) 到 x 轴 ***
    scale_x_log10(
        breaks = trans_breaks("log10", function(x) 10^x), # 主刻度设为10的幂
        labels = trans_format("log10", math_format(10^.x)), # 标签格式为 10^x
        minor_breaks = waiver() # 通常让 ggplot 自动计算合适的次刻度，它默认会产生 2*10^n, 5*10^n 等
        # 如果需要手动指定更密集的次刻度(如 2,3,4..9 * 10^n)，会复杂一些，但通常 waiver() 效果不错
    ) +
    
    # 添加刻度标记 (在底部 "b" 添加)
    annotation_logticks(sides = "b") +
    
    # 添加标签和标题
    labs(
        title = "",
        x = "log10(Prevalence)",
        y = "Frequency"
    ) +
    geom_vline(
        xintercept = median_prevalence, color = "gray70",
        linetype = "solid", linewidth = 1.2
    ) +
    annotate(
        "text", x = median_prevalence, y = Inf, 
        label = paste("Median:", round(median_prevalence, 4)),
        vjust = 1.5, hjust = -0.1, color = "black", size = 4) +
    theme_bw(base_size = 14) +
    theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white"),
        
        panel.grid.major = element_line(color = "gray70",
                                        linewidth = 0.3, linetype = "dashed"
        ),
        panel.grid.minor = element_line(
            color = "gray90", linewidth = 0.2, linetype = "dotted"
        ),
        axis.line = element_line(color = "black", linewidth = 0.6),
        axis.ticks = element_line(color = "black", linewidth = 0.6),
        
        # 文本设置
        text = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black"))


### tax sample nums
final_df <- read.csv("Data/tax_num_samples.csv", row.names = 1)
taxonomy_mapping <- fromJSON("Data/taxonomy_mapping.json") 

# 全局颜色配置
PHYLUM_COLORS <- list()
TOP_PHYLA <- c()
plot_top_taxa <- function(level, top_n = 10) {
    # 提取数据并清洗索引
    data <- final_df %>%
        select(!!paste0(level, "_Samples")) %>%
        arrange(desc(!!sym(paste0(level, "_Samples")))) %>%
        head(top_n)
    
    # 确保行名无多余空格
    row_names <- trimws(rownames(data))
    data_values <- data[[paste0(level, "_Samples")]]
    
    # 创建数据框用于ggplot
    plot_data <- data.frame(
        taxon = row_names,
        samples = data_values,
        stringsAsFactors = FALSE
    )
    
    # 颜色分配逻辑
    if (level == "Phylum") {
        # 为顶级门分配颜色
        TOP_PHYLA <<- plot_data$taxon[1:5]
        phylum_palette <- brewer.pal(n = 5, name = "Paired")
        PHYLUM_COLORS <<- setNames(phylum_palette, TOP_PHYLA)
        
        # 为所有门分配颜色
        plot_data$color <- sapply(plot_data$taxon, function(p) {
            if (p %in% names(PHYLUM_COLORS)) {
                return(PHYLUM_COLORS[[p]])
            } else {
                return("#CCCCCC")
            }
        })
    } else {
        # 为其他分类级别分配颜色
        plot_data$color <- sapply(plot_data$taxon, function(taxon) {
            clean_taxon <- tolower(trimws(taxon))
            phylum <- if (taxon %in% names(taxonomy_mapping[[level]])) {
                taxonomy_mapping[[level]][[taxon]]
            } else if (clean_taxon %in% names(taxonomy_mapping[[level]])) {
                taxonomy_mapping[[level]][[clean_taxon]]
            } else {
                "(Unassigned)"
            }
            
            if (phylum %in% names(PHYLUM_COLORS)) {
                return(PHYLUM_COLORS[[phylum]])
            } else {
                return("#CCCCCC")
            }
        })
    }
    
    # 将颜色作为因子处理，保持顺序
    plot_data$taxon <- factor(plot_data$taxon, levels = rev(plot_data$taxon))
    
    # 创建ggplot对象
    p <- ggplot(plot_data, aes(x = samples, y = taxon)) +
        geom_col(aes(fill = taxon)) +
        scale_fill_manual(values = setNames(plot_data$color, plot_data$taxon)) +
        geom_text(aes(label = taxon, x = 0.01), hjust = 0, 
                  color = "black", fontface = "bold", size = 3.5) +
        labs(
            x = "Number of Samples",
            y = level,
            title = ""
        ) +
        theme_bw(base_size=14) +
        theme(
            text = element_text(face = "bold"),
            panel.background = element_rect(fill = "white", color = "black", size = 0.8),
            panel.grid.major = element_line(linetype = "dashed", color = "gray70", size = 0.3),
            panel.grid.minor = element_line(linetype = "dotted", color = "gray40", size = 0.2),
            axis.title.x = element_text(face = "bold"),
            axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 15)),
            axis.text.y = element_blank(),
            legend.position = "none",
            plot.background = element_rect(fill = "white", color = NA)
        )
    
    # 格式化x轴
    p <- p + scale_x_continuous(
        labels = function(x) ifelse(x >= 1000, paste0(x %/% 1000, "k"), as.character(x)),
        limits = c(0, max(plot_data$samples) * 1.05)
    )
    
    return(p)
}
# 执行可视化
p4 <- plot_top_taxa("Phylum")
p5 <- plot_top_taxa("Family")
p6 <- plot_top_taxa("Order")

### top phyla abundance
BINS <- 30  # Number of bins, controls curve smoothness
LINE_WIDTH <- 1  # In ggplot2, line size is typically smaller than in matplotlib
phyla_abundance_top <- read.csv("Data/phyla_abundance_top.csv", row.names = 1)
phyla_abundance_top <- na.omit(phyla_abundance_top)
all_values <- unlist(phyla_abundance_top)
all_values <- as.vector(all_values[all_values > 0])
bin_edges <- seq(min(all_values), max(all_values), length.out = BINS + 1)

# Initialize empty dataframe for histogram data
hist_data <- data.frame()
top_phyla <- colnames(phyla_abundance_top)
# For each phylum, compute histogram and add to data frame
for (phylum in top_phyla) {
    values <- phyla_abundance_top[[phylum]][phyla_abundance_top[[phylum]] > 0]
    hist_result <- hist(values, breaks = bin_edges, plot = FALSE)
    
    # Add to dataframe
    phylum_hist <- data.frame(
        bin_center = hist_result$mids,
        count = hist_result$counts,
        phylum = phylum,
        label = paste0(phylum, " (n=", length(values), ")")
    )
    
    hist_data <- rbind(hist_data, phylum_hist)
}

# Create the plot
p7 <- ggplot(hist_data, aes(x = bin_center, y = count, color = label)) +
    geom_line(size = LINE_WIDTH) +
    scale_y_log10(
        labels = function(y) ifelse(y >= 1000, paste0(y/1000, "k"), as.character(y)),
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        minor_breaks = NULL  # Remove minor breaks
    ) +
    scale_x_continuous(
        labels = function(x) paste0(x, "%"),
        name = "Relative Abundance"
    ) +
    scale_color_manual(values = setNames(PHYLUM_COLORS[top_phyla], paste0(top_phyla, " (n=", sapply(top_phyla, function(p) sum(phyla_abundance_top[[p]] > 0)), ")"))) +
    labs(
        y = "Number of Samples",
        color = "", title = ""
    ) +
    theme_bw(base_size=14) +
    theme(text = element_text(face = "bold"), # 设置所有文本为加粗
          legend.position = "None",
          plot.title = element_text(hjust = 0.5), # 居中对齐标题
          axis.title = element_text(face = "bold"), # 加粗坐标轴标题
          axis.text = element_text(face = "bold")) 

### shannon diversity
meta_data <- read.csv("Data/meta_data_include_shannon.csv")
median_shannon <- median(meta_data$shannon)
min_shannon <- min(meta_data$shannon)
max_shannon <- max(meta_data$shannon)

p8 <- ggplot(meta_data, aes(x = shannon)) +
    geom_histogram(
        fill = "black", color = NA, alpha = 0.9) +
    geom_vline(
        xintercept = median_shannon, color = "gray70",
        linetype = "solid", linewidth = 1.2
    ) +
    annotate(
        "text", x = median_shannon, y = Inf, 
        label = paste("Median:", comma(median_shannon)),
        vjust = 1.5, hjust = -0.1, color = "black", size = 4) +
    
    scale_y_continuous(
        expand = expansion(mult = c(0, 0.1)), # 顶部留 10% 空间
        labels = comma
    ) +
    
    labs(
        x = "Shannon Diversity Index",
        y = "No. of samples",
        title = ""
    ) +
    theme_bw(base_size = 14) +
    theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white"),
        
        panel.grid.major = element_line(color = "gray70",
                                        linewidth = 0.3, linetype = "dashed"
        ),
        panel.grid.minor = element_line(
            color = "gray90", linewidth = 0.2, linetype = "dotted"
        ),
        axis.line = element_line(color = "black", linewidth = 0.6),
        axis.ticks = element_line(color = "black", linewidth = 0.6),
        
        # 文本设置
        text = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black"))

# load data
meta_data <- read.csv("/home/dongbiao/word_embedding_microbiome/all_data/gut/cell/meta_data_include_shannon.csv")
# 获取地图数据并标准化名称
world_map <- map_data("world")
# 创建区域分类标准
world_map <- world_map %>%
    mutate(Region = case_when(
        # Europe and N. America
        region %in% c(
            "Albania", "Andorra", "Austria", "Belarus", "Belgium", "Bosnia and Herzegovina",
            "Bulgaria", "Croatia", "Cyprus", "Czech Republic", "Denmark", "Estonia",
            "Faroe Islands", "Finland", "France", "Germany", "Greece", "Hungary",
            "Iceland", "Ireland", "Italy", "Latvia", "Liechtenstein", "Lithuania",
            "Luxembourg", "Malta", "Moldova", "Monaco", "Montenegro", "Netherlands",
            "North Macedonia", "Norway", "Poland", "Portugal", "Romania", "Russia",
            "San Marino", "Serbia", "Slovakia", "Slovenia", "Spain", "Sweden",
            "Switzerland", "Ukraine", "UK", "Vatican", "Kosovo", "Greenland", "USA",
            "Canada", "Guernsey", "Jersey", "Isle of Man"
        ) ~ "Europe and Northern America",
        
        # Eastern and SE Asia
        region %in% c(
            "China", "Japan", "South Korea", "North Korea", "Taiwan", "Mongolia",
            "Vietnam", "Thailand", "Malaysia", "Singapore", "Philippines", "Indonesia",
            "Cambodia", "Laos", "Myanmar", "Brunei", "Timor-Leste", "Hong Kong", "Macau"
        ) ~ "Eastern and South-Eastern Asia",
        
        # Sub-Saharan Africa
        region %in% c(
            "Angola", "Benin", "Botswana", "Burkina Faso", "Burundi", "Cameroon",
            "Cape Verde", "Central African Republic", "Chad", "Comoros",
            "Democratic Republic of the Congo", "Republic of Congo", "Djibouti",
            "Equatorial Guinea", "Eritrea", "Ethiopia", "Gabon", "Gambia", "Ghana",
            "Guinea", "Guinea-Bissau", "Ivory Coast", "Kenya", "Lesotho", "Liberia",
            "Madagascar", "Malawi", "Mali", "Mauritania", "Mauritius", "Mozambique",
            "Namibia", "Niger", "Nigeria", "Rwanda", "Sao Tome and Principe", "Senegal",
            "Seychelles", "Sierra Leone", "Somalia", "South Africa", "South Sudan",
            "Sudan", "Swaziland", "Tanzania", "Togo", "Uganda", "Zambia", "Zimbabwe"
        ) ~ "Sub-Saharan Africa",
        
        # Central and S. Asia
        region %in% c(
            "Afghanistan", "Bangladesh", "Bhutan", "India", "Kazakhstan", "Kyrgyzstan",
            "Maldives", "Nepal", "Pakistan", "Sri Lanka", "Tajikistan", "Turkmenistan",
            "Uzbekistan"
        ) ~ "Central and Southern Asia",
        
        # Australia/New Zealand
        region %in% c(
            "Australia", "New Zealand", "Papua New Guinea", "Fiji", "Solomon Islands",
            "Vanuatu", "Samoa", "Tonga", "Kiribati", "Marshall Islands", "Micronesia",
            "Palau", "Nauru", "Tuvalu", "Cook Islands", "Niue"
        ) ~ "Australia/New Zealand",
        
        # N. Africa and W. Asia
        region %in% c(
            "Algeria", "Bahrain", "Egypt", "Iran", "Iraq", "Israel", "Jordan",
            "Kuwait", "Lebanon", "Libya", "Morocco", "Oman", "Palestine", "Qatar",
            "Saudi Arabia", "Syria", "Tunisia", "Turkey", "United Arab Emirates", "Yemen",
            "Western Sahara", "Armenia", "Azerbaijan", "Georgia"
        ) ~ "Northern Africa and Western Asia",
        
        # Latin America/Caribbean
        region %in% c(
            "Argentina", "Belize", "Bolivia", "Brazil", "Chile", "Colombia",
            "Costa Rica", "Cuba", "Dominica", "Dominican Republic", "Ecuador",
            "El Salvador", "French Guiana", "Grenada", "Guatemala", "Guyana", "Haiti",
            "Honduras", "Jamaica", "Mexico", "Nicaragua", "Panama", "Paraguay", "Peru",
            "Puerto Rico", "Saint Lucia", "Saint Vincent", "Suriname", "Trinidad",
            "Tobago", "Uruguay", "Venezuela", "Barbados", "Bahamas", "Antigua", "Barbuda",
            "Cayman Islands", "Curacao", "Saint Kitts", "Nevis", "Saint Martin",
            "Turks and Caicos Islands", "Virgin Islands", "Aruba", "Anguilla",
            "Bonaire", "Sint Eustatius", "Saba", "Montserrat", "Bermuda"
        ) ~ "Latin America and the Caribbean",
        
        # 未分类的特殊区域
        TRUE ~ "Undefined"
    ))


world_map <- world_map[world_map$Region != "Undefined",]

map_colors <- c(
    "Europe and Northern America" = brewer.pal(7, "Set3")[1],
    "Eastern and South-Eastern Asia" = brewer.pal(7, "Set3")[2],
    "Sub-Saharan Africa" = brewer.pal(7, "Set3")[3],
    "Central and Southern Asia" = brewer.pal(7, "Set3")[4],
    "Australia/New Zealand" = brewer.pal(7, "Set3")[5],
    "Northern Africa and Western Asia" = brewer.pal(7, "Set3")[6],
    "Latin America and the Caribbean" = brewer.pal(7, "Set3")[7]
)


temp <- ggplot() + geom_map(
    data = world_map,
    map = world_map,
    aes(x = long, y = lat, map_id = region,fill = Region),
    color = NA,
    show.legend = F
) +
    scale_fill_manual(values = map_colors) +
    theme_minimal() +
    ylab(NULL) +
    xlab(NULL) +
    theme(
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank()
    ) 
temp

region = c("Europe and Northern America", "Eastern and South-Eastern Asia", "Central and Southern Asia",
           "Sub-Saharan Africa", "Australia/New Zealand", "Northern Africa and Western Asia", "Latin America and the Caribbean")
# Region 样本数统计
region_counts <- meta_data %>%
    group_by(region) %>%
    summarise(count = n()) %>%
    arrange(desc(count))

region_counts <- region_counts %>% na.omit()
region_counts <- region_counts[region_counts$region %in% region,]

region_name <- region_counts$region

#绘制柱状图
region_counts_plot <- ggplot(region_counts, aes(x = count, y = reorder(region, count), fill = region)) +
    geom_bar(stat = "identity") +
    geom_text(
        aes(label = region_name),
        x = 0,
        hjust = 0,
        color = "black",
        size = 12 / .pt  # 保持与主题文字大小一致（12pt）
    ) +
    scale_fill_manual(values = map_colors) +
    scale_x_continuous(
        labels = scales::label_number(scale = 1e-3, suffix = "k"),
        expand = expansion(mult = c(0, 0.1))
    ) +
    labs(x = "Samples", y = "Region") +
    theme_bw() +
    theme(
        # 全局文字设置
        text = element_text(size = 12),  # 所有文字统一12pt
        
        legend.position = "None",
        axis.title.x = element_text(face = "bold"),  # X轴标题加粗
        axis.title.y = element_text(face = "bold"),  # Y轴标题加粗
        
        # 轴标签设置
        axis.text.y = element_blank(),   # 隐藏Y轴标签
        axis.text.x = element_text(
            margin = margin(t = 2)         # 增加顶部间距
        ),
        axis.ticks.y = element_blank(),  # 隐藏Y轴刻度线
    )

# Region shannon 统计
region_shannon <- meta_data %>%
    select(region, shannon) %>%
    group_by(region) %>% filter(region %in% region)

region_shannon$region <- factor(region_shannon$region, levels = rev(region_name))

# 绘制箱线图
region_shannon <- na.omit(region_shannon)
region_shannon <- region_shannon[region_shannon$region != "Undefined",]
region_shannon_plot <- ggplot(region_shannon, aes(y = region, x = shannon)) +
    geom_violin(aes(fill = region),draw_quantiles = c( 0.5)) +
    stat_summary(fun= mean, geom = "point",
                 shape = 19, size = 2, color = "black") +
    scale_fill_manual(values = map_colors) +
    labs(y = NULL, x = "Shannon diversity") +
    theme_bw() +
    theme(
        text = element_text(size = 12),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none"
    )
region_shannon_plot


# Region reads 统计
region_reads <- meta_data %>%
    select(region, read_depth) %>%
    group_by(region) %>% filter(region %in% region)
region_reads$region <- factor(region_reads$region, levels = rev(region_name))
region_reads <- na.omit(region_reads)
region_reads_plot <- ggplot(region_reads, aes(y = region, x = read_depth)) +  # +1避免log(0)
    geom_violin(aes(fill = region), draw_quantiles = 0.5) +
    stat_summary(
        fun = mean, 
        geom = "point",
        shape = 19, 
        size = 2, 
        color = "black"
    ) +
    scale_fill_manual(values = map_colors) +
    scale_x_log10(
        name = "Read depth",
        labels = scales::trans_format("log10", scales::math_format(10^.x)),  # 指数格式
        breaks = c(1, 10, 1e2, 1e3, 1e4,1e5,1e6,1e7)  # 自定义刻度
    ) +
    labs(y = NULL) +
    theme_bw() +
    theme(
        text = element_text(size = 12),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y  = element_blank(),
        legend.position = "none"
    )

region_reads_plot

# 合并图形
merge_plot <- temp + region_counts_plot + region_shannon_plot + region_reads_plot +
    plot_layout(ncol = 4,widths = c(1, 0.5, 0.5, 0.5)) 
merge_plot


# 测序区间统计
reads_region <- meta_data %>%
    select(seq_region) 
unique_region <- unique(reads_region)
rownames(unique_region) <- 1 : nrow(unique_region)
for (i in c(4,5,6,8,13,20)){
    reads_region[reads_region == unique_region$seq_region[i]] = "V3-V4"
}
for (i in c(7,9)){
    reads_region[reads_region == unique_region$seq_region[i]] = "V1-V2"
}
for (i in c(11,19)){
    reads_region[reads_region == unique_region$seq_region[i]] = "V1-V3"
}
for (i in c(12,16,18,21)){
    reads_region[reads_region == unique_region$seq_region[i]] = "V4-V5"
}
for (i in c(14,17)){
    reads_region[reads_region == unique_region$seq_region[i]] = "V4"
}
for (i in c(14,17)){
    reads_region[reads_region == unique_region$seq_region[i]] = "V4"
}
for (i in c(23,24)){
    reads_region[reads_region == unique_region$seq_region[i]] = "V3-V5"
}
reads_region$seq_region[reads_region$seq_region == ""] <-"Unknown"
reads_region$seq_region[reads_region$seq_region == "multiple"] <-"Unknown"
reads_region$seq_region[reads_region$seq_region == "V3-V4,V5-V6"] <-"Unknown"

reads_region_count <- reads_region %>%
    group_by(seq_region) %>%
    summarise(count = n()) %>%
    arrange(desc(count))
reads_region_count$seq_region <- factor(reads_region_count$seq_region, 
                                        levels = c(unique(reads_region_count$seq_region[reads_region_count$seq_region != "Unknown"]),"Unknown"),order = TRUE)
#pie
region_colors <- c(brewer.pal(10, c("Paired")), brewer.pal(5, c("Set2")))
names(region_colors) <- unique(reads_region_count$seq_region)

# 生成带数量标签的饼图
# 计算标签位置（需先排序数据）
reads_region_count <- reads_region_count %>%
    arrange(desc(count)) %>%
    mutate(
        percentage = count / sum(count),
        cumulative = cumsum(percentage),
        angle = 90 - 360 * (cumulative - percentage/2)  # 计算标签角度（从12点方向开始）
    )

# 绘制环形图
region_pie <- ggplot(reads_region_count, aes(x = 2, y = count, fill = seq_region)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y", start = 0) +
    # 中心文字
    annotate("text", x = 0, y = 0, label = "ALL", size = 8, color = "black") +
    scale_fill_manual(values = region_colors, name=NULL) +
    theme_void(base_size=14) +
    guides(fill = guide_legend(ncol = 2)) +
    theme(text = element_text(face = "bold"))

region_pie
upp_plot <- plot_grid(p1, region_pie, p2, p3,
                         labels = c('a', 'b', 'c', 'd'),
                         align="hv",
                         scale = c(1, 1, 1, 1),
                         nrow = 1, ncol=4, plot=FALSE)

middle_plot <- plot_grid(p4, p5, p6, p7,
                         labels = c('e', 'f', 'g', 'h'),
                         align="hv",
                         scale = c(1, 1, 1, 1),
                         nrow = 1, ncol=4, plot=FALSE)

bottom_plot <- plot_grid(p8, merge_plot,
                         labels = c('i', 'j'),
                         align="hv",
                         scale = c(1, 1, 1, 1),
                         nrow = 1, ncol=2, plot=FALSE, rel_widths = c(0.25, 1))

p <- plot_grid(upp_plot, middle_plot, bottom_plot,
               align="hv",
               scale = c(1, 1, 1),
               nrow = 3, ncol = 1, plot=FALSE)

ggsave("/home/dongbiao/word_embedding_microbiome/result/gut_pretrain_data.png", p,
       width = 45, height = 25, units = "cm")


### beta diversity
map_colors <- c(
    "Europe and Northern America" = brewer.pal(7, "Set3")[1],
    "Eastern and South-Eastern Asia" = brewer.pal(7, "Set3")[2],
    "Sub-Saharan Africa" = brewer.pal(7, "Set3")[3],
    "Central and Southern Asia" = brewer.pal(7, "Set3")[4],
    "Australia/New Zealand" = brewer.pal(7, "Set3")[5],
    "Northern Africa and Western Asia" = brewer.pal(7, "Set3")[6],
    "Latin America and the Caribbean" = brewer.pal(7, "Set3")[7]
)

braytsne <- read.csv("Data/bray_TSNE.csv",header = T)
metadata <- read.csv("Data/metadata_filter.txt" ,header = T,"\t")
plot_data <- merge(braytsne, metadata, by.x = "X", by.y = "sample_id")
plot_data <- plot_data %>% mutate(region = case_when(region == "Unknown"~"unknown",
                                                     TRUE ~region))


#country散点图 p1
region_freq <- plot_data %>% group_by(region) %>% 
    summarise(count = n()) %>% arrange(desc(count))

# 获取排序后的区域列表
region_levels <- region_freq$region

# 将 region 列转换为因子，确保排序顺序
plot_data$region <- factor(plot_data$region, levels = region_levels)

# 按 region 排序数据框，确保最多样本的区域点先绘制（底层）
plot_data <- plot_data %>% arrange(region)

# 绘制散点图
p1 = ggplot(plot_data, aes(x = tsne1, y = tsne2, color = region)) + 
    geom_point(size = 0.5, show.legend = F, alpha=0.1) + 
    labs(x = "t-SNE 1", y = "t-SNE 2", title = "Region") + 
    theme_bw(base_size = 14) +
    scale_color_manual(values = map_colors) +
    theme(
        text = element_text(face = "bold"), 
        panel.border = element_rect(color = "black",               
            fill = NA, linewidth = 1),
        axis.line = element_line(color = "black", 
            linewidth = 0.5
        ),
        axis.ticks = element_line(color = "black",linewidth = 0.5
        )
    )

# p3-p9
#按区域分数据
en <- plot_data %>% filter(region=='Europe and Northern America')
es <- plot_data %>% filter(region=='Eastern and South-Eastern Asia')
sa <- plot_data %>% filter(region=='Sub-Saharan Africa')
cs <- plot_data %>% filter(region=='Central and Southern Asia')
az <- plot_data %>% filter(region=='Australia/New Zealand')
nw <- plot_data %>% filter(region=='Northern Africa and Western Asia')
lc <- plot_data %>% filter(region=='Latin America and the Caribbean')




# 创建遍历列表
data_list <- list(en, es, sa, cs, az, nw, lc)
p_vector <- paste0("p", 3:9)
limits <- c(1500, 900, 550, 320, 350, 190, 240)
title = c("Europe and Northern America",
          "Eastern and South-Eastern Asia",
          "Sub-Saharan Africa",
          "Central and Southern Asia",
          "Australia/New Zealand",
          "Northern Africa and Western Asia",
          "Latin America and the Caribbean")

# 创建一个包含标题、数据和限制的列表
variables <- list(
    t = title, 
    data = data_list, 
    limit = limits
)

# 使用purrr::pmap遍历列表
plots <- pmap(variables, function(t, data, limit) {
    ggplot(plot_data, aes(x = tsne1, y = tsne2)) +
        geom_bin_2d(fill = "#D3D3D3", color = "#D3D3D3", bins = 50) +
        geom_bin_2d(data = data, show.legend = FALSE, bins = 50) + 
        scale_fill_viridis_c(
            option = "D",
            limits = c(0, limit),
            direction = 1
        ) +
        theme_minimal() +
        labs(x = "", y = "", title = t) +
        theme(
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            axis.line = element_line(color = "black", linewidth = 0.5),
            axis.ticks = element_line(color = "black", linewidth = 0.5)
        )
})

# 为每个图形命名
names(plots) <- p_vector

# 使用list2env将列表中的对象添加到全局环境中
list2env(plots, envir = .GlobalEnv)


# 步骤1：确保所有子图具有相同的图注空间和主题参数
# ---------------------------------------------------------
# 统一所有子图的主题设置（包括图注宽度）
uniform_theme <- theme(
    legend.position = "right",  # 强制图注统一在右侧
    legend.box.margin = margin(0, 0, 0, 10),  # 图注右侧留出10pt空间
    plot.margin = margin(5, 5, 5, 5)         # 统一子图边距
)

# 为所有子图应用统一主题（即使无图注也保留空间）
p1 <- p1 + uniform_theme
p3 <- p3 + uniform_theme
p4 <- p4 + uniform_theme
p5 <- p5 + uniform_theme
p6 <- p6 + uniform_theme
p7 <- p7 + uniform_theme
p8 <- p8 + uniform_theme
p9 <- p9 + uniform_theme

# 步骤2：组合3x3布局并锁定宽度比例
# ---------------------------------------------------------
combined_plot <- p1 + p3 + p4 + p5 + p6 + p7 + p8 + p9 +
    plot_layout(
        nrow = 2,  ncol = 4,
    )

# 测序区间
metadata <- read.csv("Data/metadata_filter.txt" ,header = T,"\t")
reads_region <- metadata %>%
    select(sample_id, seq_region) 
unique_region <- as.data.frame( unique(reads_region$seq_region)) 
names(unique_region) = "seq_region"
rownames(unique_region) <- 1 : nrow(unique_region)
for (i in c(4,5,6,8,13,20)){
    reads_region[reads_region == unique_region$seq_region[i]] = "V3-V4"
}
for (i in c(7,9)){
    reads_region[reads_region == unique_region$seq_region[i]] = "V1-V2"
}
for (i in c(11,19)){
    reads_region[reads_region == unique_region$seq_region[i]] = "V1-V3"
}
for (i in c(12,16,18,21)){
    reads_region[reads_region == unique_region$seq_region[i]] = "V4-V5"
}
for (i in c(14,17)){
    reads_region[reads_region == unique_region$seq_region[i]] = "V4"
}
for (i in c(14,17)){
    reads_region[reads_region == unique_region$seq_region[i]] = "V4"
}
for (i in c(23,24)){
    reads_region[reads_region == unique_region$seq_region[i]] = "V3-V5"
}
reads_region$seq_region[reads_region$seq_region == ""] <-"Unknown"
reads_region$seq_region[reads_region$seq_region == "multiple"] <-"Unknown"
reads_region$seq_region[reads_region$seq_region == "V3-V4,V5-V6"] <-"Unknown"

reads_region <- merge(reads_region, braytsne,by.x = "sample_id", by.y = "X")
reads_region$seq_region <- factor(reads_region$seq_region, 
                                        levels = c(unique(reads_region$seq_region[reads_region$seq_region!= "Unknown"]),"Unknown"))

region_colors <- c(brewer.pal(10, c("Paired")), brewer.pal(7, c("Set2")))
reads_region <- reads_region %>% arrange(reads_region$seq_region)


# 生成 p2（注意修正未完成的 guides 参数）
p2 <- ggplot(reads_region, aes(x = tsne1, y = tsne2, color = seq_region)) + 
    geom_point(size = 0.5, alpha = 0.1) + 
    labs(x = "t-SNE 1", y = "t-SNE 2", title = "Sequence region") + 
    theme_bw(base_size = 14) +
    scale_color_manual(values = region_colors) +
    theme(
        text = element_text(face = "bold"), 
        legend.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        legend.position = "none"
    )

# 生成分面图的五个子图
filtered_data <- reads_region %>%
    filter(seq_region %in% c("V3-V4", "V3-V5", "V4-V5", "V1-V2", "V4"))

# 按分面顺序拆分数据（确保顺序与原分面一致）
split_data <- filtered_data %>%
    mutate(seq_region = factor(seq_region, levels = c("V3-V4", "V3-V5", "V4-V5", "V1-V2", "V4"))) %>%
    group_by(seq_region) %>%
    group_split()

# 生成每个子图并调整主题
plot_list <- lapply(split_data, function(sub_data) {
    ggplot(reads_region, aes(x = tsne1, y = tsne2)) +
        geom_bin_2d(fill = "#D3D3D3", color = "#D3D3D3", bins = 50) +
        geom_bin_2d(data = sub_data,show.legend = FALSE, bins = 50) + 
        scale_fill_viridis_c(option = "D", direction = 1) +
        labs(title = unique(sub_data$seq_region)) +
        theme_minimal() +
        theme(
            panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
            axis.line = element_line(color = "black", linewidth = 0.5),
            axis.ticks = element_line(color = "black", linewidth = 0.5),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.title = element_blank(),  # 隐藏坐标轴标题
            legend.position = "none"
        )
})

# 使用 patchwork 组合图形
combined_plot <- p2 + plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + plot_list[[4]] + plot_list[[5]] +
    plot_layout(nrow = 2, ncol = 3) 

# Load the data
pabundance <- read.csv("Data/phyla_abundance_top.csv", header = TRUE)
metadatashannon <- read.csv("Data/meta_data_include_shannon.csv", header = TRUE)
metadata <- read.csv("Data/metadata_filter.txt", header = TRUE,sep = "\t")
ms <- metadatashannon[,!names(metadatashannon) %in% names(metadata)]

df1 <- cbind(metadata, ms)
df2 <- merge(df1, braytsne, by.x = "sample_id", by.y = "X")
df3 <- merge(df2, pabundance, by.x = "sample_id", by.y = "X")
df4 <- df3 %>% pivot_longer(cols = c("shannon","Bacillota","Pseudomonadota","Actinomycetota","Bacteroidota","Thermodesulfobacteriota"),
                            names_to = "facet_variable",      
                            values_to = "value")
split_data <- split(df4, df4$facet_variable)
split_data <- split_data[c("shannon","Bacillota","Pseudomonadota","Actinomycetota","Bacteroidota","Thermodesulfobacteriota")]
sorted_list <- lapply(split_data, function(df) {
    df[order(df$value), ]
})
shannon_list <- sorted_list[["shannon"]]
sorted_list <- sorted_list[-which(names(sorted_list) == "shannon")]

# 散点图绘制

p1 = ggplot(shannon_list, aes(x = tsne1, y = tsne2)) +
    geom_point(aes(color = value), size = 0., alpha = 0.1, show.legend = T) +
    scale_color_viridis_c(
        option = "C",
        limits = range(shannon_list$value),  # 按当前分面数据范围设置颜色标尺
    ) +
    labs(x = "t-SNE 1", y = "t-SNE 2", title = "Shannon",color = "Shannon") +
    theme_minimal() +
    theme(
        panel.border = element_rect(      # 面板外框
            color = "black",                # 框线颜色
            fill = NA,                      # 填充透明
            linewidth = 1                   # 线宽（示例图较粗）
        ),
        axis.line = element_line(         # 坐标轴线
            color = "black", 
            linewidth = 0.5
        ),
        axis.ticks = element_line(        # 坐标轴刻度线
            color = "black",
            linewidth = 0.5
        )
    )


plot_list <- imap(sorted_list, ~ {
    # 判断是否为第一个子图（根据索引或名称调整条件）
    is_legend <- (.y == names(sorted_list)[1])  # 按列表第一个名称匹配
    # 或按位置判断：is_legend <- (.y == "1") （当列表无命名时）
    
    p <- ggplot(.x, aes(x = tsne1, y = tsne2)) +
        geom_point(
            aes(color = value),
            size = 0.,
            alpha = 0.1,
            show.legend = is_legend  # 仅目标图显示图例
        ) +
        scale_color_viridis_c(
            option = "D",
            limits = c(0, 100),
            guide = if (is_legend) "colourbar" else "none"  # 控制颜色标尺图例
        ) +
        labs(x = "t-SNE 1", y = "t-SNE 2", title = .y,color = "abundacne") +
        theme_minimal() +
        theme(
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            axis.line = element_line(color = "black", linewidth = 0.5),
            axis.ticks = element_line(color = "black", linewidth = 0.5)
        )
    return(p)
})


# 统一所有子图的主题设置（包括图注宽度）
uniform_theme <- theme(
    legend.position = "right",  # 强制图注统一在右侧
    legend.box.margin = margin(0, 0, 0, 10),  # 图注右侧留出10pt空间
    plot.margin = margin(5, 5, 5, 5)         # 统一子图边距
)  

p1 <- p1 + uniform_theme
plot_list <- lapply(plot_list, function(p) {
    p + uniform_theme
})

combined_plot <- (p1+plot_list[[1]]+plot_list[[2]])/(plot_list[[3]]+plot_list[[4]]+plot_list[[5]]) + plot_layout(
    widths = rep(1, 3), # 每列宽度相等
    guides = "collect"
    
)

