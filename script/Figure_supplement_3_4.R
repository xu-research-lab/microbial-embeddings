library(ade4)
library(ggplot2)
library(RColorBrewer)
library(vegan)
metadata<-read.delim('CRC/CRC/metadata.txt')
metadata$group[metadata$group==1]='CRC'
metadata$group[metadata$group==0]='Control'
result<-list()
for (i in 1:7) {
  tab<-read.csv(paste('CRC/CRC/fc_',i,'.csv',sep = ''),check.names = F)
  row.names(tab)<-tab[,1]
  tab<-tab[,-1]
  otu.distance <- vegdist(tab,method = 'euclidean')
  
  pcoa <- cmdscale (otu.distance,eig=TRUE)
  pc12 <- pcoa$points[,1:2]
  pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)#解释度
  pc12 <- as.data.frame(pc12)
  #给pc12添加samp1es变量
  pc12$Run <- row.names(pc12)
  
  df <- merge(pc12,metadata,by="Run")
  
  df$group<-factor(df$group,levels=c('Control','CRC'))
  
  
  colnames(df)[4]<-'Group'
  colnames(df)[5]<-'Study'
  
  
  
  Adonis <- adonis2(otu.distance~Study+Group, data=df, distance = "bray", permutations = 999)
  p1 <- df %>%
    ggplot(aes(x = V1, y = V2)) +
    geom_point(aes(shape = Group, color = Study), alpha = 0.8) +
    labs(x = paste0("PC1 ", pc[1], "%"),
         y = paste0("PC2 ", pc[2], "%"),title =paste('Model_',i,sep = '')) +
    stat_ellipse(geom = "polygon", level = 0.9,
                 linetype = 1, size = 0.8,
                 aes(fill = Group), alpha = 0.1,
                 show.legend = T) +
    theme_bw(base_size = 14) +
    theme(text = element_text(face = "bold")) +
    scale_color_brewer(palette = "Paired") +
    annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, 
             label = paste(
               "Group R:", round(Adonis$R2[2], 2), 
               ", P:", round(Adonis$`Pr(>F)`[2], 3)),
             size = 4)+
    theme(  
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold") # 设置标题居中，字体大小，加粗，颜色  
    )  
  
  
  result[[i+2]]<-p1
}



tab<-read.csv('CRC/CRC/Embedding_Transformation.csv',check.names = F)
row.names(tab)<-tab[,1]
tab<-tab[,-1]
otu.distance <- vegdist(tab,method = 'euclidean')

pcoa <- cmdscale (otu.distance,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)#解释度
pc12 <- as.data.frame(pc12)
#给pc12添加samp1es变量
pc12$Run <- row.names(pc12)

df <- merge(pc12,metadata,by="Run")

df$group<-factor(df$group,levels=c('Control','CRC'))


colnames(df)[4]<-'Group'
colnames(df)[5]<-'Study'



Adonis <- adonis2(otu.distance~Study+Group, data=df, distance = "bray", permutations = 999)
p1 <- df %>%
  ggplot(aes(x = V1, y = V2)) +
  geom_point(aes(shape = Group, color = Study), alpha = 0.8) +
  labs(x = paste0("PC1 ", pc[1], "%"),
       y = paste0("PC2 ", pc[2], "%"),title ='Embedding Transformation') +
  stat_ellipse(geom = "polygon", level = 0.9,
               linetype = 1, size = 0.8,
               aes(fill = Group), alpha = 0.1,
               show.legend = T) +
  theme_bw(base_size = 14) +
  theme(text = element_text(face = "bold")) +
  scale_color_brewer(palette = "Paired") +
  annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, 
           label = paste(
             "Group R:", round(Adonis$R2[2], 2), 
             ", P:", round(Adonis$`Pr(>F)`[2], 3)),
           size = 4)+
  theme(  
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold") # 设置标题居中，字体大小，加粗，颜色  
  )  


result[[2]]<-p1




biom<-read.delim2('CRC/CRC/table_from_biom.txt')

row.names(biom)<-biom[,1]
biom<-biom[,-1]
biom<-as.data.frame(t(biom))
biom[] <- lapply(biom, function(x) as.numeric(x))

biom<-as.data.frame(t(biom))
biom<-rownames_to_column(biom)

biom<-biom %>%
  mutate(across(-rowname,~ round(.x / sum(.x), 10)))
row.names(biom)<-biom[,1]
biom<-biom[,-1]

biom<-as.data.frame(t(biom))
otu.distance <- vegdist(biom,method = 'euclidean')

pcoa <- cmdscale (otu.distance,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)#解释度
pc12 <- as.data.frame(pc12)
#给pc12添加samp1es变量
pc12$Run <- row.names(pc12)

df <- merge(pc12,metadata,by="Run")

df$group<-factor(df$group,levels=c('Control','CRC'))


colnames(df)[4]<-'Group'
colnames(df)[5]<-'Study'



Adonis <- adonis2(otu.distance~Study+Group, data=df, distance = "bray", permutations = 999)
p1 <- df %>%
  ggplot(aes(x = V1, y = V2)) +
  geom_point(aes(shape = Group, color = Study), alpha = 0.8) +
  labs(x = paste0("PC1 ", pc[1], "%"),
       y = paste0("PC2 ", pc[2], "%"),title ='Abundance table') +
  stat_ellipse(geom = "polygon", level = 0.9,
               linetype = 1, size = 0.8,
               aes(fill = Group), alpha = 0.1,
               show.legend = T) +
  theme_bw(base_size = 14) +
  theme(text = element_text(face = "bold")) +
  scale_color_brewer(palette = "Paired") +
  annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, 
           label = paste(
             "Group R:", round(Adonis$R2[2], 2), 
             ", P:", round(Adonis$`Pr(>F)`[2], 3)),
           size = 4)+
  theme(  
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold") # 设置标题居中，字体大小，加粗，颜色  
  )  


result[[1]]<-p1  

library(patchwork)#
combined_plot <- wrap_plots(result, ncol = 3, nrow = 3) +  
  plot_layout(guides = 'collect') +
  plot_annotation(   
    theme = theme(plot.title = element_text(size = 17, hjust = 0.5)))  
ggsave("CRC_1.pdf", plot = combined_plot, device = "pdf", width = 11, height = 10) 



metadata<-read.delim('CRC/IBD/IBD/metadata.txt')
metadata$group[metadata$group==1]='IBD'
metadata$group[metadata$group==0]='Control'
result_IBD<-list()
for (i in 1:7) {
  tab<-read.csv(paste('CRC/IBD/IBD/fc_',i,'.csv',sep = ''),check.names = F)
  row.names(tab)<-tab[,1]
  tab<-tab[,-1]
  otu.distance <- vegdist(tab,method = 'euclidean')
  
  pcoa <- cmdscale (otu.distance,eig=TRUE)
  pc12 <- pcoa$points[,1:2]
  pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)#解释度
  pc12 <- as.data.frame(pc12)
  #给pc12添加samp1es变量
  pc12$'sample' <- row.names(pc12)
  
  df <- merge(pc12,metadata,by="sample")
  df$'V4'<-pc
  df$group<-factor(df$group,levels=c('Control','IBD'))
  
  
  colnames(df)[7]<-'Group'
  colnames(df)[4]<-'Study'
  
  write.csv(df,paste('fc_',i,'.txt',sep = ''))
  
  
  Adonis <- adonis2(otu.distance~Study+Group, data=df, distance = "bray", permutations = 999)
  write.csv(Adonis,paste('Adonis_','fc_',i,'.txt',sep ='' ))
  
}





tab<-read.csv('CRC/IBD/IBD/Embedding_Transformation.csv',check.names = F)
row.names(tab)<-tab[,1]
tab<-tab[,-1]
otu.distance <- vegdist(tab,method = 'euclidean')

pcoa <- cmdscale (otu.distance,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)#解释度
pc12 <- as.data.frame(pc12)
#给pc12添加samp1es变量
pc12$'sample'<- row.names(pc12)

df <- merge(pc12,metadata,by="sample")
df$'V4'<-pc
df$group<-factor(df$group,levels=c('Control','IBD'))


colnames(df)[7]<-'Group'
colnames(df)[4]<-'Study'

write.csv(df,'Embedding_Transformation_IBD.txt')

Adonis <- adonis2(otu.distance~Study+Group, data=df, distance = "bray", permutations = 999)
write.csv(Adonis,'Adonis_Embedding_Transformation_IBD.txt')





biom<-read.delim2('CRC/IBD/IBD/table.from_biom.txt',check.names = F)

row.names(biom)<-biom[,1]
biom<-biom[,-1]
biom<-as.data.frame(t(biom))
biom[] <- lapply(biom, function(x) as.numeric(x))

biom<-as.data.frame(t(biom))
biom<-rownames_to_column(biom)

biom<-biom %>%
  mutate(across(-rowname,~ round(.x / sum(.x), 10)))
row.names(biom)<-biom[,1]
biom<-biom[,-1]

biom<-as.data.frame(t(biom))
otu.distance <- vegdist(biom,method = 'euclidean')

pcoa <- cmdscale (otu.distance,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)#解释度
pc12 <- as.data.frame(pc12)
#给pc12添加samp1es变量
pc12$'sample' <- row.names(pc12)

df <- merge(pc12,metadata,by="sample")
df$'V4'<-pc
df$group<-factor(df$group,levels=c('Control','IBD'))


colnames(df)[7]<-'Group'
colnames(df)[4]<-'Study'


write.csv(df,'biom_IBD.txt')

Adonis <- adonis2(otu.distance~Study+Group, data=df, distance = "bray", permutations = 999)

write.csv(Adonis,'Adonis_biom_IBD.txt')
#plot
for (i in 1:7) {
  df<-read.csv2(paste('fc_',i,'.txt',sep = ''),sep = ',',check.names = F)
  Adonis<-read.csv2(paste('Adonis_fc_',i,'.txt',sep = ''),sep = ',',check.names = F)
  df$V1<-as.numeric(df$V1)
  df$V2<-as.numeric(df$V2)
  Adonis$R2 <-as.numeric(Adonis$R2)
  Adonis$`Pr(>F)`<-as.numeric(Adonis$`Pr(>F)`)
  
  p1 <- df %>%
    ggplot(aes(x = V1, y = V2)) +
    geom_point(aes(shape = Group, color = Study), alpha = 0.8) +
    labs(x = paste0("PC1 ", df$V4[1], "%"),
         y = paste0("PC2 ", df$V4[2], "%"),title =paste('Model_',i,sep = '')) +
    stat_ellipse(geom = "polygon", level = 0.9,
                 linetype = 1, size = 0.8,
                 aes(fill = Group), alpha = 0.1,
                 show.legend = T) +
    theme_bw(base_size = 14) +
    theme(text = element_text(face = "bold")) +
    scale_color_brewer(palette = "Paired") +
    annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, 
             label = paste(
                           "Group R:", round(Adonis$R2[2], 2), 
                           ", P:", round(Adonis$`Pr(>F)`[2], 3)),
             size = 4)+
    theme(  
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold") # 设置标题居中，字体大小，加粗，颜色  
    )  
  
  
  result_IBD[[i+2]]<-p1
}


df<-read.csv2('Embedding_Transformation_IBD.txt',sep = ',',check.names = F)
Adonis<-read.csv2('Adonis_Embedding_Transformation_IBD.txt',sep = ',',check.names = F)

Adonis$R2<-as.numeric(Adonis$R2)
Adonis$`Pr(>F)`<-as.numeric(Adonis$`Pr(>F)`)

df$V1<-as.numeric(df$V1)
df$V2<-as.numeric(df$V2)
p1 <- df %>%
  ggplot(aes(x = V1, y = V2)) +
  geom_point(aes(shape = Group, color = Study), alpha = 0.8) +
  labs(x = paste0("PC1 ", df$V4[1], "%"),
       y = paste0("PC2 ", df$V4[2], "%"),title ='Abundance table') +
  stat_ellipse(geom = "polygon", level = 0.9,
               linetype = 1, size = 0.8,
               aes(fill = Group), alpha = 0.1,
               show.legend = T) +
  theme_bw(base_size = 14) +
  theme(text = element_text(face = "bold")) +
  scale_color_brewer(palette = "Paired") +
  annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, 
           label = paste(
             "Group R:", round(Adonis$R2[2], 2), 
             ", P:", round(Adonis$`Pr(>F)`[2], 3)),
           size = 4)+
  theme(  
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold") # 设置标题居中，字体大小，加粗，颜色  
  )  


result_IBD[[2]]<-p1  








df<-read.csv2('biom_IBD.txt',sep = ',',check.names = F)
Adonis<-read.csv2('Adonis_biom_IBD.txt',sep = ',',check.names = F)

Adonis$R2<-as.numeric(Adonis$R2)
Adonis$`Pr(>F)`<-as.numeric(Adonis$`Pr(>F)`)

df$V1<-as.numeric(df$V1)
df$V2<-as.numeric(df$V2)
p1 <- df %>%
  ggplot(aes(x = V1, y = V2)) +
  geom_point(aes(shape = Group, color = Study), alpha = 0.8) +
  labs(x = paste0("PC1 ", df$V4[1], "%"),
       y = paste0("PC2 ", df$V4[2], "%"),title ='Abundance table') +
  stat_ellipse(geom = "polygon", level = 0.9,
               linetype = 1, size = 0.8,
               aes(fill = Group), alpha = 0.1,
               show.legend = T) +
  theme_bw(base_size = 14) +
  theme(text = element_text(face = "bold")) +
  scale_color_brewer(palette = "Paired") +
  annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, 
           label = paste(
             "Group R:", round(Adonis$R2[2], 2), 
             ", P:", round(Adonis$`Pr(>F)`[2], 3)),
           size = 4)+
  theme(  
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold") # 设置标题居中，字体大小，加粗，颜色  
  )  


result_IBD[[1]]<-p1  

combined_plot <- wrap_plots(result_IBD, ncol = 3, nrow = 3) +  
  plot_layout(guides = 'collect') +
  plot_annotation(   
    theme = theme(plot.title = element_text(size = 17, hjust = 0.5)))  
ggsave("IBD_1.pdf", plot = combined_plot, device = "pdf", width = 11, height = 10) 

