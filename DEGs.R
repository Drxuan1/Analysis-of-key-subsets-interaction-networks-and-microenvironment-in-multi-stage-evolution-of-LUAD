library(tidyverse)
library(Seurat)
library(ggplot2)
library(ggrepel)

setwd("")

pbmc<-readRDS("pbmcnew.rds")

logFCfilter=1       

adjPvalFilter=0.05    

pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
pbmc.markers <- pbmc.markers %>%mutate(type = ifelse(avg_log2FC >=1,"Up","Down"))%>%mutate(type2 = ifelse(p_val_adj < 0.05,"adjust Pvalue < 0.05","adjust Pvalue >= 0.05"))

write.csv(pbmc.markers,"pbmc.markers.csv")

pbmc.markers <- read.csv("pbmc.markers.csv", row.names = 1)

cell <-unique(pbmc.markers$cluster)

back.data<- data.frame()

for(n in 1:length(cell))
{
  tmp <- pbmc.markers %>%filter(cluster==cell[n])
  new.tmp <- data.frame(cluster = cell[n],min = min(tmp$avg_log2FC) - 0.2,max = max(tmp$avg_log2FC) + 0.2)
  back.data <- rbind(back.data,new.tmp)
}
#
cell <-unique(pbmc.markers$cluster)

up.top<- data.frame()

for(n in 1:length(cell))
{
  tmp <- pbmc.markers %>% filter(cluster==cell[n]) %>% filter(avg_log2FC>0)%>% arrange(desc(avg_log2FC)) %>%head(5)
  up.top <- rbind(up.top,tmp)
}

down.top<- data.frame()

for(n in 1:length(cell)){
  tmp <- pbmc.markers %>% filter(cluster==cell[n]) %>% filter(avg_log2FC<0)%>% arrange(avg_log2FC) %>%head(5)
  down.top <- rbind(down.top,tmp)
}
#
ggplot(pbmc.markers, aes(cluster, avg_log2FC)) +
  geom_jitter(aes(color = type)) +
  geom_col(data = back.data,aes(x = cluster,y = min),fill="grey93",alpha=0.5) +
  geom_col(data = back.data,aes(x = cluster,y = max),fill="grey93",alpha=0.5) +
  scale_color_manual(values=c(Down="#0099CC",Up="#CC3333"))+
  theme_classic(base_size = 24) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.7,0.9),
        legend.title = element_blank(),
        legend.background = element_blank()) +
  xlab('Clusters') + ylab('Average log2FoldChange') +
  guides(color = guide_legend(override.aes = list(size = 8)))+
  geom_tile(aes(x = cluster,y = 0,fill = cluster),color = 'black',height = 1,alpha = 0.3,show.legend = F)+
  scale_fill_manual(values = c("AT2" = "red","AT1" = "#C73E3A",  "CLUB" = "#F46D43","CILIATED" = "#A35E47"))+ 
  geom_text(data=back.data,aes(x = cluster,y = 0,label = cluster),size=12,color="white") +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  geom_text_repel(data = up.top,aes(x = cluster,y = avg_log2FC,label = gene),max.overlaps = 50, size = 8)+
  geom_text_repel(data = down.top,aes(x = cluster,y = avg_log2FC,label = gene),max.overlaps = 50, size = 8)
#
