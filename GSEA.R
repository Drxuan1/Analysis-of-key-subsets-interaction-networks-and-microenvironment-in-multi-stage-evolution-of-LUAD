library(Seurat)
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(aPEAR)
library(cowplot)

setwd("")

set.seed(2024) 

options(scipen=100)

pbmc<-readRDS("pbmcnew.rds")

logFCfilter=0         

adjPvalFilter=0.05      

pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]

pcDEG <- sig.markers

marker<-function(y){
  results<-list()
  for (i in names(table(y$cluster))) {
    result<-y[which(y$cluster==i),]
    results[[i]]<-result
  }
  return(results)
}

markers<-marker(pcDEG)

AT1<-as.data.frame(markers$AT1)

AT2<-as.data.frame(markers$AT2)

CLUB<-as.data.frame(markers$CLUB)

CILIATED<-as.data.frame(markers$CILIATED)

geneMap <- bitr(AT1$gene, 
                fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")
AT1 <- merge(AT1,geneMap,by.x = "gene",by.y = "SYMBOL")

geneMap <- bitr(AT2$gene, 
                fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")
AT2 <- merge(AT2,geneMap,by.x = "gene",by.y = "SYMBOL")

geneMap <- bitr(CLUB$gene, 
                fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")
CLUB <- merge(CLUB,geneMap,by.x = "gene",by.y = "SYMBOL")

geneMap <- bitr(CILIATED$gene, 
                fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")
CILIATED <- merge(CILIATED,geneMap,by.x = "gene",by.y = "SYMBOL")

AT1<-arrange(AT1, desc(AT1$avg_log2FC))

AT2<-arrange(AT2, desc(AT2$avg_log2FC))

CLUB<-arrange(CLUB, desc(CLUB$avg_log2FC))

CILIATED<-arrange(CILIATED, desc(CILIATED$avg_log2FC))

geneListAT1 = AT1$avg_log2FC

names(geneListAT1) = AT1$ENTREZID

geneListAT2 = AT2$avg_log2FC

names(geneListAT2) = AT2$ENTREZID

geneListCLUB = CLUB$avg_log2FC

names(geneListCLUB) = CLUB$ENTREZID

geneListCILIATED = CILIATED$avg_log2FC

names(geneListCILIATED) = CILIATED$ENTREZID
#
enrich1 <- gseGO(geneListAT1,OrgDb = org.Hs.eg.db, ont = 'BP',eps = 0,pvalueCutoff = 1)

dim(enrich1@result)

p1 <- enrichmentNetwork(enrich1@result[1:100,],repelLabels = TRUE, drawEllipses = TRUE);p1

data1 <- enrich1@result[1:100,]

p1 + scale_x_continuous(limits = c(-10, 10))+ 
  
  scale_y_continuous(limits = c(-10, 10))

saveRDS(enrich1,"enrich1_AT1_BP.rds")

saveRDS(data1, "data1.rds")

saveRDS(p1, "p1.rds")
#
enrich2 <- gseGO(geneListAT2,OrgDb = org.Hs.eg.db, ont = 'BP',eps = 0,pvalueCutoff = 1)

dim(enrich2@result)

p2<-enrichmentNetwork(enrich2@result[1:100,],repelLabels = TRUE, drawEllipses = TRUE);p2

data2 <- enrich2@result[1:100,]

p2 + scale_x_continuous(limits = c(-10, 10))+
  
  scale_y_continuous(limits = c(-10, 10))

saveRDS(enrich2,"enrich2_AT2_BP.rds")

saveRDS(data2, "data2.rds")

saveRDS(p2, "p2.rds")
#
enrich3 <- gseGO(geneListCLUB,OrgDb = org.Hs.eg.db, ont = 'BP',eps = 0,pvalueCutoff = 1)

dim(enrich3@result)

p3<-enrichmentNetwork(enrich3@result[1:100,],repelLabels = TRUE, drawEllipses = TRUE);p3

data3 <- enrich3@result[1:100,]
p3 + scale_x_continuous(limits = c(-10, 10))+
  
  scale_y_continuous(limits = c(-10, 10))

saveRDS(enrich3,"enrich3_CLUB_BP.rds")

saveRDS(data3, "data3.rds")

saveRDS(p3, "p3.rds")
#
enrich4 <- gseGO(geneListCILIATED,OrgDb = org.Hs.eg.db, ont = 'BP',eps = 0,pvalueCutoff = 1)

dim(enrich4@result)

p4<-enrichmentNetwork(enrich4@result[1:100,],repelLabels = TRUE, drawEllipses = TRUE);p4

data4 <- enrich4@result[1:100,]

p4 + scale_x_continuous(limits = c(-10, 10))+
  
  scale_y_continuous(limits = c(-10, 10))

saveRDS(enrich4,"enrich4_CILIATED_BP.rds")

saveRDS(data4, "data4.rds")

saveRDS(p4, "p4.rds")

plot_grid(p1,p2,p3,p4,ncol = 2)
