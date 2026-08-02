library(Seurat)
library(pheatmap)
library(raster)
library(corrplot)

setwd("")

pbmc<-readRDS("pbmcMACAT2CD4CD8CAF1B_END.rds")

pbmc$Group<-ifelse(substr(pbmc$orig.ident,10,12)=="AIS","AIS","IAC")

pbmcAIS<-subset(pbmc,Group %in% c("AIS"))

pbmcIAC<-subset(pbmc,Group %in% c("IAC"))

saveRDS(pbmcAIS,"pbmcAIS.rds")

saveRDS(pbmcIAC,"pbmcIAC.rds")
#
av.expAIS <- AverageExpression(pbmcAIS)$RNA#pbmcIAC

features <- names(tail(sort(apply(av.expAIS, 1, sd)),nrow(av.expAIS)))#av.expIAC

av.expAIS <- as.data.frame(av.expAIS[which(row.names(av.expAIS)%in% features),])#av.expIAC

av.expAIS <- cor(av.expAIS, method= "spearman")#av.expIAC

testRes <- cor.mtest(av.expAIS,conf.level=0.95)#av.expIAC

corrplot(av.expAIS, #av.expIAC
         method = "number", 
         type = "lower",
         tl.pos = "lt",
         tl.col = "black"
) 
corrplot(av.expAIS, #av.expIAC
         method = "color", 
         type = "upper",
         add = T,
         tl.pos = "n",
         cl.pos = "n",
         diag = F,
         p.mat = testRes$p,
         sig.level = c(0.001,0.01,0.05),
         pch.cex = 1.5,
         insig = "label_sig",
         tl.col = "black"
) 
