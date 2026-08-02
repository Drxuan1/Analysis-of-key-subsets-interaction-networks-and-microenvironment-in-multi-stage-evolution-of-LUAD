library(phangorn)
library(scMetabolism)
library(ggplot2)
library(rsvd)
library(AUCell)
library(GSEABase)
library(GSVA)

setwd("")

pbmc <- readRDS("pbmcnew.rds")

if (!require('R.utils')) install.packages('R.utils')
R.utils::setOption( "clusterProfiler.download.method",'auto' )

pbmc@meta.data$ident<-pbmc@active.ident

countexp.Seurat <- sc.metabolism.Seurat(obj = pbmc, 
                                      method = "AUCell", #VISION, AUCell, ssgsea，gsva
                                      imputation = F, ncores = 2, 
                                      metabolism.type = "KEGG")

saveRDS(countexp.Seurat, "countexp.Seurat.rds")

metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score

#MedBioInfoCloud: rownames(countexp.Seurat@assays[["METABOLISM"]][["score"]])[1:6]

DimPlot.metabolism(obj = countexp.Seurat, 
                   pathway = "Pentose phosphate pathway", 
                   dimention.reduction.type = "umap", 
                   dimention.reduction.run = F, size = 1)

DimPlot.metabolism(obj = countexp.Seurat, 
                   pathway = "Citrate cycle (TCA cycle)", 
                   dimention.reduction.type = "umap", 
                   dimention.reduction.run = F, size = 1)

input.pathway<-rownames(countexp.Seurat@assays[["METABOLISM"]][["score"]])[1:20]

DotPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, 
                   phenotype = "ident", norm = "y")

