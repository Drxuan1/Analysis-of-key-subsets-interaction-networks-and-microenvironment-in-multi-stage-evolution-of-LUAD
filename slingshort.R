library(slingshot)
library(Seurat)
library(devtools)
library(ggplot2)
library(Matrix)
library(dplyr)
library(RColorBrewer)
library(DelayedMatrixStats)
library(scales)
library(paletteer) 
library(viridis)

setwd("")

pbmc <- readRDS("pbmcnew.rds")

pbmc$celltype<-pbmc@active.ident

pbmc$group <- ifelse(substr(pbmc@meta.data$orig.ident,10,12) =="AIS" , "AIS", "IAC")

sce <- as.SingleCellExperiment(pbmc, assay = "RNA")

sce_slingshot1 <- slingshot(sce,     
                            reducedDim = 'UMAP',  
                            clusterLabels = sce$celltype)

SlingshotDataSet(sce_slingshot1) 

cell_pal <- function(cell_vars, pal_fun, ...) {
  if (is.numeric(cell_vars)) {
    if (is.function(pal_fun)) {
      pal <- pal_fun(100, ...)
      return(pal[cut(cell_vars, breaks = 100, include.lowest = TRUE)])
    } else {
      stop("For numeric data, pal_fun must be a color palette function")
    }
  } else {
    categories <- sort(unique(cell_vars))
    n_categories <- length(categories)
    if (is.function(pal_fun)) {
      pal <- setNames(pal_fun(n_categories, ...), categories)
    } else {
      color_vector <- pal_fun
      if (n_categories <= length(color_vector)) {
        pal <- setNames(color_vector[1:n_categories], categories)
      } else {
        pal <- setNames(rep(color_vector, ceiling(n_categories / length(color_vector)))[1:n_categories], categories)
      }
    }
    return(pal[cell_vars])
  }
}

cell_colors <- cell_pal(sce_slingshot1$celltype, c( "#F8766D", "#7CAE00","#00BFC4",  "#C77CFF"))

plot(reducedDims(sce_slingshot1)$UMAP, col = cell_colors, pch=16, asp = 1, cex = 0.6)

lines(SlingshotDataSet(sce_slingshot1), lwd=3, col='black')

celltype_label <- pbmc@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(celltype = pbmc@meta.data$celltype) %>%
  group_by(celltype) %>%
  summarise(UMAP1 = median(UMAP_1),
            UMAP2 = median(UMAP_2))

for (i in 1:8) {
  text(celltype_label$celltype[i], x=celltype_label$UMAP1[i], y=celltype_label$UMAP2[i]+2)
}
