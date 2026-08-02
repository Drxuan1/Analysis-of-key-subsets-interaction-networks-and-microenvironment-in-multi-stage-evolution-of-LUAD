#devtools::install_github("campbio/decontX")
library(Seurat)
library(SingleR)
library(stringr)
library(dplyr)
library(magrittr)
library(umap)
library(ggplot2)
library(decontX)

setwd("")

pbmc<-readRDS("pbmc.rds")

counts <- pbmc@assays$RNA@counts

decontX_results <- decontX(counts)

pbmc$Contamination <- decontX_results$contamination

FeaturePlot(pbmc, 
            features = 'Contamination', 
            raster=FALSE     
) + 
  scale_color_viridis_c()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  xlab('scVI_UMAP_1')+
  ylab('scVI_UMAP_2')

low_con_scobj <- pbmc[,pbmc$Contamination < 0.20] 

UMAPPlot(object <- pbmc,pt.size = 0.01, label = TRUE,label.box=F)  

UMAPPlot(object <- low_con_scobj,pt.size = 0.01, label = TRUE,label.box=F)   

saveRDS(low_con_scobj,"low_con_scobj.rds")

table(low_con_scobj@active.ident)
#
counts1 <- low_con_scobj@assays$RNA@counts

low_con_scobj <- CreateSeuratObject(counts = counts1,project = "seurat", min.cells=3, min.features=50, names.delim = "_")

low_con_scobj[["percent.mt"]] <- PercentageFeatureSet(object = low_con_scobj, pattern = "^MT-")

pdf(file="01.featureViolin.pdf", width=10, height=6)

VlnPlot(object = low_con_scobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

dev.off()

low_con_scobj <- subset(x = low_con_scobj, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt <30)  

data30 <- as.data.frame(table(low_con_scobj@meta.data$orig.ident))

rownames(data30) <- data30$Var1

pdf(file="01.featureCor.pdf",width=10,height=6)

plot1 <- FeatureScatter(object = low_con_scobj, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)

plot2 <- FeatureScatter(object = low_con_scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)

CombinePlots(plots = list(plot1, plot2))

dev.off()

low_con_scobj <- NormalizeData(object = low_con_scobj, normalization.method = "LogNormalize", scale.factor = 10000)

low_con_scobj <- FindVariableFeatures(object = low_con_scobj, selection.method = "vst", nfeatures = 2000)

top10 <- head(x = VariableFeatures(object = low_con_scobj), 10)

pdf(file="01.featureVar.pdf",width=10,height=6)

plot1 <- VariableFeaturePlot(object = low_con_scobj)

plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

CombinePlots(plots = list(plot1, plot2))

dev.off()
#
low_con_scobj <- ScaleData(low_con_scobj)          

low_con_scobj <- RunPCA(object= low_con_scobj,npcs = 50,pc.genes=VariableFeatures(object = low_con_scobj))     

pdf(file="02.pcaGene.pdf",width=10,height=8)

VizDimLoadings(object = low_con_scobj, dims = 1:4, reduction = "pca",nfeatures = 20)

dev.off()

pdf(file="02.PCA.pdf",width=6.5,height=6)

DimPlot(object = low_con_scobj, reduction = "pca")

dev.off()

pdf(file="02.pcaHeatmap.pdf",width=10,height=8)

DimHeatmap(object = low_con_scobj, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)

dev.off()

low_con_scobj <- JackStraw(object = low_con_scobj, num.replicate = 100)

low_con_scobj <- ScoreJackStraw(object = low_con_scobj, dims = 1:20)

pdf(file="02.pcaJackStraw.pdf",width=8,height=6)

JackStrawPlot(object = low_con_scobj, dims = 1:20)

dev.off()

ElbowPlot(low_con_scobj,ndims=50)
#
pcSelect=19

low_con_scobj <- FindNeighbors(object = low_con_scobj, dims = 1:pcSelect)     

low_con_scobj <- FindClusters(object = low_con_scobj, resolution = 1)   

low_con_scobj <- RunUMAP(object = low_con_scobj, dims = 1:pcSelect)

pdf(file="03.UMAP.pdf",width=6.5,height=6)

UMAPPlot(object = low_con_scobj, pt.size = 0.005, label = TRUE)    

dev.off()

#write.table(pbmc$seurat_clusters,file="03.tsneCluster.txt",quote=F,sep="\t",col.names=F)

logFCfilter=2 

adjPvalFilter=0.05

low_con_scobj.markers <- FindAllMarkers(object = low_con_scobj,
                                        only.pos = FALSE,
                                        min.pct = 0.25,
                                        logfc.threshold = logFCfilter)
sig.markers <- low_con_scobj.markers[(abs(as.numeric(as.vector(low_con_scobj.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(low_con_scobj.markers$p_val_adj))<adjPvalFilter),]

write.table(sig.markers,file="03.clusterMarkers.txt",sep="\t",row.names=F,quote=F)

top10 <- low_con_scobj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

pdf(file="03.tsneHeatmap.pdf",width=12,height=9)

DoHeatmap(object = low_con_scobj, features = top10$gene) + NoLegend()

dev.off()

showGenes=c("SPP1")

pdf(file="03.markerViolin.pdf",width=10,height=6)

VlnPlot(object = low_con_scobj, features = showGenes)

dev.off()

showGenes=c("SPP1") 

pdf(file="03.markerScatter.pdf",width=10,height=6)

FeaturePlot(object = low_con_scobj, features = showGenes, cols = c("green", "red"))

dev.off()

pdf(file="03.markerBubble.pdf",width=12,height=6)

cluster10Marker=showGenes

DotPlot(object = low_con_scobj, features = cluster10Marker)

dev.off()
#
counts <- low_con_scobj@assays$RNA@counts

clusters <- low_con_scobj@meta.data$seurat_clusters

#ref=get(load("ref_Human_all.RData"))

ref <- celldex::HumanPrimaryCellAtlasData()

singler <- SingleR(test=counts, ref =ref,
                labels=ref$label.main, clusters = clusters)#(labels=ref$label.fine)

clusterAnn <- as.data.frame(singler)

clusterAnn <- cbind(id=row.names(clusterAnn), clusterAnn)

clusterAnn <- clusterAnn[,c("id", "labels")]

write.table(clusterAnn,file="04.clusterAnn.txt",quote=F,sep="\t", row.names=F)

singler2 <- SingleR(test=counts, ref =ref, 
                 labels=ref$label.main)#(labels=ref$label.fine)

cellAnn <- as.data.frame(singler2)

cellAnn <- cbind(id=row.names(cellAnn), cellAnn)

cellAnn <- cellAnn[,c("id", "labels")]

write.table(cellAnn, file="04.cellAnn.txt", quote=F, sep="\t", row.names=F)

newLabels <- singler$labels

names(newLabels) <- levels(low_con_scobj)

low_con_scobj <- RenameIdents(low_con_scobj, newLabels)

#pdf(file="04.TSNE.pdf",width=10,height=10)

#TSNEPlot(object = pbmc, pt.size = 0.005, label = TRUE)    

#dev.off()

pdf(file="04.UMAP.pdf",width=10,height=10)

UMAPPlot(object = low_con_scobj, pt.size = 0.005, label = TRUE)   

dev.off()

saveRDS(low_con_scobj,"low_con_scobjend.rds")

table(low_con_scobj@active.ident)
#
low_con_scobj <- readRDS("low_con_scobjend.rds")

low_con_scobj$celltype <- low_con_scobj$seurat_clusters

table(low_con_scobj$celltype)

genes_to_check <- c("CD3D","CD3E","CD8A",#T cells
                   "KLRD1","KLRF1","FGFBP2",#NK cells
                   "CD19","CD79A","MS4A1","JCHAIN",#B cells
                   "EPCAM","KRT19","SCGB3A1","SCGB3A2",#epithelial cells
                   "CD68","CD86","CD80","CD163","MRC1","CD14","CSF1R",#macrophage cells
                   "VWF","PECAM1","CLDN5","RAMP2",#endothelial cells
                   "FGF7","COL1A2","DCN","ACTA2",#fibroblast cells
                   "S100A8","S100A9","FCGR3B","MME"#neutrophil cells
)

genes_to_check <- str_to_upper(unique(genes_to_check))

genes_to_check

th <- theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))

p_all_markers <- DotPlot(low_con_scobj, features = genes_to_check,
                         assay='RNA' ,group.by = 'celltype' ,cols = c("black","#E54924"))  + coord_flip()+th
p_all_markers

data <- p_all_markers$data

colnames(data)

colnames(data) <- c("AverageExpression_unscaled","Precent Expressed","Features","celltype","Average Expression")

unique(data$`Precent Expressed`)
#
p <- ggplot(data,aes(celltype,Features,size = `Precent Expressed` ))+
  geom_point(shape=21,aes(fill= `Average Expression`),position =position_dodge(0))+
  theme_minimal()+xlab(NULL)+ylab(NULL) +
  scale_size_continuous(range=c(1,10))+theme_bw()+
  scale_fill_gradient(low = "#498EA4", high = "#E54924")+
  theme(legend.position = "right",legend.box = "vertical",
        legend.margin=margin(t= 0, unit='cm'),
        legend.spacing = unit(0,"in"),
        axis.text.x  = element_text(color="black",size=16,angle = 45, 
                                    vjust = 0.5, hjust=0.5),
        axis.text.y  = element_text(color="black",size=12),
        legend.text = element_text(size =12,color="black"),
        legend.title = element_text(size =12,color="black"),
        axis.title.y=element_text(vjust=1,  
                                  size=18),
        element_blank()
  )+labs(x="Celltype",y = "Features");p
#
low_con_scobj@active.ident <- low_con_scobj$seurat_clusters

table(low_con_scobj@active.ident)

newLabels <- c("T cells","B cells","T cells","T cells","T cells","NK cells","T cells","T cells","Macrophage cells","NK cells","Macrophage cells",#0-10
            "Epithelial cells","B cells","Epithelial cells","Macrophage cells","T cells","B cells","Epithelial cells","Epithelial cells",#11-18
            "Macrophage cells","Fibroblast cells","T cells","Neutrophil cells","Endothelial cells","Epithelial cells","Epithelial cells","Double cells","Epithelial cells","Fibroblast cells","Epithelial cells","T cells",#19-30
            "B cells")#31

names(newLabels) <- levels(low_con_scobj)

low_con_scobj <- RenameIdents(low_con_scobj, newLabels)

low_con_scobj$celltype <- low_con_scobj@active.ident

saveRDS(low_con_scobj,"low_con_scobjendann1.rds")

low_con_scobj <- readRDS("low_con_scobjendann1.rds")

low_con_scobj <- subset(low_con_scobj, celltype %in% c("T cells","B cells","T cells","T cells","T cells","NK cells","T cells","T cells","Macrophage cells","NK cells","Macrophage cells",
                                                       "Epithelial cells","B cells","Epithelial cells","Macrophage cells","T cells","B cells","Epithelial cells","Epithelial cells",
                                                       "Macrophage cells","Fibroblast cells","T cells","Neutrophil cells","Endothelial cells","Epithelial cells","Epithelial cells","Epithelial cells","Fibroblast cells","Epithelial cells","T cells",
                                                       "B cells"))

low_con_scobj$celltype <- low_con_scobj@active.ident

saveRDS(low_con_scobj,"low_con_scobjendann2.rds")

low_con_scobj$celltype <- low_con_scobj@active.ident

genes_to_check

th <- theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))

p_all_markers <- DotPlot(low_con_scobj, features = genes_to_check,
                         assay='RNA' ,group.by = 'celltype' ,cols = c("black","#E54924"))  + coord_flip()+th
p_all_markers

data <- p_all_markers$data

colnames(data)

colnames(data) <- c("AverageExpression_unscaled","Precent Expressed","Features","celltype","Average Expression")

unique(data$`Precent Expressed`)

celltypes <- c("T cells", "NK cells", "B cells", "Macrophage cells", "Neutrophil cells", "Epithelial cells", "Endothelial cells", "Fibroblast cells")

data$celltype <- factor(data$celltype, levels = celltypes)

p1 <- ggplot(data,aes(celltype, Features,size = `Precent Expressed` ))+
  geom_point(shape=21,aes(fill= `Average Expression`),position =position_dodge(0))+
  theme_minimal()+xlab(NULL)+ylab(NULL) +
  scale_size_continuous(range=c(1,10))+theme_bw()+
  scale_fill_gradient2(low = "#4682B4", mid = "grey",high = "#E54924")+
  theme(legend.position = "right",legend.box = "vertical",
        legend.margin=margin(t= 0, unit='cm'),
        legend.spacing = unit(0,"in"),
        axis.text.x  = element_text(color="black",size=24,angle = 45, 
                                    vjust = 0.5, hjust=0.5),
        axis.text.y  = element_text(color="black",size=24),
        legend.text = element_text(size =24,color="black"),
        legend.title = element_text(size =24,color="black"),
        axis.title.y=element_text(vjust=1,  
                                  size=24),
        element_blank()
  )+labs(x="Celltype",y = "Features");p1

