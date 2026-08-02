library(Seurat)
library(stringr)
library(dplyr)
library(magrittr)
library(monocle)
library(umap)
library(ggsci)
library(data.table)
library(ggpubr)

setwd("")

pbmc <- readRDS("pbmcnew.rds")

expr_matrix <- as(as.sparse(pbmc@assays$RNA@counts), 'sparseMatrix')

monocle.matrix=as.sparse(pbmc@assays$RNA@data)

p_data <- pbmc@meta.data 

p_data$celltype <- pbmc@active.ident  

f_data <- data.frame(gene_short_name = row.names(pbmc),row.names = row.names(pbmc))

pd <- new('AnnotatedDataFrame', data = p_data) 

fd <- new('AnnotatedDataFrame', data = f_data)

cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)

cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)

print(head(fData(cds)))

expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10))

expressed_genes <- VariableFeatures(pbmc)

cds <- setOrderingFilter(cds, expressed_genes)

plot_ordering_genes(cds)

diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~celltype",cores=3) 

deg <- subset(diff, qval < 0.01) 

deg <- deg[order(deg$qval,decreasing=F),]

write.table(deg,file="train.monocle.DEG.xls",col.names=T,row.names=F,sep="\t",quote=F)

ordergene <- rownames(deg) 

saveRDS(ordergene, "ordergene.rds")

cds <- setOrderingFilter(cds, ordergene)  

plot_ordering_genes(cds)

ordergene <- row.names(deg)[order(deg$qval)][1:30]

cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

cds <- orderCells(cds)

cds <- orderCells(cds, root_state = 2)

cds$Group<-ifelse(substr(cds$orig.ident,10,12)=="AIS","AIS Cancer","IAC Cancer")

p1<-plot_cell_trajectory(cds,color_by="Pseudotime", cell_size=2,show_backbone=TRUE) + theme(text = element_text(size = 20));p1

p2<-plot_cell_trajectory(cds,color_by="celltype", cell_size=2,show_backbone=TRUE) + theme(text = element_text(size = 20));p2

saveRDS(cds,"cds.rds")

df <- pData(cds) 

saveRDS(df,"dfstatedata.rds")

p3 <- ggplot(df, aes(Pseudotime, colour = celltype, fill=celltype)) + geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2();p3

keygenes <- head(ordergene,30)

cds_subset <- cds[keygenes,]

Time_diff <- differentialGeneTest(cds[ordergene,], cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")

Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)] 

write.csv(Time_diff, "Time_diff_all.csv", row.names = F)

Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()

clusters <- cutree(p$tree_row, k = 4)

clustering <- data.frame(clusters)

clustering[,1] <- as.character(clustering[,1])

colnames(clustering) <- "Gene_Clusters"

table(clustering)

write.csv(clustering, "Time_clustering_all.csv", row.names = F)

Time_genes <- top_n(Time_diff, n = 30, sort(qval,decreasing = T)) %>% pull(gene_short_name) %>% as.character()

p4 <- plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T,return_heatmap=T,add_annotation_col = NULL);p4

BEAM_res <- BEAM(cds[Time_genes,], branch_point = 1, progenitor_method = "duplicate")

BEAM_res <- BEAM_res[order(BEAM_res$qval),] 

write.csv(BEAM_res, "BEAM_res.csv", row.names = F)

BEAM_res <- read.csv("BEAM_res.csv")

BEAM_genes <- BEAM_res$gene_short_name

p5 <-plot_genes_branched_heatmap(cds[BEAM_genes,],
                               #branch_labels = c("Cell fate 1", "Cell fate 2","3","4"),
                               branch_point = 1, #绘制的是哪个分支
                               num_clusters = 3, #分成几个cluster，根据需要调整
                               show_rownames = T,#热图上现实基因名
                               return_heatmap = T);p5
saveRDS(p5,"p5.rds")
