library(Seurat)
library(ggplot2)
library(infercnv)
library(tidyverse)
library(ggpubr)
library(rjags)
library(reticulate)
library(ComplexHeatmap)
library(RColorBrewer)

setwd("")

rm(list=ls())

options(stringsAsFactors = F)

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix="dataepien.txt",
                                    annotations_file="groupepien.txt",
                                    gene_order_file="geneFile.txt",
                                    ref_group_names=c("Endothelial cells"))  

infercnv_obj2 <- infercnv::run(infercnv_obj,
                              cutoff=0.1, 
                              out_dir=  "infercnv_out" ,  
                              cluster_by_groups=T,   
                              hclust_method="ward.D2",denoise=F, HMM=F,write_expr_matrix = T,write_phylo = T,plot_chr_scale = T)
#
infercnv_obj <- readRDS("run.final.infercnv_obj")

gene_order <- infercnv_obj@gene_order

gene_order$gene <- rownames(gene_order) 

expr <- infercnv_obj@expr.data 

normal_loc <- infercnv_obj@reference_grouped_cell_indices 

test_loc <- infercnv_obj@observation_grouped_cell_indices 

anno.df <- data.frame(  
              idents=c(colnames(expr)[normal_loc$`Endothelial cells`], 
               colnames(expr)[test_loc$AT1],
               colnames(expr)[test_loc$AT2],
               colnames(expr)[test_loc$CILIATED],
               colnames(expr)[test_loc$CLUB]),
            Celltype=c(rep("Endothelial cells",length(normal_loc$`Endothelial cells`)),          
                    rep("AT1",length(test_loc$AT1)),          
                    rep("AT2",length(test_loc$AT2)),          
                    rep("CILIATED",length(test_loc$CILIATED)),
                    rep("CLUB",length(test_loc$CLUB))))

head(anno.df)

rownames(anno.df) <- anno.df[,1]

anno.df <- anno.df[,-1,drop = F]

genenames <- rownames(expr)

geneInfor<-read.table("geneFile.txt",row.names = 1)

sub_geneInfor <-  geneInfor[intersect(genenames,rownames(geneInfor)),]

expr <- expr[intersect(genenames,rownames(geneInfor)),]

head(sub_geneInfor,4)

expr[1:4,1:4]

anno.df$Class <- ifelse(anno.df$Celltype %in% c("Endothelial cells"),"Reference","Observe")

top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = paste0("chr",1:22),labels_gp = gpar(cex = 1.5),labels_rot = 90))

saveRDS(top_anno,"top_anno.rds")

left_anno <- rowAnnotation(Class = anno.df$Class, 
                          Celltype = anno.df$Celltype, 
                          col = list(Class = c("Observe" = "#328D5CFF","Reference" = "#B326B9FF"),
                                     Celltype=c("Endothelial cells"="#B3DE69","AT1"="#C73E3A","AT2" = "red","CILIATED" = "#A35E47","CLUB" = "#F46D43")))

p <- Heatmap(t(expr)[rownames(anno.df),],
             col = colorRamp2(c(0.85,1,1.15),c("#8DD3C7","white","red4")),            
             cluster_rows = F,
             cluster_columns = F,
             show_column_names = F,
             show_row_names = F,             
             column_split = factor(sub_geneInfor$V2, paste("chr",1:22,sep = "")), 
             column_gap = unit(2, "mm"), 
             heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.85,1,1.15),legend_height = unit(3, "cm")),             
             top_annotation = top_anno,             
             left_annotation = left_anno,
             row_title = NULL,
             column_title = "inferCNV",
             column_title_gp = gpar(fontsize = 40),
             border_gp = gpar(col = "black", lty = 1),
             row_split = c(rep("Group1", max(which(anno.df$Celltype %in% c("Endothelial cells")))), rep("Group2", nrow(anno.df)-max(which(anno.df$Celltype %in% c("Endothelial cells")))))
)

draw(p, heatmap_legend_side = "right")

saveRDS(p,"p.rds")
