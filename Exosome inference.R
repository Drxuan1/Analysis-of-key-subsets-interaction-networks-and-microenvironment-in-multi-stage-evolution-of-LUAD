#devtools::install_github("Danko-Lab/BayesPrism/BayesPrism")
library(BayesPrism)
library(Seurat)
library(reshape)
library(ggplot2)
library(stringr)
library(tidyverse)
library(GeneOverlap)
library(rtracklayer)

setwd("")

bk.datAIS<-readRDS("bk.datAIS.rds")#bk.datAIAC

sc.dat<-readRDS("sc.dat.rds")#sc.dat.IAC

cell.type.labels<-readRDS("cell.type.labels_AIS_sample.rds")#cell.type.labels_IAC_sample

sc.stat <- plot.scRNA.outlier(
  input=sc.dat, 
  cell.type.labels=cell.type.labels,
  species="mm", 
  return.raw=TRUE 

)

head(sc.stat)

sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs", 
                                  gene.group=c("other_Rb","chrM","chrX","chrY","Rb","Mrp","act","hb"),
                                  exp.cells=3)

plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = bk.datAIS
)

sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered,
                                         gene.type = "protein_coding")

cell.state.labels<-cell.type.labels

myPrism <- new.prism(
  reference=sc.dat.filtered.pc, 
  mixture=bk.datAIS,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.state.labels,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

saveRDS(myPrism,"myPrism_AIS_sample.rds")

bp.res <- run.prism(prism = myPrism, n.cores=2)

bp.res

slotNames(bp.res)

save(bp.res, file="bp.res_AIS_sample.rdata")

theta <- get.fraction(bp=bp.res,
                      which.theta="final",
                      state.or.type="type")

head(theta)

write.csv(theta,file="theta_AIS_sample.csv")

theta.cv <- bp.res@posterior.theta_f@theta.cv

head(theta.cv)

ratio <- as.data.frame(theta)

ratio <- t(ratio)

ratio <- as.data.frame(ratio)

ratio <- tibble::rownames_to_column(ratio)

ratio <- melt(ratio)

colourCount = length(ratio$rowname)

colnames(ratio)<-c("Celltype","Counts","Ratio")

ratio<-ratio[order(ratio$Ratio),]

#saveRDS(ratio,"ratio.rds")

ratio$celltype <- sapply(strsplit(ratio$Celltype, "_"), "[", 1)

ratio$sample <- sapply(strsplit(ratio$Celltype, "_"), "[", 2)

saveRDS(ratio,"ratio_AIS_sample.rds")

ratio <- readRDS("ratio_AIS_sample.rds")

cell_proportions_sample <- ratio %>%

  group_by(celltype) %>%

  mutate(
    group_total = sum(Ratio, na.rm = TRUE)
  ) %>%

  mutate(sample_count = n()) %>%

  ungroup()  %>%  
  arrange(group_total)  

total_labels <- cell_proportions_sample %>%
  distinct(celltype, group_total) %>%
  mutate(
    celltype = fct_inorder(celltype),  
    y = as.integer(celltype),         
    x = group_total * 1.05
  )

ais_samples <- unique(cell_proportions_sample$sample)
group_colors <- c("#3B9AB2")

ais_pal <- colorRampPalette(c(group_colors, "#A6D8F0"))(length(ais_samples))

ggplot(cell_proportions_sample, 
       aes(x = celltype, y = Ratio, fill = sample))+
  
  geom_bar(stat = "identity", position = "stack")+
  
  coord_flip()+
  
  scale_x_discrete(limits = unique(cell_proportions_sample$celltype))+

  scale_fill_manual(values = ais_pal, name = "Sample")+

  annotate("text", 
           x = total_labels$y, 
           y = total_labels$x, 
           label = paste0(round(total_labels$group_total, 3), "%"), 
           size = 4, fontface = "bold")+
  labs(x="Celltype")+

  expand_limits(y = max(total_labels$x) * 1.1)+
  
  theme_minimal()

