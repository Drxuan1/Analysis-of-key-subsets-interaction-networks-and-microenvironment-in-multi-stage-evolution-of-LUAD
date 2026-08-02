library(limma)
library(NMF)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(CellChat)
library(Seurat)

setwd("")        

pbmc<-readRDS("pbmcMACAT2CD4CD8CAF1B_END.rds")

pbmc$labels<- Idents(pbmc)

pbmc$Group=ifelse(grepl("IAC",pbmc@meta.data$orig.ident),"IAC","AIS")

pbmcAIS<-subset(pbmc,Group %in% c("AIS"))

pbmcIAC<-subset(pbmc,Group %in% c("IAC"))

cellchat <- createCellChat(normalizeData(pbmcAIS[['RNA']]@counts), 
                           meta = pbmcAIS[[]],
                           group.by = 'labels')

#cellchat <- createCellChat(normalizeData(pbmcIAC[['RNA']]@counts), 
#                           meta = pbmcIAC[[]],
#                           group.by = 'labels')

groupSize <- as.numeric(table(cellchat@idents))    

CellChatDB <- CellChatDB.human      

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")

cellchat@DB <- CellChatDB.use

pdf(file="COMM01.DatabaseCategory.pdf", width=7, height=5)

showDatabaseCategory(CellChatDB)

dev.off()

cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)  

cellchat <- identifyOverExpressedInteractions(cellchat)    

cellchat <- projectData(cellchat, PPI.human)  

cellchat <- computeCommunProb(cellchat)

cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net=subsetCommunication(cellchat)

write.table(file="COMM02.Comm.network.xls", df.net, sep="\t", row.names=F, quote=F)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

pdf(file="COMM03.cellNetworkCount.pdf", width=6, height=6)

netVisual_circle(cellchat@net$count, vertex.weight = groupSize,sources.use = c("SPP1.MBMS.M2"),targets.use = c("AT2","CD4.Exhaust","CD8.Exhaust","B cells","Endothelial cells","CAF1"),weight.scale = T, label.edge= F, title.name = "Number of interactions")

dev.off()

pdf(file="COMM04.cellNetworkWeight.pdf", width=6, height=6)

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,sources.use = c("SPP1.MBMS.M2"),targets.use = c("AT2","CD4.Exhaust","CD8.Exhaust","B cells","Endothelial cells","CAF1"), weight.scale = T, label.edge= F, title.name = "Interaction strength")

dev.off()

pdf(file=paste0("COMM05.bubble.pdf"), width=5, height=5)

netVisual_bubble(cellchat,sources.use = c("SPP1.MBMS.M2"),targets.use = c("AT2","CD4.Exhaust","CD8.Exhaust","B cells","Endothelial cells","CAF1"),signaling = "SPP1",remove.isolate = FALSE, angle.x = 45)

dev.off()

