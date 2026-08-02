library(nichenetr)
library(Seurat)
library(tidyverse)
library(cowplot)

setwd("")

#ligand_target_matrix <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
#lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))

ligand_target_matrix <- readRDS("ligand_target_matrix.rds")

lr_network <- readRDS("lr_network.rds")

weighted_networks <- readRDS("weighted_networks.rds")

pbmc <- readRDS("pbmcMACAT2CD4CD8CAF1B_END.rds")

pbmc$Group <- ifelse(substr(pbmc$orig.ident,10,12)=="AIS","AIS","IAC")

table(pbmc@active.ident)

nichenet_output <- nichenet_seuratobj_aggregate(seurat_obj = pbmc, 
                                               top_n_ligands = 100,
                                               receiver = c("AT2"), #"CD4.Exhaust"，"CD8.Exhaust"，"CAF1"，"B cells"，"Endothelial cells"
                                               sender = c("SPP1.MBMS.M2"),
                                               condition_colname = "Group", 
                                               condition_oi = "IAC", 
                                               condition_reference = "AIS", 
                                               ligand_target_matrix = ligand_target_matrix, 
                                               lr_network = lr_network,
                                               weighted_networks = weighted_networks 
)

p <- nichenet_output$ligand_target_heatmap+
  xlab("Target(AT2)") + #"CD4.Exhaust"，"CD8.Exhaust"，"CAF1"，"B cells"，"Endothelial cells"
  ylab("Ligand(SPP1.MBMS.M2)")

data <- p$data

data <- data[which(data$y%in%c("SPP1")),]

data <- data[data$score>0,]

data <- data[order(data$score,decreasing = T),]

p$data <- data[1:10,]

p

saveRDS(data, "data_AT2.rds")#"CD4.Exhaust"，"CD8.Exhaust"，"CAF1"，"B cells"，"Endothelial cells"
#
data_AT2 <- readRDS("data_AT2.rds")#"CD4.Exhaust"，"CD8.Exhaust"，"CAF1"，"B cells"，"Endothelial cells"

data_AT2_TOP <- data_AT2[1:10, ]

data_AT2_TOP <- data_AT2_TOP[order(data_AT2_TOP$score),]

names(data_AT2_TOP)[3] <- c("Potential")

data_AT2_TOP$Rank <- as.numeric(rownames(data_AT2_TOP))

p1 <- ggplot(data_AT2_TOP,aes(y,x,size = Rank))+#"CD4.Exhaust"，"CD8.Exhaust"，"CAF1"，"B cells"，"Endothelial cells"
  
  geom_point(shape=21,aes(fill= `Potential`),position =position_dodge(0))+
  
  scale_y_discrete(limits = data_AT2_TOP$x)+ 
  
  theme_minimal()+xlab("SPP1(SPP1.MDMs.M2)")+ylab("Target Gene(AT2)")+
  
  scale_size_continuous(range=c(0,10),labels = function(x) as.integer(x)) + theme_bw()+ 
  
  scale_fill_gradient2(low = "#4682B4", mid = "grey",high = "#E54924",midpoint = median(data_AT2_TOP$`Potential`))+#mid = "grey",
  
  theme(legend.position = "right",legend.box = "vertical",
        legend.margin=margin(t= 0, unit='cm'),
        legend.spacing = unit(0,"in"),
        axis.text.x  = element_blank(),
        axis.text.y  = element_text(color="black",size=24),
        legend.text = element_text(size =24,color="black"),
        legend.title = element_text(size =24,color="black"),
        axis.title.y=element_text(vjust=1,size=24),
        element_blank()
  )+
  
  labs(x = "SPP1(SPP1.MDMs.M2)",y = "Target Gene(AT2)");p1
