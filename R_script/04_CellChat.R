library(CellChat)
library(tidyverse)
library(ggalluvial)
library(Seurat)
library(data.table)
library(ggsci)


scobj <- readRDS('output/data/hamony_seurat_annotaion.rds')
scobj <- subset(scobj,orig.ident %in% c('Wound1','Wound2','Wound3'))

scobj <- scobj[,Idents(scobj) %in% c('Fibroblast I','Fibroblast II','Fibroblast III','Myofibroblast','Macrophage I','Macrophage II')]

data.input <- GetAssayData(scobj, assay = "RNA", slot = "data")

identity <- subset(scobj@meta.data, select = "celltype")
identity$celltype <- as.character(identity$celltype)

cellchat <- createCellChat(object = data.input, meta = identity,group.by = 'celltype')


CellChatDB <- CellChatDB.mouse

CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor")

cellchat@DB <- CellChatDB.use 

cellchat <- subsetData(cellchat)
future::plan("multicore", workers = 4)

cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- projectData(cellchat, PPI.mouse)

unique(cellchat@idents)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE)


cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

levels(cellchat@idents)             
vertex.receiver = c(2,6)            
cellchat@netP$pathways              
pathways.show <- "THBS"            


par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling ="THBS", layout = "circle")
