library(Seurat)
library(tidyverse)
library(patchwork)
setwd(here::here())
source("R/compute_module_score.R")

seu <- readRDS("data/ifnb_pbmc.seurat.rds")
seu <- subset(seu, celltype %in% c("Melanocyte", "Spinous cell",'T cell',
                                   'Hair follicle II','Hair follicle III','Langerhans cell',
                                   'Skeletal muscle','Basal cell','Hair follicle I',
                                   'Endothelial cell'), invert = TRUE)
DimPlot(seu, group.by = "celltype", reduction = "umap")+ggsci::scale_color_d3("category20")

celltype.levels <- c("Fibroblast I", "Fibroblast II", "Fibroblast III", 'Myofibroblast',
                     "Macrophage I", "Macrophage II")
seu$celltype <- factor(seu$celltype, levels = celltype.levels)


regulons <- clusterProfiler::read.gmt("output/02-ifnb_pbmc.regulons.gmt")
head(regulons)

rg.names <- unique(regulons$term)
regulon.list <- lapply(rg.names, function(rg) {
  subset(regulons, term == rg)$gene
})

names(regulon.list) <- sub("[0-9]+g", "\\+", rg.names)
summary(sapply(regulon.list, length))
print(regulon.list[1])

saveRDS(regulon.list, "output/03-1.ifnb_pbmc.regulons.rds")

## AUCell
## RAS = regulon activity score
seu <- ComputeModuleScore(seu, gene.sets = regulon.list, min.size = 10, cores = 10)
seu
DefaultAssay(seu) <- "AUCell"
rownames(seu)

p1 <- FeaturePlot(seu, features = "Creb3l1(+)", split.by = "group")
p2 <- FeaturePlot(seu, features = "Creb3l1", split.by = "group")
(p1 / p2) & scale_color_viridis_c()

VlnPlot(seu, group.by = "celltype", features = "Creb3l1(+)", pt.size = 0,
        split.by = "group", split.plot = TRUE, cols = c("blue", "red")) + ylab("TF activity")
VlnPlot(seu, group.by = "celltype", features = "Creb3l1", pt.size = 0,
        split.by = "group", split.plot = TRUE, cols = c("blue", "red"))
        
seu <- RunUMAP(object = seu,
               features = rownames(seu),
               metric = "correlation", # 注意这里用correlation效果最好
               reduction.name = "umapRAS",
               reduction.key = "umapRAS_")

p1 <- DimPlot(seu, reduction = "umap", group.by = "celltype") + ggsci::scale_color_d3("category20") + NoLegend()
p2 <- DimPlot(seu, reduction = "umap", group.by = "group") + NoLegend()

p3 <- DimPlot(seu, reduction = "umapRAS", group.by = "celltype") + ggsci::scale_color_d3("category20")
p4 <- DimPlot(seu, reduction = "umapRAS", group.by = "group")

(p1 + p3) / (p2 + p4)

DefaultAssay(seu) <- "AUCell"
seu <- ScaleData(seu)
seu <- RunPCA(object = seu,
              features = rownames(seu),
              reduction.name = "pcaRAS",
              reduction.key = "pcaRAS_")

p5 <- DimPlot(seu, reduction = "pcaRAS", group.by = "celltype") + ggsci::scale_color_d3("category20")
p6 <- DimPlot(seu, reduction = "pcaRAS", group.by = "group")
p5 + p6


## PC1 encoding the regulons related to cell type
## PC2 encoding the regulons affected by sample group
## The injury induced transcriptome shift is orthogonal to the cell identity transcriptional programs.

VlnPlot(seu, group.by = "celltype", features = "pcaRAS_1", pt.size = 0,
        split.by = "group", split.plot = TRUE, cols = c("blue", "red"))

VlnPlot(seu, group.by = "celltype", features = "pcaRAS_2", pt.size = 0,
        split.by = "group", split.plot = TRUE, cols = c("blue", "red"))

qs::qsave(seu, "output/03-2.ifnb_pbmc.seurat.aucell.qs")

