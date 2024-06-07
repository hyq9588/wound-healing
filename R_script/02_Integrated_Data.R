library(Seurat)
library(tidyverse)
library(patchwork)
library(harmony)
library(SeuratWrappers)
library(ggplot2)

path <- "output/data" 
filenames <- list.files(path, pattern = "\\.rds$") 
full_paths <- file.path(path, filenames)
data <- lapply(full_paths, readRDS)

scobj <- merge(x=data[[1]], y = data[-1])

rm(list = ls(pattern="data.*"))

scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

scobj <- subset(scobj, subset = nFeature_RNA > 300 & nFeature_RNA < 7500 & 
                  nCount_RNA > 1000 )

VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(scobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

scobj <- NormalizeData(scobj)
scobj <- FindVariableFeatures(scobj, selection.method = "vst", nfeatures = 2000)

scobj <- ScaleData(scobj, features = rownames(scobj))
scobj <- RunPCA(scobj, features = VariableFeatures(object = scobj),reduction.name = "pca")

ElbowPlot(scobj, ndims=50, reduction="pca") 

scobj <- RunUMAP(scobj,reduction = "pca", dims = 1:15, reduction.name = "umap_naive")

scobj <- RunHarmony(scobj,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")

scobj <- RunUMAP(scobj, reduction = "harmony", dims = 1:15,reduction.name = "umap")
p1 <- DimPlot(scobj, reduction = "umap_naive",group.by = "orig.ident")+ggsci::scale_color_d3("category20")
p2 <- DimPlot(scobj, reduction = "umap",group.by = "orig.ident")+ggsci::scale_color_d3("category20")
p1+p2

scobj <- FindNeighbors(scobj,dims = 1:15)
scobj <- FindClusters(scobj, resolution = seq(0.2,1.2,0.1))

metadata <- scobj@meta.data
library(clustree)
clustree(scobj)

scobj@meta.data$seurat_clusters <- scobj@meta.data$RNA_snn_res.0.3
Idents(scobj) <- "seurat_clusters"
DimPlot(scobj, reduction = "umap")+ggsci::scale_color_d3("category20")
DimPlot(scobj, reduction = "umap", group.by = "orig.ident")
DimPlot(scobj, reduction = "umap", split.by = "orig.ident")

scobj@assays$RNA@scale.data <- matrix()
scobj@reductions$umap_naive <- NULL
saveRDS(scobj,file = "output/data/hamony_seurat_unannotaion.rds")

