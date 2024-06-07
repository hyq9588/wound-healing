library(Seurat)
library(tidyverse)
library(patchwork)
library(harmony)
library(SeuratWrappers)

scobj <- readRDS(file = "output/data/hamony_seurat_unannotaion.rds")
DimPlot(scobj, reduction = "umap", label = T)+ggsci::scale_color_d3("category20")
          
all_markers <- FindAllMarkers(object = scobj)
saveRDS(all_markers,file = "output/data/Seurat_markers.rds")
all_markers <- readRDS('output/data/Seurat_markers.rds')

top_markers <- all_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC))%>%
  dplyr::slice(1:10) %>%
  ungroup() 
top_markers$gene[1:10]

head(Idents(scobj))
Idents(scobj) <- "seurat_clusters"
scobj <- RenameIdents(scobj,
                      "0"= "Fibroblast I",
                      "1"= "Basal cell", 
                      "2"= "Spinous cell", 
                      "3"= "Hair follicle I", 
                      "4"= "Fibroblast II", 
                      "5"= "Macrophage I",
                      "6"= "Macrophage II", 
                      "7"= "T cell", 
                      "8"= "Fibroblast III",
                      "9"= "Myofibroblast", 
                      "10"= "Hair follicle II", 
                      "11"= "Hair follicle III",
                      "12"= "Langerhans cell", 
                      "13"= "Endothelial cell",
                      "14"= "Fibroblast III",
                      "15"= "Melanocyte",
                      "16"= "Skeletal muscle"
)
head(Idents(scobj))

metadata <- scobj@meta.data
scobj@meta.data$celltype = Idents(scobj)
table(scobj$celltype)
subset_df <- metadata[metadata$orig.ident %in% c('Wound1','Wound2','Wound3'),]
subset_df2 <- metadata[!metadata$orig.ident %in% c('Wound1','Wound2','Wound3'),]
table(subset_df$celltype)
table(subset_df2$celltype)

DimPlot(scobj, reduction = "umap")+ggsci::scale_color_d3("category20")
saveRDS(scobj,file = "output/data/hamony_seurat_annotaion.rds")

