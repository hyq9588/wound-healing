library(Seurat)
library(tidyverse)
library(patchwork)
library(harmony)
library(SingleCellExperiment)

# data come from 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142471

scdata <- Read10X(data.dir = "data/Wound1/")
scobj <- CreateSeuratObject(counts = scdata, 
                            project = "Wound1", 
                            min.cells = 3, 
                            min.features = 300)
scobj@meta.data$group = "Wound"
saveRDS(scobj,file = "output/data/Wound1_seurat.rds")

scdata <- Read10X(data.dir = "data/Wound2/")
scobj <- CreateSeuratObject(counts = scdata, 
                            project = "Wound2", 
                            min.cells = 3, 
                            min.features = 300)
scobj@meta.data$group = "Wound"
saveRDS(scobj,file = "output/data/Wound2_seurat.rds")

scdata <- Read10X(data.dir = "data/Wound3/")
scobj <- CreateSeuratObject(counts = scdata, 
                            project = "Wound3", 
                            min.cells = 3, 
                            min.features = 300)
scobj@meta.data$group = "Wound"
saveRDS(scobj,file = "output/data/Wound3_seurat.rds")

scdata <- Read10X(data.dir = "data/Unwound1/")
scobj <- CreateSeuratObject(counts = scdata, 
                            project = "Unwound1", 
                            min.cells = 3, 
                            min.features = 300)
scobj@meta.data$group = "Unwound"
saveRDS(scobj,file = "output/data/Unwound1_seurat.rds")

scdata <- Read10X(data.dir = "data/Unwound2/")
scobj <- CreateSeuratObject(counts = scdata, 
                            project = "Unwound2", 
                            min.cells = 3, 
                            min.features = 300)
scobj@meta.data$group = "Unwound"
saveRDS(scobj,file = "output/data/Unwound2_seurat.rds")



