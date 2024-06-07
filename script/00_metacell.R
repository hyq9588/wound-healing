library(Seurat)
library(tidyverse)
setwd(here::here())
source("R/makeMetaCells.R")
if (!dir.exists("output")) {
  dir.create("output")
}

# ifnb_pbmc.seurat.rds is hamony_seurat_annotaion.rds
seu <- readRDS("data/ifnb_pbmc.seurat.rds")
DimPlot(seu, group.by = "celltype", split.by = "group") & ggsci::scale_color_d3("category20")

table(seu$orig.ident)
seu.list <- SplitObject(seu, split.by = "orig.ident")
seu.list[[1]]@project.name <- "pbmc_wound1"
seu.list[[2]]@project.name <- "pbmc_wound2"
seu.list[[3]]@project.name <- "pbmc_wound3"
seu.list[[4]]@project.name <- "pbmc_unwound1"
seu.list[[5]]@project.name <- "pbmc_unwound2"

metacells.list <- lapply(seq_along(seu.list), function(ii) {
  makeMetaCells(
    seu       = seu.list[[ii]],
    min.cells = 10,
    reduction = "umap",
    dims      = 1:2,
    k.param   = 10,
    cores     = 20)
})

mc.mat <- lapply(metacells.list, function(mc) mc$mat) %>% Reduce(cbind, .)
mc.cellmeta <- lapply(metacells.list, function(mc) mc$metadata) %>% Reduce(rbind, .)
saveRDS(mc.mat, "output/00-1.mc.mat.rds")

motif2tfs <- data.table::fread("cisTarget_db/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl")
TFs <- sort(unique(motif2tfs$gene_name))
writeLines(TFs, "cisTarget_db/hsa_hgnc_tfs.motifs-v10.txt")

mc.mat <- readRDS("output/00-1.mc.mat.rds")
expr.in.cells <- rowSums(mc.mat > 0)
mc.mat <- mc.mat[expr.in.cells >= 5, ]

cisdb <- arrow::read_feather("cisTarget_db/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
genes.use <- intersect(colnames(cisdb), rownames(mc.mat))
length(genes.use)
dim(mc.mat)
mc.mat <- mc.mat[genes.use, ]
dim(mc.mat)

loom <- SCopeLoomR::build_loom(
  file.name         = "output/00-2.mc_mat_for_step1.loom",
  dgem              = mc.mat,
  default.embedding = NULL
)
loom$close()

rm(loom)
gc()

