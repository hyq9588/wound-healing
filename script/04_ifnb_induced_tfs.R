library(Seurat)
library(tidyverse)
library(patchwork)
library(ggsci)
setwd(here::here())

seu <- qs::qread("output/03-2.ifnb_pbmc.seurat.aucell.qs")

p1 <- DimPlot(seu, reduction = "pcaRAS", group.by = "celltype") + ggsci::scale_color_d3("category20")
p2 <- DimPlot(seu, reduction = "pcaRAS", group.by = "group")
p1 + p2

VlnPlot(seu, group.by = "celltype", features = "pcaRAS_2", split.by = "group",
        pt.size = 0, split.plot = TRUE, sort = F, cols = c("blue", "red"))

VlnPlot(seu, group.by = "celltype", features = "pcaRAS_2", split.by = "group",
        pt.size = 0, split.plot = TRUE, sort = F)+scale_fill_lancet() 

VlnPlot(seu, group.by = "celltype", features = "pcaRAS_1", split.by = "group",
        pt.size = 0, split.plot = TRUE, sort = F)+scale_fill_lancet() 

Loadings(seu, reduction = "pcaRAS")

VizDimLoadings(seu, reduction = "pcaRAS", dims = 2, balanced = T, projected = F)
VlnPlot(seu, group.by = "celltype", features = "Gatad1(+)", split.by = "group",
        pt.size = 0, split.plot = TRUE, sort = F, cols = c("blue", "red")) + ylab("TF activity")


#### Method 2: Variance Decomposition ####
source("R/variance_decompose.R")
DefaultAssay(seu) <- "AUCell"

vd.vars <- c("celltype", "group")
meta.data <- seu@meta.data[, vd.vars]
ras.data <- FetchData(seu, vars = rownames(seu))

vd.res <- VarDecompose(data = ras.data, meta.data = meta.data, vd.vars = vd.vars, cores = 10)
saveRDS(vd.res, "output/04-1.VD_res.rds")

vd.res$PC1.loadings <- Loadings(seu, reduction = "pcaRAS")[vd.res$gene, 1]
vd.res$PC2.loadings <- Loadings(seu, reduction = "pcaRAS")[vd.res$gene, 2]

to.label <- subset(vd.res, celltype > 0.35 & group > 0.05)

ggplot(vd.res, aes(celltype, group, color = abs(PC1.loadings))) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.05, color = "blue", linetype="dashed") +
  geom_vline(xintercept = 0.35, color = "blue", linetype="dashed") +
  ggrepel::geom_text_repel(inherit.aes = F, data = to.label, aes(celltype, group, label=gene)) +
  xlab("Fraction of variance by cell type") +
  ylab("Fraction of variance by group") +
  scale_color_viridis_c() +
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank())
  
to.label2 <- subset(vd.res, celltype > 0.25 & group > 0.15)
ggplot(vd.res, aes(celltype, group, color = abs(PC2.loadings))) + 
  geom_point(size = 2) +
  xlab("Fraction of variance by cell type") +
  ylab("Fraction of variance by group") +
  scale_color_viridis_c() +
  theme_bw(base_size = 15)+
  geom_hline(yintercept = 0.15, color = "blue", linetype="dashed") +
  geom_vline(xintercept = 0.5, color = "blue", linetype="dashed") +
  ggrepel::geom_text_repel(inherit.aes = F, data = to.label2, aes(celltype, group, label=gene))+
  theme(panel.grid = element_blank())

ggplot(vd.res, aes(group)) +
  geom_density()

to.label <- subset(vd.res, group > 0.4)
# to.label <- subset(vd.res, celltype > 0.35 & group > 0.05)
ggplot(vd.res, aes(PC2.loadings, group)) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(inherit.aes = F, data = to.label, aes(PC2.loadings, group, label=gene)) +
  ylab("Fraction of variance by group") +
  theme_bw(base_size = 15)



VlnPlot(seu, group.by = "celltype", features = "Irf7(+)", split.by = "group",
        pt.size = 0, split.plot = TRUE, sort = F) + ylab("TF activity")+ggsci::scale_color_d3("category20")
FeaturePlot(seu, reduction = "umap", features = "Irf7(+)", split.by = "group")
