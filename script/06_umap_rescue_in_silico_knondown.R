library(Seurat)

seu <- qs::qread("output/03-2.ifnb_pbmc.seurat.aucell.qs")

#### 1. remove module D ####

regulon.cluster <- readRDS("output/05-1.regulon_clusters_network.rds")
head(regulon.cluster)

regulons.F <- subset(regulon.cluster, cluster == "D")$TF

regulons.F <- paste0(regulons.F, "(+)")

features.use <- setdiff(rownames(seu), regulons.F)


seu <- RunUMAP(object = seu,
               features = features.use,
               metric = "correlation", # 注意这里用correlation效果最好
               reduction.name = "umapRASrmF",
               reduction.key = "umapRASrmF_")


p1 <- DimPlot(seu, reduction = "umap", group.by = "celltype") + ggsci::scale_color_d3("category20") + NoLegend()
p2 <- DimPlot(seu, reduction = "umap", group.by = "group") + NoLegend()


p3 <- DimPlot(seu, reduction = "umapRAS", group.by = "celltype") + ggsci::scale_color_d3("category20") + NoLegend()
p4 <- DimPlot(seu, reduction = "umapRAS", group.by = "group") + NoLegend()

p5 <- DimPlot(seu, reduction = "umapRASrmF", group.by = "celltype") + ggsci::scale_color_d3("category20")
p6 <- DimPlot(seu, reduction = "umapRASrmF", group.by = "group")

(p1 | p3 | p5) / (p2 | p4 | p6)


regulon.cluster <- readRDS("output/05-2.regulon_clusters_activity_similarity.rds")

regulons.rm <- subset(regulon.cluster, cluster %in% paste0("M",c(1,6,7,9)))$regulon
features.use <- setdiff(rownames(seu), regulons.rm)
seu <- RunUMAP(object = seu,
               features = features.use,
               metric = "correlation", # 注意这里用correlation效果最好
               reduction.name = "umapRASrm",
               reduction.key = "umapRASrm_")

p5 <- DimPlot(seu, reduction = "umapRASrm", group.by = "celltype") + ggsci::scale_color_d3("category20")
p6 <- DimPlot(seu, reduction = "umapRASrm", group.by = "group")

(p1 | p3 | p5) / (p2 | p4 | p6)



#### 3. remove group specific regulons ####
vd.res <- readRDS("output/04-1.VD_res.rds")
regulons.rm <- subset(vd.res, group > 0.1)$gene
features.use <- setdiff(rownames(seu), regulons.rm)

seu <- RunUMAP(object = seu,
               features = features.use,
               metric = "correlation", # 注意这里用correlation效果最好
               reduction.name = "umapRASrm",
               reduction.key = "umapRASrm_")

## UMAP on RAS (remove group > 0.1)
p5 <- DimPlot(seu, reduction = "umapRASrm", group.by = "celltype") + ggsci::scale_color_d3("category20")
p6 <- DimPlot(seu, reduction = "umapRASrm", group.by = "group")

(p1 | p3 | p5) / (p2 | p4 | p6)

#(p1 | p5) / (p2 |  p6)

#### 4. remove celltype specific regulons ####
regulons.rm <- subset(vd.res, celltype > 0.1)$gene
features.use <- setdiff(rownames(seu), regulons.rm)


seu <- RunUMAP(object = seu,
               features = features.use,
               metric = "correlation", 
               reduction.name = "umapRASrmct",
               reduction.key = "umapRASrmct_")

## UMAP on RAS (remove celltype > 0.1)
p7 <- DimPlot(seu, reduction = "umapRASrmct", group.by = "celltype") + ggsci::scale_color_d3("category20")
p8 <- DimPlot(seu, reduction = "umapRASrmct", group.by = "group")

(p1 | p3 | p5 | p7) / (p2 | p4 | p6 | p8)
(p1 | p3 | p7) / (p2 | p4 | p8)
