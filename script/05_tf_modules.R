library(tidyverse)

tf2target <- clusterProfiler::read.gmt("output/02-ifnb_pbmc.regulons.gmt")
head(tf2target)
tf2target$TF <- sub("\\([0-9]+g\\)", "", tf2target$term)
colnames(tf2target)[1:2] <- c("regulon", "target")
tf2target$id <- paste0(tf2target$TF, "-", tf2target$target)
head(tf2target)

adj <- data.table::fread("output/01-step1_adj.tsv", sep = "\t", header = T)
head(adj)
adj$id <- paste0(adj$TF, "-", adj$target)
adj <- subset(adj, id %in% tf2target$id)
head(adj)

length(setdiff(tf2target$id, adj$id))
length(setdiff(adj$id, tf2target$id))

data.use <- left_join(adj, tf2target, by = c("id", "TF", "target"))


hist(data.use$importance, breaks = 100)
summary(data.use$importance)

data.use <- subset(data.use, importance >= 1) # cut.off by importance


regulators <- data.use %>%
  group_by(TF) %>%
  summarise(num_target = n()) %>%
  arrange(desc(num_target)) %>%
  as.data.frame()
rownames(regulators) <- regulators$TF
head(regulators)

hub.genes <- regulators$TF ## TFs
edge.list <- data.use %>% subset(target %in% hub.genes & TF %in% hub.genes)
edge.igraph <- igraph::make_undirected_graph(
  edges = mapply(c, edge.list$TF, edge.list$target, SIMPLIFY = F) %>% Reduce(f=c)
)
edge.igraph


set.seed(1024)
fr.layout <- igraph::layout_with_fr(edge.igraph)
head(fr.layout)
nodes.in.grpah <- igraph::get.vertex.attribute(edge.igraph)[[1]]
nodes.in.grpah
regulators[nodes.in.grpah, "fr1"] = fr.layout[,1]
regulators[nodes.in.grpah, "fr2"] = fr.layout[,2]
head(regulators)

edge.list$x_start = regulators[edge.list$TF, "fr1"]
edge.list$y_start = regulators[edge.list$TF, "fr2"]
edge.list$x_end = regulators[edge.list$target, "fr1"]
edge.list$y_end = regulators[edge.list$target, "fr2"]

set.seed(1024)
clusters <- igraph::cluster_louvain(edge.igraph, resolution = 1)
regulators[clusters$names, "Louvain_cluster"] <- LETTERS[clusters$membership]
table(regulators$Louvain_cluster)

## install python packages
# reticulate::py_install("igraph", pip = TRUE)
# reticulate::py_install("leidenalg", pip = TRUE)
# clusters <- leiden::leiden(edge.igraph, resolution_parameter = 1, seed = 1024)
# regulators[nodes.in.grpah, "Leiden_cluster"] <- LETTERS[clusters]
# table(regulators$Leiden_cluster)


regulators$cluster <- regulators$Louvain_cluster 
regulators


source("R/network_plot.R")
## TF调控模块柱状图，同一个色块代表一个cluster
RegModuleBar(regulators = regulators, topN = 5)

## Combine the variance decomposition results
source("R/vd_plot.R")
vd.res <- readRDS("output/04-1.VD_res.rds")
vd.res
vd.res$TF <- sub("\\(\\+\\)", "", vd.res$gene)
head(vd.res)
vd.res$cluster <- regulators[vd.res$TF, ]$cluster
head(vd.res)
rownames(vd.res) <- vd.res$TF

plot.list <- lapply(LETTERS[1:8], function(clu) {
  VDPlot(vd.res, y = "group", group.by = "cluster", group.show = clu)
})
cowplot::plot_grid(plotlist = plot.list, ncol = 4)


library(Seurat)
seu <- qs::qread("output/03-2.ifnb_pbmc.seurat.aucell.qs")
VlnPlot(seu, group.by = "celltype", features = "Irf7(+)", split.by = "group",
        pt.size = 0, split.plot = TRUE, sort = F)


regulators$var.group <- vd.res[rownames(regulators), ]$group
RegNetVis(regulators = regulators, edge.list = edge.list, size.by = "var.group", topN = 5)
RegNetVis(regulators = regulators, edge.list = edge.list, size.by = "num_target", topN = 5)
RegModuleBar(regulators = regulators, topN = 5)


# Example
nodes.show <- subset(regulators, cluster == "D")$TF
RegNetVis(regulators = regulators, edge.list = edge.list, nodes.show = nodes.show, topN = 10)

nodes.show <- subset(vd.res, group > 0.05)$TF # 响应IFNB的regulon
RegNetVis(regulators = regulators, edge.list = edge.list, nodes.show = nodes.show, topN = 10)

saveRDS(regulators, "output/05-1.regulon_clusters_network.rds")

regulators <- readRDS('output/05-1.regulon_clusters_network.rds')

#### Method 2: TF activity similarity based clustering 
library(Seurat)
source("R/tf_module_utils.R")
seu <- qs::qread("output/03-2.ifnb_pbmc.seurat.aucell.qs")
### 获取380个TF的打分
rasMat <- seu[["AUCell"]]@data
dim(rasMat)
rasMat <- t(rasMat)
pccMat <- cor(rasMat) 

csiMat <- pbapply::pblapply(rownames(pccMat), function(i) sapply(colnames(pccMat), function(j) CSI(i, j)))
csiMat <- do.call(rbind, csiMat)
rownames(csiMat) <- rownames(pccMat)
dim(csiMat)
csiMat[1:5,1:5]

library(dendextend)
library(ggsci)
h = 8
row_dend = as.dendrogram(hclust(dist(pccMat), method = "complete"))
clusters <- dendextend::cutree(row_dend, h = h) # dendextend::cutree()
row_dend = color_branches(row_dend, h = h, col = pal_d3("category20")(20))
plot(row_dend)

library(ComplexHeatmap)
library(circlize)

col_range = c(0.7, 1)
col_fun <- colorRamp2(col_range, c("#FCF8DE", "#253177"))

ht <- Heatmap(
  matrix = csiMat,
  col = col_fun,
  name = "ht1",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_heatmap_legend = FALSE
)

lgd <- Legend(
  col_fun = col_fun,
  title = "",
  at = col_range,
  labels = c("low", "high"),
  direction = "horizontal",
  legend_width = unit(1, "in"),
  border = FALSE
)

{
  draw(ht, heatmap_legend_list = list(lgd), heatmap_legend_side = c("bottom"))
  decorate_heatmap_body("ht1", {
    tree = column_dend(ht)
    ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
    first_index = function(l) which(l)[1]
    last_index = function(l) { x = which(l); x[length(x)] }
    clusters <- names(table(ind))
    x1 = sapply(clusters, function(x) first_index(ind == x)) - 1
    x2 = sapply(clusters, function(x) last_index(ind == x))
    x1 = x1/length(ind)
    x2 = x2/length(ind)
    grid.rect(x = x1, width = (x2 - x1), y = 1-x1, height = (x1 - x2),
              hjust = 0, vjust = 0, default.units = "npc",
              gp = gpar(fill=NA, col="#FCB800", lwd=3))
    grid.text(label = paste0("M",clusters),
              x = x2-length(clusters)/length(ind), y = 1-x1-(x2-x1)/2,
              default.units = "npc",
              hjust = 1, vjust = 0.5,
              gp = gpar(fontsize=12, fontface="bold"))
  })
  decorate_column_dend("ht1", {
    tree = column_dend(ht)
    ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
    first_index = function(l) which(l)[1]
    last_index = function(l) { x = which(l); x[length(x)] }
    clusters <- names(table(ind))
    x1 = sapply(clusters, function(x) first_index(ind == x)) - 1
    x2 = sapply(clusters, function(x) last_index(ind == x))
    grid.rect(x = x1/length(ind), width = (x2 - x1)/length(ind), just = "left",
              default.units = "npc", gp = gpar(fill = pal_d3("category20")(20), alpha=.5, col = NA))
  })
}

tree <- column_dend(ht)
ind <- cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
clusters <- names(table(ind))
regulon.clusters <- data.frame(regulon=names(ind), cluster=paste0("M",ind))
table(regulon.clusters$cluster)
regulon.clusters

k = length(clusters)
cell.info <- seu@meta.data
moduleRasMat <- lapply(paste0("M",1:k), function(x){
  regulon.use <- subset(regulon.clusters, cluster == x)$regulon
  rowMeans(rasMat[, regulon.use, drop=FALSE])
})
names(moduleRasMat) <- paste0("M",1:k)
moduleRasMat <- do.call(cbind, moduleRasMat)
cell.info <- cbind(cell.info, moduleRasMat[rownames(cell.info), ])
cell.info <- cbind(cell.info, FetchData(seu, vars = paste0("umap_", 1:2)))

p.list <- lapply(paste0("M",1:k), function(module){
  data.use <- cell.info
  expression.color <- c("darkblue", "lightblue", "green", "yellow", "red")
  max.val <- quantile(data.use[, module], 0.99)
  low.val <- quantile(data.use[, module], 0.1)
  data.use[, module] <- ifelse(data.use[, module] > max.val, max.val, data.use[, module])
  ggplot(data.use, aes(umap_1, umap_2, color=get(module))) +
    geom_point(size=0.05) +
    theme_bw(base_size = 15) +
    ggtitle(module) +
    facet_wrap(~group) +
    scale_color_gradientn(name = NULL, colors = expression.color) +
    theme(legend.position = "right",
          legend.title = element_blank(),
          plot.title = element_text(hjust = .5, face = "bold", size = 20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black")
    )
})

cowplot::plot_grid(plotlist = p.list, ncol = 3)

## VD plot
vd.res <- readRDS("output/04-1.VD_res.rds")
vd.res[regulon.clusters$regulon, "cluster"] <- regulon.clusters$cluster
plot.list <- lapply(paste0("M", 1:9), function(clu) {
  VDPlot(vd.res, y = "group", group.by = "cluster", group.show = clu)
})
cowplot::plot_grid(plotlist = plot.list, ncol = 3)

saveRDS(regulon.clusters, "output/05-2.regulon_clusters_activity_similarity.rds")
