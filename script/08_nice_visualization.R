################### 08 Network visualization #######################
## Aims:
## => 1. Fantasy plots showing the TF-target relationships.
library(Seurat)
library(tidyverse)
library(ggraph)
library(tidygraph)
source("R/IO.R")

data <- LoadpySCENICOutput(regulon.gmt = "output/02-ifnb_pbmc.regulons.gmt",
                           adj.mat.file = "output/01-step1_adj.tsv")
head(data)
summary(data$importance)
data <- subset(data, importance > 1)

### 展示调控IFN-gamma通路的转录调控网络
hallmarks <- clusterProfiler::read.gmt("resource/mh.all.v2023.2.Mm.symbols.gmt")
genes <- subset(hallmarks, term == "HALLMARK_INFLAMMATORY_RESPONSE")$gene
source("R/network_plot.R")

# 红色标注代表通路中的基因，大的橘黄色的点代表转录因子
# 灰色的代表TF对应的Target，灰色的外围的圈代表它是该通路的基因
RegulonGraphVis(data, tf.show = c("Spi1", "Irf8", "Irf5", "Irf7", "Rel", "Stat5a", "Ikzf1"),
                targets.show = genes)

RegulonGraphVis(data, tf.show = c("Spi1", "Irf8", "Irf5", "Irf7", "Rel", "Stat5a", "Ikzf1"),
                targets.show = genes, layout = "circle")

RegulonGraphVis(data, tf.show = c("Spi1", "Irf8", "Irf5", "Irf7", "Rel", "Stat5a", "Ikzf1"),
                targets.show = genes, prop = 0.01)

RegulonGraphVis(data, tf.show = c("STAT2", "STAT1", "IRF7", "IRF2", "IRF8", "ETV7", "ELF1"),
                targets.show = genes, prop = NULL, n = 20)
