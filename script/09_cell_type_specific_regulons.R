library(Seurat)
library(tidyverse)
library(data.table)
library(ggrepel)
# The philentropy package is used to calculate JS divergence,
# and a detailed description of this package is provided here:http://blog.fens.me/r-entropy/
library(philentropy) # install.packages("philentropy")
source("R/regulon_specificity.R")

#### Calculation of RSS Matrix (Regulon Specificity Score)

## The RSS Example
a <- c(0.7, 0.5, 0.1, 0.1, 0) # TF activity score
## cell type: a,a,b,b,c
b1 <- c(1, 1, 0, 0, 0) # cell type a
b2 <- c(0, 0 ,1, 1, 0) # cell type b
b3 <- c(0, 0, 0, 0, 1) # cell type c

1 - philentropy::JSD(rbind(a, b1), unit = 'log2', est.prob = "empirical")
1 - philentropy::JSD(rbind(a, b2), unit = 'log2', est.prob = "empirical")
1 - philentropy::JSD(rbind(a, b3), unit = 'log2', est.prob = "empirical")
cor(a, b1, method = )
 
seu <- qs::qread("output/03-2.ifnb_pbmc.seurat.aucell.qs")
table(seu$group)
seu <- subset(seu, group == "Unwound")
DefaultAssay(seu) <- "AUCell"

# rows are cells and columns are features (regulons)
rasMat <- t(seu[["AUCell"]]@data)
dim(rasMat)

## 2.2 cell type indicate matrix
ctMat <- calIndMat(seu$celltype)
dim(ctMat)

## 2.3 RSS matrix
rssMat <- calRSSMat(rasMat, ctMat)
dim(rssMat)
write.csv(rssMat,file = 'Figure/Figure4/rssMat.csv',row.names = T)


## Overall
plot.list <- lapply(levels(seu$celltype), function(xx) {
  PlotRegulonRank(rssMat, xx)
})
cowplot::plot_grid(plotlist = plot.list, ncol = 3)

