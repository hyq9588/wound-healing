library(tidyverse)
library(Seurat)
library(clusterProfiler)
library(enrichplot)
library(GseaVis)
library(ggrepel)
library(ggplot2)

seu <- qs::qread("output/03-2.ifnb_pbmc.seurat.aucell.qs")
DefaultAssay(seu) <- "RNA"

seu$celltype.stim <- paste(seu$celltype, seu$group, sep = "_")
Idents(seu) <- "celltype.stim"

table(seu$celltype.stim)
interferon.response <- FindMarkers(seu,
                                   ident.1 = "Macrophage I_Wound",
                                   ident.2 = "Macrophage I_Unwound",
                                   logfc.threshold = 0)
head(interferon.response, n = 15)

gene_df <- interferon.response
geneList <- gene_df$avg_log2FC
names(geneList) =  rownames(gene_df)
geneList = sort(geneList, decreasing = TRUE)

head(geneList,10)

hallmarks <- read.gmt("resource/mh.all.v2023.2.Mm.symbols.gmt")
y <- GSEA(geneList, TERM2GENE = hallmarks)
yd <- as.data.frame(y)
dotplot(y,showCategory=30,split=".sign") + facet_grid(~.sign)


geneSetID = c('HALLMARK_INFLAMMATORY_RESPONSE',
              'HALLMARK_INTERFERON_GAMMA_RESPONSE',
              'HALLMARK_TNFA_SIGNALING_VIA_NFKB')
gseaNb(object = y, geneSetID = geneSetID,addPval = T)


mygene <- c('Irf8','Irf7')
gseaNb(object = y,
       geneSetID = 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',addPval = T)

gseaNb(object = y,
       geneSetID = 'HALLMARK_TNFA_SIGNALING_VIA_NFKB',addPval = T)

gseaNb(object = y,
       geneSetID = 'HALLMARK_INFLAMMATORY_RESPONSE',addPval = T)

gseaNb(object = y,
       geneSetID = 'HALLMARK_INTERFERON_GAMMA_RESPONSE',addPval = T)

gseaNb(object = y,
       geneSetID = 'HALLMARK_TNFA_SIGNALING_VIA_NFKB',addPval = T)

gseaNb(object = y,
       geneSetID = 'HALLMARK_COMPLEMENT',addPval = T)


## activated
gseaplot2(y, "HALLMARK_INTERFERON_GAMMA_RESPONSE",color = "red", pvalue_table = T)

## repressed
gseaplot2(y, "HALLMARK_INTERFERON_ALPHA_RESPONSE",color = "red", pvalue_table = T)


regulon.list <- readRDS("output/03-1.ifnb_pbmc.regulons.rds")
regulon.list

sapply(regulon.list, length) %>% summary()

regulons <- read.gmt("output/02-ifnb_pbmc.regulons.gmt")
head(regulons)

## query 1: HALLMARK_INTERFERON_GAMMA_RESPONSE (activated)
genes <- subset(hallmarks, term == "HALLMARK_INTERFERON_GAMMA_RESPONSE")$gene

e.regulon <- enricher(gene = genes, TERM2GENE = regulons, minGSSize = 10, maxGSSize = 3500)
dotplot(e.regulon)+ ggtitle("Interferon Gamma Response")+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid = element_blank())

DefaultAssay(seu) <- "AUCell"
VlnPlot(seu, group.by = "celltype", features = c("Irf7(+)"), split.by = "group",
        split.plot = TRUE, pt.size = 0, cols = c("blue", "red")) + ylab("TF activity")
VlnPlot(seu, group.by = "celltype", features = c("Spi1(+)"), split.by = "group",
        split.plot = TRUE, pt.size = 0, cols = c("blue", "red")) + ylab("TF activity")
VlnPlot(seu, group.by = "celltype", features = c("Irf8(+)"), split.by = "group",
        split.plot = TRUE, pt.size = 0, cols = c("blue", "red")) + ylab("TF activity")


genes <- subset(hallmarks, term == "HALLMARK_INTERFERON_GAMMA_RESPONSE")$gene

e.regulon <- enricher(gene = genes, TERM2GENE = regulons, minGSSize = 10, maxGSSize = 5000)
?dotplot
dotplot(e.regulon)+ ggtitle("HALLMARK_INTERFERON_GAMMA_RESPONSE")+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid = element_blank())

DefaultAssay(seu) <- "AUCell"
VlnPlot(seu, group.by = "celltype", features = c("SPI1(+)"), split.by = "group",
        split.plot = TRUE, pt.size = 0, cols = c("blue", "red")) + ylab("TF activity")
VlnPlot(seu, group.by = "celltype", features = c("ILF2(+)"), split.by = "group",
        split.plot = TRUE, pt.size = 0, cols = c("blue", "red")) + ylab("TF activity")
VlnPlot(seu, group.by = "celltype", features = c("IRF3(+)"), split.by = "group",
        split.plot = TRUE, pt.size = 0, cols = c("blue", "red")) + ylab("TF activity")




table(seu$celltype.stim)
interferon.response <- FindMarkers(seu,
                                   ident.1 = "Fibroblast II_Wound",  
                                   ident.2 = "Fibroblast II_Unwound", 
                                   logfc.threshold = 0)
head(interferon.response, n = 15)


gene_df <- interferon.response

geneList <- gene_df$avg_log2FC

names(geneList) =  rownames(gene_df)
geneList = sort(geneList, decreasing = TRUE)

head(geneList,10)

hallmarks <- read.gmt("resource/mh.all.v2023.2.Mm.symbols.gmt")
y <- GSEA(geneList, TERM2GENE = hallmarks)
yd <- as.data.frame(y)

dotplot(y,showCategory=30,split=".sign") + facet_grid(~.sign)

mygene <- c('Irf8','Irf7')
gseaNb(object = y,
       geneSetID = 'HALLMARK_TNFA_SIGNALING_VIA_NFKB',addPval = T)

gseaNb(object = y,
       geneSetID = 'HALLMARK_INFLAMMATORY_RESPONSE',addPval = T)

gseaNb(object = y,
       geneSetID = 'HALLMARK_INTERFERON_GAMMA_RESPONSE',addPval = T)




## activated
gseaplot2(y, "HALLMARK_INTERFERON_GAMMA_RESPONSE",color = "red", pvalue_table = T)

## repressed
gseaplot2(y, "HALLMARK_INTERFERON_ALPHA_RESPONSE",color = "red", pvalue_table = T)



regulon.list <- readRDS("output/03-1.ifnb_pbmc.regulons.rds")
regulon.list

sapply(regulon.list, length) %>% summary()

regulons <- read.gmt("output/02-ifnb_pbmc.regulons.gmt")
head(regulons)

genes <- subset(hallmarks, term == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")$gene

e.regulon <- enricher(gene = genes, TERM2GENE = regulons, minGSSize = 10, maxGSSize = 3500)
dotplot(e.regulon)+ ggtitle("Epithelial mesenchymal transition")+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid = element_blank())



