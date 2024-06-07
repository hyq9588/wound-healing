library(Seurat)
library(tidyverse)
library(patchwork)
library(harmony)
library(SeuratWrappers)
library(scRNAtoolVis)
library(ggplot2)
library(ggsci)
library(scales) 
library(RColorBrewer)
library(scCustomize)
library(ggalluvial)
library(ggSCvis)
library(grid)

# ifnb_pbmc.seurat.rds is hamony_seurat_annotaion.rds
scobj <- readRDS(file = 'data/ifnb_pbmc.seurat.rds')
celltype.levels <- c("Fibroblast I", "Fibroblast II", "Fibroblast III", 'Myofibroblast',
                     "Macrophage I", "Macrophage II",'Hair follicle I','Hair follicle II',
                     "Hair follicle III",'Basal cell','Spinous cell','Langerhans cell','T cell',
                     'Endothelial cell','Melanocyte','Skeletal muscle')
scobj@meta.data$celltype <- factor(scobj@meta.data$celltype, levels = celltype.levels)

DimPlot(scobj, reduction = "umap",group.by = 'celltype',split.by = 'group')+ggsci::scale_color_d3("category20")+ 
  theme(axis.line = element_blank(),  
        axis.text = element_blank(),axis.ticks = element_blank())+labs(x = "", y = "")+ggtitle('')
       

table(scobj$group)  
prop.table(table(Idents(scobj)))
table(Idents(scobj), scobj$group)
Cellratio <- prop.table(table(Idents(scobj), scobj$group), margin = 2)
Cellratio <- as.data.frame(Cellratio)

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

colnames(Cellratio)[1] <- 'Cell Type'
Cellratio$`Cell Type` <- factor(Cellratio$`Cell Type`,levels = celltype.levels)

ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = `Cell Type`),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Cell Ratio')+
  scale_fill_manual(values = allcolour)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))

ggplot(Cellratio, aes(x="", y=Freq, fill=`Cell Type`)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  facet_wrap(~ Var2, nrow=1) +
  theme_void() +
  theme(legend.position = "bottom")+
  scale_fill_manual(values = allcolour)
  

top_markers<- genes <- c(
  "Col1a2", "Col1a1", "Dcn", "Ptx3", "Gpx3", 
  "Col1a2", "Crabp1", "Mfap4", "Col1a2","Itgam", "Acta2", 
  "Rgs5", "Tagln", "Lyz2", "Pf4", "Cd14", 
  "Cd68", "Cd74", "H2-Ab1", "Itgax", "Cd86", 
  "Ifitm1", "Krt17", "Krt79", "Fst", "Sostdc1", 
  "Defb6", "Krt79", "Krt17", "Cst6", "Krt17", 
  "Krt79", "Ube2c", "Krt14", "Krt15", "Krt5", "Krt10", 
  "Krt1", "Krtdap", "Cd207", "H2-M2", "Cd74", "H2-Ab1", 
  "H2-Eb1", "Cd3g", "Nkg7", "Gzmc", "Aqp1", 
  "Plvap", "Cldn5", "Fabp4", "Dct", "Plp1", 
  "Acta1", "Mylpf", "Tnnc2", "Myl1"
)


scobj <- ScaleData(scobj, features = rownames(scobj))
AverageHeatmap(object = scobj,markerGene = top_markers,fontsize = 8,cluster.order=celltype.levels)


library(RColorBrewer)
display.brewer.all()
mycol <- brewer.pal(13, "Set3")
AverageHeatmap(object = scobj,
               markerGene = top_markers,
               annoCol = TRUE,
               myanCol = mycol,fontsize = 9,
               row_title = '') 
