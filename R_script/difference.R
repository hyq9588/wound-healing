library(tidyverse)
library(Seurat)
library(clusterProfiler)
library(enrichplot)
library(GseaVis)
library(ggrepel)
library(ggplot2)
library(cowplot)  
library(ggpubr)  
library(gridExtra)  

seu <- qs::qread("output/03-2.ifnb_pbmc.seurat.aucell.qs")
DefaultAssay(seu) <- "RNA"

seu$celltype.stim <- paste(seu$celltype, seu$group, sep = "_")
Idents(seu) <- "celltype.stim"

table(seu$celltype.stim)

DEGs_PMN0 <- FindMarkers(seu,  ident.1 = "Macrophage I_Wound",
                               ident.2 = "Macrophage I_Unwound",
                               logfc.threshold = 0)

DEGs_PMN0$difference <- DEGs_PMN0$pct.1 - DEGs_PMN0$pct.2
DEGs_PMN0_sig <- DEGs_PMN0[which(DEGs_PMN0$p_val_adj<0.05 & abs(DEGs_PMN0$avg_log2FC) >0.25),]
DEGs_PMN0_sig$label <- rownames(DEGs_PMN0_sig)

ggplot(DEGs_PMN0, aes(x = difference, y = avg_log2FC)) + 
  geom_point(size = 2, color = "grey60") +
  geom_text(data = DEGs_PMN0[which(rownames(DEGs_PMN0) == "Thbs1"), ],
            aes(label = 'Thbs1'), hjust = 1.5, vjust = -1, size = 5, color = "black") +
  geom_point(data = DEGs_PMN0[which(DEGs_PMN0$p_val_adj < 0.05 & DEGs_PMN0$avg_log2FC > 1), ],
             aes(x = difference, y = avg_log2FC), size = 2, color = "red") +
  geom_point(data = DEGs_PMN0[which(DEGs_PMN0$p_val_adj < 0.05 & DEGs_PMN0$avg_log2FC < -1), ],
             aes(x = difference, y = avg_log2FC), size = 2, color = "blue") +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 15),
        axis.line = element_line(color = "black", size = 1)) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1) +
  geom_vline(xintercept = 0, lty = 2, lwd = 1) +
  labs(x = "Delta Percent", y = "Log-fold Change")+ggtitle('Macrophage I')+theme(plot.title = element_text(hjust = 0.5))

# ylab("Log-fold Change") +
# xlab("Delta Percent")

seu <- readRDS("data/ifnb_pbmc.seurat.rds")
mc.mat <-  readRDS('output/00-1.mc.mat.rds')
seu2 <- CreateSeuratObject(mc.mat)
seu2 <- NormalizeData(seu2)

exprSet <- seu2@assays[['RNA']]@data
exprSet <- as.data.frame(t(as.matrix(exprSet)))

library(ggplot2)
library(ggpubr)
library(ggExtra)

xlab="Thbs1"
ylab='Irf7'
x=as.numeric(exprSet[,xlab])
y=as.numeric(exprSet[,ylab])
df1=as.data.frame(cbind(x,y))

library(ggsci)
p1=ggplot(df1, aes(x, y)) + 
  xlab("Thbs1") + ylab('Irf7')+
  geom_point(color='black') + geom_smooth(method="lm",formula = y ~ x,color='blue') + theme_bw()+
  stat_cor(method = 'spearman', aes(x =x, y =y))
ggMarginal(p1, type="density", xparams=list(fill = "#F8766D"), yparams=list(fill = "#00BFC4"))

ggscatter(df1, x = "x", y = "y",
          size = 1.5,
          add = "reg.line",  # 添加回归线
          add.params = list(color = "#77C034", fill = "#C5E99B", size = 1),  # 自定义回归线的颜色
          conf.int = TRUE  # 添加置信区间
) +
  stat_cor(method = "spearman", label.x = 2.5, label.y = 2.0, label.sep = "\n") +
  xlab("Thbs1 expression") +
  ylab("Irf7 expression")
