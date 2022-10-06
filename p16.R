rm(list = ls())
library(Seurat)
library(dplyr)
library(ggplot2)
pbmc <- readRDS('pbmc_tutorial.rds')
pbmc.markers <- FindAllMarkers(pbmc,only.pos = TRUE,min.pct = 0.1,logfc.threshold = 0.25)
top5 <- pbmc.markers%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono","NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc,new.cluster.ids)

p16 <- DoHeatmap(pbmc,features = top5$gene)
png(
  filename = 'p16.png',
  width = 2000,
  height = 800
)
plot(p16)
dev.off()