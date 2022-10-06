rm(list=ls())
library(Seurat)
library(pheatmap)
pbmc <- readRDS('pbmc_tutorial.rds')
table(pbmc$seurat_clusters)

av <- AverageExpression(pbmc,
                        group.by = 'seurat_clusters',
                        assays = 'RNA')
av = av[[1]]
head(av)
cg <- names(tail(sort(apply(av,1,sd)),1000))
View(av[cg,])
View(cor(av[cg,],method = 'spearman'))

data <- cor(av[cg,],method = 'spearman')
colnames(data) <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono","NK", "DC", "Platelet")
rownames(data) <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono","NK", "DC", "Platelet")

windows()
p18 <- pheatmap(data)



