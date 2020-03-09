# Neuts
Isolation of neutrophils and DE analysis 
groups <- sample(c("healthy", "adjacent", "contralateral", "tumor"), size = 25633, replace = TRUE)
names(groups) <- colnames(neutrophils.seu)
neutrophils.seu <- AddMetaData(object = neutrophils.seu, metadata = groups, col.name = "group")
obj.list <- SplitObject(neutrophils.seu, split.by = "group")

reference.list <- obj.list[c("healthy", "adjacent", "contralateral", "tumor")]
neutrophil.anchors <-FindIntegrationAnchors(object.list = reference.list, dims = 1:50)

library(ggplot2)
library(cowplot)
# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(neutrophil.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
neutrophil.integrated <- ScaleData(neutrophil.integrated, verbose = FALSE)
neutrophil.integrated <- RunPCA(neutrophil.integrated, npcs = 30, verbose = FALSE)
neutrophil.integrated <- FindNeighbors(neutrophil.integrated, dims = 1:10)
neutrophil.integrated <- FindClusters(neutrophil.integrated, resolution = 0.5)
neutrophil.integrated <- RunTSNE(neutrophil.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(neutrophil.integrated, reduction = "tsne", group.by = "group")
plot_grid(p1)



#FindMarkers for Differential Expression Analysis

tumor.markers = FindMarkers(neutrophil.integrated, ident.1 = 'tumor', assay = 'RNA', group.by = 'group')
