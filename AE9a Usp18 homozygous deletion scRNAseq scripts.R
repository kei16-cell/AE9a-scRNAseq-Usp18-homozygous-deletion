library(dplyr)
library(Seurat)
library(cowplot)
library(viridis)

ffoil.data <- Read10X(data.dir = "/Users/karim/Desktop/ff_oil_filtered_feature_bc_matrix")
fftam.data <- Read10X(data.dir = "/Users/karim/Desktop/ff_tam_filtered_feature_bc_matrix")

ffoil <- CreateSeuratObject(counts = ffoil.data, project = "ff_oil", min.cells = 3)
ffoil[["percent.mt"]] <- PercentageFeatureSet(ffoil, pattern = "^mt-")
ffoil[["percent.ribo"]] <- PercentageFeatureSet(ffoil, pattern = "^Rp[sl][[:digit:]]")

fftam <- CreateSeuratObject(counts = fftam.data, project = "ff_tam", min.cells = 3)
fftam[["percent.mt"]] <- PercentageFeatureSet(fftam, pattern = "^mt-")
fftam[["percent.ribo"]] <- PercentageFeatureSet(fftam, pattern = "^Rp[sl][[:digit:]]")

VlnPlot(ffoil, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 3)
VlnPlot(fftam, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 3)

ffoil$tam <- "oil"
ffoil <- subset(ffoil, subset = nFeature_RNA > 200 & percent.mt < 15)
ffoil <- ScaleData(object = ffoil, vars.to.regress = c("percent.ribo","Rn45s"))
ffoil <- NormalizeData(ffoil, verbose = FALSE)
ffoil <- FindVariableFeatures(ffoil, selection.method = "vst", nfeatures = 2000)

fftam$tam <- "tam"
fftam <- subset(fftam, subset = nFeature_RNA > 200 & percent.mt < 15)
fftam <- ScaleData(object = fftam, vars.to.regress = c("percent.ribo","Rn45s"))
fftam <- NormalizeData(fftam, verbose = FALSE)
fftam <- FindVariableFeatures(fftam, selection.method = "vst", nfeatures = 2000)


ff.anchors <- FindIntegrationAnchors(object.list = list(ffoil, fftam), dims = 1:20)
to_integrate <- Reduce(intersect, lapply(ff.anchors@object.list, rownames))
ff.combined <- IntegrateData(anchorset = ff.anchors, features.to.integrate = to_integrate, dims = 1:20)

DefaultAssay(ff.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
ff.combined <- ScaleData(ff.combined, verbose = FALSE)
ff.combined <- RunPCA(ff.combined, npcs = 30, verbose = FALSE)

ElbowPlot(ff.combined)

library(reticulate)
use_python(python = "C:/Users/karim/Anaconda3", required = TRUE)

# UMAP and Clustering
ff.combined <- RunUMAP(ff.combined, reduction = "pca", dims = 1:19)
ff.combined <- FindNeighbors(ff.combined, reduction = "pca", dims = 1:19)
ff.combined <- FindClusters(ff.combined, resolution = 0.75)

# Visualization
p1 <- DimPlot(ff.combined, reduction = "umap", group.by = "tam")
p2 <- DimPlot(ff.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
DimPlot(ff.combined, reduction = "umap", split.by = "tam", label = TRUE)


table(Idents(ff.combined), ff.combined$tam)

# t-SNE and Clustering
ff.combined <- RunTSNE(ff.combined, reduction = "pca", dims = 1:19)
ff.combined <- FindNeighbors(ff.combined, reduction = "pca", dims = 1:19)
ff.combined <- FindClusters(ff.combined, resolution = 0.75)

# Visualization
p1 <- TSNEPlot(ff.combined, group.by = "tam")
p2 <- TSNEPlot(ff.combined, label = TRUE)
plot_grid(p1, p2)
TSNEPlot(ff.combined, split.by = "tam",label = TRUE)


table(Idents(ff.combined), ff.combined$tam)


FeaturePlot(ff.combined, features = c("Mki67"), min.cutoff = "q9")
FeaturePlot(ff.combined, features = c("Mki67"), split.by = "tam", max.cutoff = 3, cols = c("gray", "red"))

ff.combined.markers <- FindAllMarkers(ff.combined, min.pct = 0.5, logfc.threshold = 0.25)
ff.combined.markers %>% group_by(cluster)

write.csv(ff.combined.markers, "Cluster_markers_fc0.25 190822 15 19 0.75 del ribo.csv")




#Average gene expresssion
gene.exp.average <- AverageExpression(ff.combined)
write.csv(gene.exp.average[["RNA"]], "Average_gene_exp_ff.csv")
write.csv(gene.exp.average[["integrated"]], "Average_gene_exp_ff_integrated.csv")

DoHeatmap(ff.combined, features = "Isg15", angle = 0) + NoL
