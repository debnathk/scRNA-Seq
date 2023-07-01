library(dplyr)
library(Seurat)
library(patchwork)
library(pacman)

# Load the lungs dataset
# Dataset link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132771
lungs.data <- Read10X(data.dir = "C:/Users/Bumba/Desktop/STUDY/phd@vcu/scRNA-Seq/data/NML1")
# Initialize the Seurat object with the raw (non-normalized data).
lungs <- CreateSeuratObject(counts = lungs.data, project = "lungs", min.cells = 3, min.features = 200)
lungs

# Pre-processing ====
# QC ====
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
lungs[["percent.mt"]] <- PercentageFeatureSet(lungs, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(lungs, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(lungs, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(lungs, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

lungs <- subset(lungs, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalization ====
# lungs <- NormalizeData(lungs, normalization.method = "LogNormalize", scale.factor = 10000)
lungs <- NormalizeData(lungs)

# Identification of highly variable features
lungs <- FindVariableFeatures(lungs, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(lungs), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(lungs)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot1 + plot2

# Scaling ====
all.genes <- rownames(lungs)
lungs <- ScaleData(lungs, features = all.genes)

# Linear Dimension reduction ====
lungs <- RunPCA(lungs, features = VariableFeatures(object = lungs))
# Examine and visualize PCA results a few different ways
print(lungs[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(lungs, dims = 1:2, reduction = "pca")
DimPlot(lungs, reduction = "pca")

# Heatmap ====
DimHeatmap(lungs, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(lungs, dims = 1:15, cells = 500, balanced = TRUE)

# Determining dimensionality of the datset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
lungs <- JackStraw(lungs, num.replicate = 100)
lungs <- ScoreJackStraw(lungs, dims = 1:20)
JackStrawPlot(lungs, dims = 1:15)

ElbowPlot(lungs)

# Clustering ====
lungs <- FindNeighbors(lungs, dims = 1:10)
lungs <- FindClusters(lungs, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(lungs), 5)

# UMAP/tSNE ====
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
lungs <- RunUMAP(lungs, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(lungs, reduction = "umap")

# Save object ====
saveRDS(lungs, file = "C:/Users/Bumba/Desktop/STUDY/phd@vcu/scRNA-Seq/output_files/lungs_tutorial.rds")

# Finding differentially expressed features (cluster biomarkers) ====
# find all markers of cluster 2
cluster2.markers <- FindMarkers(lungs, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 20)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(lungs, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
lungs.markers <- FindAllMarkers(lungs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
lungs.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
cluster0.markers <- FindMarkers(lungs, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(lungs, features = c("ADIRF", "NGFRAP1"))
# you can plot raw counts as well
VlnPlot(lungs, features = c("NUPR1", "MARCO"), slot = "counts", log = TRUE)

FeaturePlot(lungs, features = c("SFTA2", "SFTPB", "NAPSA", "HOPX", "PEBP4", "RGS2", 
                               "AREG", "CORO1A", "IL1R2"))

lungs.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(lungs, features = top10$gene) + NoLegend()

# Clear up ====
rm(list = ls())
p_unload(dplyr, Seurat, patchwork)
q()
