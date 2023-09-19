library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset 
pbmc.data <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# QC and cell selection

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

## Visualize QC metrics as violin plots
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## Feature scatter Plot
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalizing the data (Using LogNormalize method)

pbmc <- NormalizeData(pbmc)

# Feature selection

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

## Identify top 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

## Plot variable features with and without features
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = T)
plot1 + plot2

# Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Linear Dimension Reduction (PCA)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

# Determining the dimensionality of the dataset
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)

# Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

head(Idents(pbmc), 5)

# Non-Linear dimensional reduction (UMAP/tSNE)
pbmc <- RunUMAP(pbmc, dims =  1:10)
DimPlot(pbmc, reduction = "umap")

## Saving the object
saveRDS(pbmc, file = "pbmc_tutorial.rds")

# Finding differentially expressed features (cluster biomarkers)

## find all marker of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

## find all clusters distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

## find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% 
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster0.markers, n = 5)

## VlnPlot(pbmc.markers, features = c("RPS12, RPS6"))

# Assigning cell type identity to clusters
