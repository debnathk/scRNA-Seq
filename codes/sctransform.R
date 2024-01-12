# sctransform:
# - A modelling framework for normalization and variance stabilization of
# molecular count data from scRNA-Seq experiments, diminishing the need for
# heuristic steps including pseudocount addition, log-transformation and improves
# common downstream analytical tasks such as variable gene selection, dimensional
# reduction, and differential gene expression.

# Load libraries
library(Seurat)
library(ggplot2)
library(sctransform)

# Load data and create seurat object
pbmc_data <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19")
pbmc <- CreateSeuratObject(counts = pbmc_data)


# Apply sctransform normalization
# Store mitochondrial percentage in object meta data
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")

# Run sctransform
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = F)


# Perform dimensionality reduction using UMAP and PCA embedding
# These are now standard steps in the Seurat workflow for visualization and clustering
pbmc <- RunPCA(pbmc, verbose = F)
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = F)

pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = F)
pbmc <- FindClusters(pbmc, verbose = F)
DimPlot(pbmc, label = T)

