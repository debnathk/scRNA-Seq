# sctransform:
# - A modelling framework for the normalization and variance stabilization of
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