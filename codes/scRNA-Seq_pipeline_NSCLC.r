# Tutorial : https://biostatsquid.com/scrnaseq-preprocessing-workflow-seurat/

# Step 1: Set up your script and load the data
# Setting up environment ===================================================

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(Seurat)
library(Matrix)

# Set input path
path <- "C:/Users/Bumba/Desktop/STUDY/phd@vcu/scRNA-Seq"
setwd(path)

list.files(path)
set.seed(42)

# Import data ===================================================
nsclc_sm <- Read10X_h5(filename = "data/40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_count_raw_feature_bc_matrix.h5")
str(nsclc_sm) # Check the multiple modalities (list of matrixes) - we're interested in Gene expression
cts <- nsclc_sm$`Gene Expression`

#Step 2: Create your Seurat object
# Initialise Seurat object (raw, non-normalised data)
nsclc_seu <- CreateSeuratObject(counts = cts, project = 'NSCLC', min.cells = 3, min.features = 200)
str(nsclc_seu)

#Step 3: Quality control
# QC ===================================================
## perc_mt -----------------------
nsclc_seu[['percent_mt']] <- PercentageFeatureSet(nsclc_seu, pattern = '^MT-')
View(nsclc_seu@meta.data)
#colnames(nsclc_seu@meta.data)
VlnPlot(nsclc_seu, features = c("nCount_RNA", "nFeature_RNA", "percent_mt"), ncol = 3)
FeatureScatter(nsclc_seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

## Filtering -----------------------
nsclc_seu <- subset(nsclc_seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent_mt < 5)
View(nsclc_seu@meta.data)

#Step 4: Normalisation
## Normalisation -----------------------
#nsclc_seu <- NormalizeData(nsclc_seu, normalization.method = 'LogNormalize', scale.factor = 10000)
nsclc_seu <- NormalizeData(nsclc_seu)
str(nsclc_seu@meta.data)

#Step 5: Identification of highly-variable genes
## Identify highly-variable features ---------------------------
nsclc_seu <- FindVariableFeatures(nsclc_seu, selection.method =  'vst', nfeatures = 2000)

# Identify the top 10 HVGs
top10 <- head(VariableFeatures(nsclc_seu), 10)
top10_plot <- VariableFeaturePlot(nsclc_seu)
LabelPoints(plot = top10_plot, points = top10, repel = TRUE, xnudge =0, ynudge = 0)

#Step 6: Scaling
# Scaling ==================================================
all_genes <- rownames(nsclc_seu)
nsclc_seu <- ScaleData(nsclc_seu, features = all_genes)
View(nsclc_seu@assays$RNA)

#Step 7: PCA
# PCA ===================================================
nsclc_seu <- RunPCA(nsclc_seu, features = VariableFeatures(nsclc_seu))
print(nsclc_seu[['pca']], dims = 1:5, nfeatures = 5)
DimHeatmap(nsclc_seu, dims = 1, cells = 500, balanced = TRUE)
DimPlot(nsclc_seu, reduction = "pca") + NoLegend()
# determine dimensionality of the data
ElbowPlot(nsclc_seu)

#Step 8: Clustering
# Clustering ===================================================
nsclc_seu <- FindNeighbors(nsclc_seu, dims = 1:15)
nsclc_seu <- FindClusters(nsclc_seu, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
View(nsclc_seu@meta.data)
DimPlot(nsclc_seu, group.by = 'RNA_snn_res.1', label = TRUE)
Idents(nsclc_seu) <- 'RNA_snn_res.0.1' # set identity of clusters

#Step 9: UMAP
# UMAP ===================================================
nsclc_seu <- RunUMAP(nsclc_seu, dims = 1:15)
DimPlot(nsclc_seu, reduction = 'umap')

#Step 10: Export and save your Seurat object
# Save it! ===================================================
saveRDS(nsclc_seu, file = 'nsclc_seu.RDS')


