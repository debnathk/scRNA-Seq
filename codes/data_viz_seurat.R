# install.packages("remotes")
# remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
# remotes::install_github("satijalab/seurat-data")
# SeuratData::InstallData("pbmc3k")

# Load libraries
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)

pbmc3k.final <- LoadData("pbmc3k", type = "pbmc3k.final")
pbmc3k.final$groups <- sample(c("group1", "group2"), size = ncol(pbmc3k.final),
                              replace = TRUE)
features <- c("CCL5", "IL32", "LYZ", "PTPRCAP", "FCGR3A", "PF4")
pbmc3k.final


# Five visualizations of marker feature expression
# Ridge plots - to visualize single-cell expression distributions in each cluster
RidgePlot(pbmc3k.final, features = features, ncol = 2)

# Violin plots - to visualize single-cell expression distributions in each cluster
VlnPlot(pbmc3k.final, features = features)

# Feature plot - to visualize feature expression in low-dimensional space
FeaturePlot(pbmc3k.final, features = features)




