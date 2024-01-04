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
# The above code randomly assigns cells to a new meta.data column to use as
# an example in demonstrating plotting functions that can split or group plots by
# another valid meta.data variable.
features <- c("CCL5", "IL32", "LYZ", "PTPRCAP", "FCGR3A", "PF4")
pbmc3k.final


# Five visualizations of marker feature expression
# Ridge plots - to visualize single-cell expression distributions in each cluster
RidgePlot(pbmc3k.final, features = features, ncol = 2)


# Violin plots - to visualize single-cell expression distributions in each cluster
VlnPlot(pbmc3k.final, features = features)


# Feature plot - to visualize feature expression in low-dimensional space
FeaturePlot(pbmc3k.final, features = features)
# Plot a legend to map colors to expression levels
FeaturePlot(pbmc3k.final, features = "MS4A1")
# Adjust the contrast in the plot
FeaturePlot(pbmc3k.final, features = "MS4A1", min.cutoff = 1, max.cutoff = 3)
# Calculate feature-specific contrast levels based on quantiles of non-zero expression
# Particularly useful when plotting multiple markers
FeaturePlot(pbmc3k.final, features = c("MS4A1", "PTPRCAP"), min.cutoff = "q10",
            max.cutoff = "q90")
# Visualize co-expression of two features simultaneously
FeaturePlot(pbmc3k.final, features = c("MS4A1", "CD79A"), blend = TRUE)
# Split visualization to view expression by groups (replaces FeatureHeatMap)
FeaturePlot(pbmc3k.final, features = c("MS4A1", "CD79A"), split.by = "groups")

# Dot plot - the size of the dot corresponds to the percentage of cells expressing
# the feature in each cluster. The color represents the average expression level
DotPlot(pbmc3k.final, features = features) + RotatedAxis()


# Single-cell heatmap of feature expression
DoHeatmap(subset(pbmc3k.final, downsample = 100), features = features, size = 3)



# Updated and expanded visualization functions
# Violin plots can also be split on some variable, by adding the splitting variable
# to object metadata and pass it to the split.by argument
VlnPlot(pbmc3k.final, features = "percent.mt", split.by = "groups")

# Dot plot with split by groups
DotPlot(pbmc3k.final, features = features, split.by = "groups") + RotatedAxis()

# Dimplot replaces TSNEplot, PCAplot - it will plot umap, tsne or pca by default 
# in that order
DimPlot(pbmc3k.final)
# Plotting PCA by nullifying UMAP
pbmc3k.final.no.umap <- pbmc3k.final
pbmc3k.final.no.umap[['umap']] <- NULL
DimPlot(pbmc3k.final.no.umap) + RotatedAxis()

#DoHeatMap with a grouping bar, splitting the heatmap into groups or clusters
DoHeatmap(pbmc3k.final, features = VariableFeatures(pbmc3k.final)[1:100], 
          cells = 1:500, size = 4, angle = 90) + NoLegend()






