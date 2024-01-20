# https://satijalab.org/seurat/articles/multimodal_vignette

# Dive into the scientific groove where single-cell genomics becomes the ultimate
# multitasker! Picture this: cells revealing their secrets in a symphony of data types,
# dancing together in a sensational performance of multimodal analysis.
# It's not just science; it's a cellular carnival of discovery! 🧬🎉 #GenomicsParty

# Example: CITE-seq enables simultaneous measurement of transcriptomes and cell-surface
# proteins from the same cell.

# Cell Hashing is a groundbreaking technique that tackles challenges in single-cell
# genomics head-on! Imagine using special oligo-tagged antibodies that act like
# personalized ID tags for cells. These tags, targeting universally expressed surface
# proteins, give each cell a unique label. Now, here's the magic: you can pool these
# cells together and, by sequencing the tags alongside their transcriptome, unveil
# their individual identities. It's like giving each cell a passport for its origin!
# This not only helps overcome batch effects and detect multiplets but also
# supercharges commercial droplet-based systems, making experiments more cost-effective.
# It's genomics with a touch of spy intrigue – Cell Hashing, where cells reveal their
# secrets with style! 🕵️‍♂️🔍🧬 #CellHashingRevolution

library(Seurat)
library(ggplot2)
library(patchwork)

# Load in the RNA UMI matrix
# Note that this dataset also contains ~5% of mouse cells, which we can use as negative
# controls for the protein measurements. For this reason, the gene expression matrix has
# HUMAN_ or MOUSE_ appended to the beginning of each gene.
cbmc.rna <- as.sparse(read.csv(file = "C:/Users/debnathk/Desktop/study/scRNA-Seq/data/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz",
                               sep = ",", header = T, row.names = 1))

# To make life a bit easier going forward, we're going to discard all but the top 100 most
# highly expressed mouse genes, and remove the 'HUMAN_' from the CITE-seq prefix
cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna)

# Load in the ADT UMI matrix
cbmc.adt <- as.sparse(read.csv(file = "C:/Users/debnathk/Desktop/study/scRNA-Seq/data/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz",
                                           sep = ",", header = T, row.names = 1))

# Note that since measurements were made in the same cells, the two matrices have identical
# column names
all.equal(colnames(cbmc.rna), colnames(cbmc.adt))


# Setup a Seurat object, add the RNA and protein data
# creates a Seurat object based on the scRNA-seq data
cbmc <- CreateSeuratObject(counts = cbmc.rna)

# We can see that by default, the cbmc object contains an assay storing RNA measurement
Assays(cbmc)

# create a new assay to store ADT information
adt_assay <- CreateAssay5Object(counts = cbmc.adt)

# add this assay to the previously created Seurat object
cbmc[["ADT"]] <- adt_assay

# Validate that the object now contains multiple assays
Assays(cbmc)

# Cluster cells on the basis of their scRNA-seq profiles
# perform visualization and clustering steps
cbmc <- NormalizeData(cbmc)
cbmc <- FindVariableFeatures(cbmc)
cbmc <- ScaleData(cbmc)
cbmc <- RunPCA(cbmc, verbose = F)
cbmc <- FindNeighbors(cbmc, dims = 1:30)
cbmc <- FindClusters(cbmc, resolution = 0.8, verbose = F)
cbmc <- RunUMAP(cbmc, dims = 1:30)
DimPlot(cbmc, label = T)


# Visualize multiple modalities side-by-side
# Normalize ADT data,
cbmc <- NormalizeData(cbmc, normalization.method = 'CLR', margin = 2, assay = 'ADT')

# Now, we will visualize CD19 levels for RNA and protein By setting the default assay, we can
# visualize one or the other
DefaultAssay(cbmc) <- 'ADT'
p1 <- FeaturePlot(cbmc, "CD19", cols = c('lightgreen', 'darkgreen')) + ggtitle("CD19 protein")
DefaultAssay(cbmc) <- 'RNA'
p2 <- FeaturePlot(cbmc, "CD19") + ggtitle("CD19 RNA")

# Show plots side-by-side
p1 | p2

# Alternate way - using Assay key
p1 <- FeaturePlot(cbmc, 'adt_CD19', cols = c("lightgreen", "darkgreen")) + ggtitle("CD19 protein")
p2 <- FeaturePlot(cbmc, 'rna_CD19') + ggtitle("CD19 RNA")
p1 | p2
