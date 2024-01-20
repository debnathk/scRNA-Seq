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

