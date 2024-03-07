rm(list = ls())
BiocManager::install("fgsea")
library(fgsea)
library(dplyr)
library(clusterProfiler)
library("DOSE")
setwd("/Users/lodimk2/Documents/schwartz_analysis")

# Read in the rank_table as a csv
rank_table <- read.csv("Volcano_data/Volcano_sig/Titan_TCPS_sig.csv")

# Sort the table by adjusted pvalue
rank_table <- rank_table %>% arrange((padj))

# Use padj cutoff of 0.01, this is customizable based on specific requirements. 
rank_table <- rank_table %>% filter(padj < 0.01)

# convert gene symbols to Entrez ID
# NOTE: species is rat, so convert from human Db to rat Db. 
eg = bitr(rank_table$Gene.names, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Rn.eg.db")

# Run enrichGO for overrepresented pathways test. 
ego <- enrichGO(gene          = eg$ENTREZID,
                OrgDb = "org.Rn.eg.db",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

# Create the figure 

overrepresented_figure <- barplot(ego, showCategory=10, title = "Overrepresented GO BPs in Titan_TCPS")

# Save figure
ggsave("Titan_TCPS_overrepresented_figure.png", plot = overrepresented_figure, width = 10, height = 6, units = "in", dpi = 300)


