rm(list = ls())
BiocManager::install("fgsea")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Rn.eg.db")
library(fgsea)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library("DOSE")
setwd("C:/Users/debnathk/Desktop/study/scRNA-Seq/Volcano_data/Volcano_sig")

files = list.files(pattern = "\\.csv$")

for (file in files) {
  # Read in the rank_table as a csv
  # file = "Titan_TCPS_sig.csv"
  tryCatch({
    rank_table <- read.csv(file)
    
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
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
    
    # Create the figure 
    file_sub = sub(".csv", "", file)
    overrepresented_figure <- barplot(ego, showCategory=10, title = paste("Overrepresented GO MFs in ", file_sub, sep = ""))
    
    # Save figure
    file_name <- paste("../results/", file_sub, "_overrepresented_MF_figure.png", sep = "")
    ggsave(file_name, plot = overrepresented_figure, width = 10, height = 6, units = "in", dpi = 300)
    
  }, error = function(e) {
    # Print the error message without stopping the loop
    cat("Error:", conditionMessage(e), "\n")
  })
  rm(list=ls())
  # Clear console
  cat("\014")
}