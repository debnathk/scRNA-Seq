library(igraph)
library(minet)
library(energy)
library(snowfall)

writeLines("\nExecuting this R script\n")
args = commandArgs(trailingOnly = TRUE)
if (length(args)==0){
 stop("Atleast one argument must be supplied (input file).n",call.=FALSE)
}

writeLines("Loading the expression matrix")
inputData <- read.delim(args[1], header = TRUE, sep = "\t")

mim<-build.mim(inputData, estimator="spearman")

writeLines("5. Executing GENIE algorithm...\n")
source("GENIE3_R/GENIE3_R/genie3.R")
expr.matrix <- read.expr.matrix(args[1], form="rows.are.samples")
weight.matrix <- get.weight.matrix(expr.matrix)
link.list <- get.link.list(weight.matrix)
writeLines("...Generating output file - GENIE3.txt")
write.table(link.list,"GENIE3_predictions.txt", quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE);


