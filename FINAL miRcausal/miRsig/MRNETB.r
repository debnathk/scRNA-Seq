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

writeLines("4. Executing MRNETB algorithm...")
mrnetB_net<-mrnetb(mim)
mrnetB_g <- graph.adjacency(mrnetB_net,weighted=TRUE)
edgeList <- get.data.frame(mrnetB_g)
edgeListSorted <- edgeList[with(edgeList, order(-weight)),]

writeLines("...Generating output file - MRNETB predictions")
write.table(edgeListSorted,"MRNETB_predictions.txt", quote = FALSE,sep="\t", row.names = FALSE, col.names = FALSE);

