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

writeLines("3. Executing CLR algorithm...")
CLR_net <- clr(mim)
CLR_g <- graph.adjacency(CLR_net, weighted=TRUE)
edgeList <- get.data.frame(CLR_g)
edgeListSorted <- edgeList[with(edgeList, order(-weight)),]

writeLines("...Generating output file - CLR_predictions.txt")
write.table(edgeListSorted, "CLR_predictions.txt", quote = FALSE,sep="\t", row.names = FALSE, col.names = FALSE);

