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

writeLines("4. Executing MRNETB algorithm...")
mrnetB_net<-mrnetb(mim)
mrnetB_g <- graph.adjacency(mrnetB_net,weighted=TRUE)
edgeList <- get.data.frame(mrnetB_g)
edgeListSorted <- edgeList[with(edgeList, order(-weight)),]

writeLines("...Generating output file - MRNETB predictions")
write.table(edgeListSorted,"MRNETB_predictions.txt", quote = FALSE,sep="\t", row.names = FALSE, col.names = FALSE);

writeLines("5. Executing GENIE algorithm...\n")
source("GENIE3_R/GENIE3_R/genie3.R")
expr.matrix <- read.expr.matrix(args[1], form="rows.are.samples")
weight.matrix <- get.weight.matrix(expr.matrix)
link.list <- get.link.list(weight.matrix)
writeLines("...Generating output file - GENIE3.txt")
write.table(link.list,"GENIE3_predictions.txt", quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE);


writeLines("6. Executing Distance-Correlation predictions...\n")
p <- ncol(inputData) # p stands for the number of genes

#compute distance correlations for each pair of genes
#optional: use snowfall package for parallel computing if the gene expression file is big

sfInit(parallel=T,cpus=8)
sfLibrary(energy)
sfExport("inputData")
sfExport("p")

dcor2<-function(i){
  dis<-matrix(0,ncol=p)
  for(j in (i+1):p){
    dis[j]<-dcor(inputData[,i],inputData[,j])
  }
  return(dis)
}

res<-sfSapply(1:(p-1),dcor2)
res<-cbind(res,rep(0,p))
res<-res+t(res)

colnames(res) <- colnames(inputData) # replacing the column headers of the matrix with gene names

#res will be a matrix whose entry(i,j) stands for the distance correlation of gene i and gene j. 
#res is the ouput adjacency matrix from Distance Correlation algo.
g <- graph.adjacency(res, weighted=TRUE)
edgeList <- get.data.frame(g)
edgeListSorted <- edgeList[with(edgeList, order(-weight)),]

writeLines("\n...Generating output file - DC_predictions.txt\n")
write.table(edgeListSorted,"DC_predictions.txt", quote = FALSE,sep="\t", row.names = FALSE, col.names = FALSE);

sfStop()
writeLines("** R script completed ***\n")

