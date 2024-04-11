library(dplyr)
library(igraph)

# Put file-path here 
setwd("D:/scRNA-Seq_data/final_schwartz_networks")
files = list.files(pattern = "\\_final.tsv")

for (file in files){
  network_graph <- read.csv(file, sep = "\t", header = FALSE)
  colnames(network_graph) <- c("Edge1", "Edge2", "Weight")
  processed_graph <- network_graph %>% filter(Weight > 0.975)
  
  edge_list <- processed_graph %>% select(Edge1, Edge2)
  
  # Creates the igraph object
  igraph_obj <- graph_from_edgelist(as.matrix(edge_list), directed = TRUE)
  
  # Saves the plot
  file_sub = sub(".tsv", "", file)
  plot_name = paste("D:/scRNA-Seq_data/final_schwartz_networks/", file_sub, "_0975.png", sep = "")
  png(plot_name, width = 4800, height = 3200)
  # Creates the plot
  plot(igraph_obj, edge.arrow.size=.5, vertex.color="gold", vertex.size=3,
       vertex.frame.color="black", vertex.label.color="black",
       vertex.label.cex=1, vertex.label.dist=3, edge.curved=0.2)
  # Close the device to save the plot
  dev.off()

  # Create sorted degrees dataframe
  degrees <- as.data.frame(degree(igraph_obj))
  colnames(degrees) <- "degree"
  sorted_degree <- degrees %>% arrange(desc(degree))
  degree_file_sub = sub(".tsv", "", file)
  degree_file_name = paste("D:/scRNA-Seq_data/final_schwartz_networks/sorted_degrees_",
                           file_sub, "_0975.csv", sep = "")
  write.csv(sorted_degree, degree_file_name)
}