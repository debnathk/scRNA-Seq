library(dplyr)
library(igraph)
library(plotrix)

# Put file-path here 
setwd("D:/scRNA-Seq_data/final_schwartz_networks")
network_files = list.files(pattern = "\\_final.tsv")
node_files = list.files(pattern = "\\_selected.csv")

implants <- list('mSLA', 'OM', 'PT', 'TCPS', 'Titan')

# Vizualize network for top 10 high degree nodes for each implant
for (nw_file in network_files){
  for (no_file in node_files){
    if (grepl('mSLA', nw_file) & grepl('mSLA', no_file)) {
      network_graph <- read.csv(nw_file, sep = "\t", header = FALSE)
      node_graph <- read.csv(no_file, header = TRUE)
      colnames(network_graph) <- c("Edge1", "Edge2", "Weight")
      top_5_nodes = head(node_graph[, 1], 5)
      # print(top_5_nodes)
      processed_graph <- subset(network_graph, Edge1 %in% top_5_nodes)
      edge_list <- processed_graph %>% select(Edge1, Edge2)

      # Creates the igraph object
      igraph_obj <- graph_from_edgelist(as.matrix(edge_list), directed = TRUE)
      
      # Set labels for top 5 nodes
      # V(igraph_obj)$label <- ifelse(V(igraph_obj)$name %in% top_5_nodes, V(igraph_obj)$name, NA)
      
      # Assign different colors to the top 5 nodes
      node_colors <- rainbow(length(top_5_nodes))
      
      # Set colors for all nodes in the graph
      V(igraph_obj)$color <- ifelse(V(igraph_obj)$name %in% top_5_nodes,
                                    node_colors[match(V(igraph_obj)$name, top_5_nodes)], "gold")
      V(igraph_obj)$label <- NA
      # Set the size of the top nodes to 10
      V(igraph_obj)$size <- ifelse(V(igraph_obj)$name %in% top_5_nodes, 10, 5)
      
      # Saves the plot
      file_sub = sub(".tsv", "", nw_file)
      plot_name = paste("D:/scRNA-Seq_data/final_schwartz_networks/", file_sub, "_top_5.png", sep = "")
      png(plot_name, width = 4800, height = 3200)
      # Creates the plot
      plot(igraph_obj, edge.arrow.size=.5,
           vertex.frame.color="black",edge.curved=0.2)
      # Close the device to save the plot
      # Create legend for top 5 nodes with their names
      legend("topright", legend = top_5_nodes, fill = node_colors, title = "Top Nodes", cex = 5)
      dev.off()
    }
  }
}
