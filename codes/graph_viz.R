library(dplyr)
library(igraph)

# Put file-path here 
full_network_graph <- read.csv("/Users/lodimk2/Documents/alcoholism_miRNA/final_network/mirsig_network_final.tsv", sep = "\t", header = FALSE)
colnames(full_network_graph) <- c("Edge1", "Edge2", "Weight")
processed_graph <- full_network_graph %>% filter(Weight > 0.9)

edge_list <- processed_graph %>% select(Edge1, Edge2)

# Creates the igraph object
igraph_obj <- graph_from_edgelist(as.matrix(edge_list), directed = TRUE)

# Saves the plot
png("full_net_viz.png", width = 1200, height = 800)
# Creates the plot
plot(igraph_obj, edge.arrow.size=.5, vertex.color="gold", vertex.size=10, vertex.frame.color="black", vertex.label.color="black", vertex.label.cex=0.8, vertex.label.dist=3, edge.curved=0.2, edge.arrow.size=.4)
# Close the device to save the plot
dev.off()


# Create sorted degrees dataframe
degrees <- as.data.frame(degree(igraph_obj))
colnames(degrees) <- "degree"
sorted_degree <- degrees %>% arrange(desc(degree))

write.csv(sorted_degree, "sorted_degrees_{sample_name_here}.csv")
