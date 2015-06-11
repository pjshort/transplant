# functions for clustering and cutting gene expression data

library(WGCNA)

build_tree <- function(expr_data, soft_thresh_power){
  
  # build a dendrogram from expression data using topological overlap matrix
  # returns dendrogram and distance matrix (one minus topological overlap)
    
  # calculate the adjacency matrix
  write("Calculating gene-gene adjacency matrix.", stderr())
  adjacency = adjacency(expr_data, power = soft_thresh_power);
  
  # adjacency matrix -> topological overlap matrix
  write("Converting adjacency matrix to Topological Overlap Matrix.", stderr())
  
  TOM = TOMsimilarity(adjacency)
  dissTOM = 1-TOM # TOM dissimilarity
  
  # create dendrogram with hclust
  write("Creating dendrogram from gene-gene TOM.", stderr())
  
  geneTree = hclust(as.dist(dissTOM), method = "average")
  
  return(list("tree" = geneTree, "dissTOM" = dissTOM))
  
}