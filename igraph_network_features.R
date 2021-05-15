coExpressNetwork <- function(corMat){
  
  corNet = igraph::graph_from_adjacency_matrix(corMat, mode = "lower", weighted = TRUE,diag = F, add.colnames=T, add.rownames=T)
  corNet = set_vertex_attr(corNet,'name',value = rownames(corMat))
  # corNet = delete_vertex_attr(corNet,'TRUE')
  corNet = igraph::simplify(corNet, remove.multiple = T, remove.loops = T)
  return(corNet)
}

# input network, output annotation of network

networkNodeAnno <- function(corNet){
  
  degree <- igraph::degree(corNet)  

  betweenness <- round(centralization.betweenness(corNet)$res,3) 

  closeness <- round(centralization.closeness(corNet)$res,3)   
  
  netAnno = data.frame(node=names(degree), degree, betweenness, closeness)
  return(netAnno)
  
}