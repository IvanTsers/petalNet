
#' Finds all maximal cliques
#'  
#' This function calls the maximal.cliques function from igraph version 0.7
#'
#' @param adjMatrix a network matrix with 0 and 1 indicating none- and existing links, respectively
#' 
#' @return a list of containing numeric vetors of vertex ids; each list element is a clique
#'
#' @export

findMaximalCliques <- function(adjMatrix){
   library("igraph")
   smallG = graph.adjacency(adjMatrix,mode="undirected", weighted=NULL, diag=FALSE, add.colnames=NULL)
   MC = maximal.cliques(smallG)
   return(MC)
}


