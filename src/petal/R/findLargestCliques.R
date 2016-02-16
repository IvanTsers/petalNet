
#' Finds all largest cliques
#'
#' This function calls the largest.cliques function from igraph version 0.7
#'
#' @param adjMatrix a network matrix with 0 and 1 indicating none- and existing links, respectively
#'
#' @return a numeric matrix; each column is a largest clique containing the vertex ids
#'
#' @export

findLargestCliques <- function(adjMatrix){
	library("igraph")
	smallG = graph.adjacency(adjMatrix,mode="undirected", weighted=NULL, diag=TRUE, add.colnames=NULL)
	LC = largest.cliques(smallG)

	LCM = NULL
	for(i in 1:length(LC)){
		x   = as.vector(LC[[i]])
		LCM = cbind(LCM, x)	
	}			
	return(LCM)
}


