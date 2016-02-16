
#' Calculating all vertices' cluster coefficient and degree
#'
#' Based on a network matrix, for each vertex of the network/row of the matrix, its cluster coeffient and connectivity is calculated and return in table format 
#'
#' @param adjMatrix binary adjacency matrix, diagonal is assumed to be 1
#' @param totxtFile whether to right the table to a separate text file, boolean variable TRUE or FALSE(default)
#'
#' @return matrix with colnames "Index", "Gene","LocalCC","k"
#'
#' @export

getVertexStats <-function(adjMatrix,totxtFile=FALSE){

	#converting adjacency matrix to a graph format for igraph to use it
	library("igraph")
	GraphAdjM = graph.adjacency(adjMatrix, mode="undirected", weighted=NULL, diag=FALSE, add.colnames=NULL)
	LocalCC = round(transitivity(GraphAdjM, type="local"), digits=3)
	rm(GraphAdjM)

	## node stats
	Index=seq(1:dim(adjMatrix)[1])
	k_net=rowSums(adjMatrix)-1

	nodeStatsT=cbind(Index,rownames(adjMatrix),LocalCC,k_net)
	colnames(nodeStatsT) = c("Index", "Gene","LocalCC","k")

	if(totxtFile){write.table(nodeStatsT,file="NetworkVertexStats.txt",sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)}

	return(nodeStatsT)
}

