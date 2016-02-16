
#' Internal function of downstreamAnalysis, calculates the density of a network
#'
#' The density of the passed network is calculated (number of edges in the graph versus number of all possible edges). 
#'
#' @param FriendsM subnetwork in matrix format; diagonal must be set to 1
#' @param FileOut name of file to write to if desired
#' @param outReturn whether to return the calculated density, boolean variable TRUE(default) or FALSE
#'
#' @return decimal between 0 and 1; zero indicating there are no edges in the network, and a network with density 1 is a clique (completely connected network).

densityM <- function(FriendsM, FileOut=NULL, outReturn=TRUE){
	n = dim(FriendsM)[1]
	numEdges = (sum(rowSums(FriendsM))-n)/2		# number of edges in the matrix
	cliqueEdges = ((n-1)*n)/2			# number of edges in matching perfect clique
	
	if(length(FileOut)!=0){
	toFile = paste("This subnetwork has a density of: ", round(numEdges/cliqueEdges,digits=4), sep="")
   	writeToFile(toFile, FileOut, TRUE)
	rm(toFile)

	 if(isTRUE((cliqueEdges-numEdges)!=0)){
		if(isTRUE((cliqueEdges-numEdges)==1)){
		      toFile = paste("In the ", n, "-node subnetwork ", cliqueEdges-numEdges, " edge is missing to be a clique", sep="")
   	              writeToFile(toFile, FileOut, TRUE)
		}else{toFile = paste("In the ", n, "-node subnetwork ", cliqueEdges-numEdges, " edges are missing to be a clique", sep="")
   	              writeToFile(toFile, FileOut, TRUE)
		}
	 }
	}
	if(outReturn){return(numEdges/cliqueEdges)}
}


