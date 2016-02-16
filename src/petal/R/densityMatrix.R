
#' Calculates the density of a provided network
#'
#' The density of the passed network is calculated (number of edges in the graph versus number of all possible edges). 
#'
#' @param FriendsM a network matrix with 0 and 1 indicating none- and existing links, respectively; diagonal must be set to 1
#' @param moreInfo whether to return more information, boolean variable TRUE or FALSE(default)
#'
#' @return decimal between 0 and 1; zero indicating there are no edges in the network, and a network with density 1 is a clique (completely connected network).
#'
#' @export

densityMatrix <- function(FriendsM, moreInfo=FALSE){
	n = dim(FriendsM)[1]
	numEdges = (sum(rowSums(FriendsM))-n)/2		# number of edges in the matrix
	cliqueEdges = ((n-1)*n)/2			# number of edges in corresponding clique

	if(moreInfo){
	 if(isTRUE((cliqueEdges-numEdges)!=0)){
		if(isTRUE((cliqueEdges-numEdges)==1)){
		      	print(paste("In the ", n, "-node subnetwork ", cliqueEdges-numEdges, " edge is missing to be a clique", sep=""))
		}else{
			print(paste("In the ", n, "-node subnetwork ", cliqueEdges-numEdges, " edges are missing to be a clique", sep=""))
		}
	 }
	}

	return(numEdges/cliqueEdges)
}

