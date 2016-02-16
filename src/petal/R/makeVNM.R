
#' Creating a vicinity network martrix 
#'
#' Creats a vicinity network based on one or multiple node. If created on multiple nodes the input nodes must be a clique in the network
#'
#' @param currentSubgroup a single node identifier or a string of node identifiers (rownames of adjMartrix)
#' @param adjMatrix boolean network matrix
#'
#' @return vicinity network of input node(s)
#'
#' @export

makeVNM <- function(currentSubgroup, adjMatrix){

   nodeIndex = match(currentSubgroup, rownames(adjMatrix))
   tmpM      = adjMatrix[,nodeIndex]

   if(isTRUE(length(nodeIndex)==1)){
   	 Friends   = as.vector(which(tmpM==1))
   }else{Friends   = as.vector(which(rowSums(tmpM)==length(nodeIndex)))}

   FriendsM  = adjMatrix[Friends,Friends]
	
   if(length(Friends)==0){
	   print("The gene group entered is not a completely connected subgraph, try to analyze the genes separately or clique them first.")
   }else{  return(FriendsM)}

}

