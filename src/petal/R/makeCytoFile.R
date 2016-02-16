
#' Creates an network file for Cytoscape
#'
#' Creates a txt file which can be loaded as is into Cytoscape to plot the network based on the indicated threshold and metric.
#'
#' @param threshold numeric value indicating the threshold building a network model
#' @param metric charater (i.e. "SP", "PE", "EU", etc) indicating the metric used to define association between vertex IDs
#' @param orderedMM measure matrix including all pairwise comparisions between vertex IDs based on a particular metric defining assosiction; association measured are sorted from best to worst
#'
#' @return None; output file is generated
#'
#' @export

makeCytoFile <-function(threshold, metric, orderedMM){
	if(colnames(orderedMM)[3]!=metric) stop(paste("The provided measure matrix is based on ", colnames(orderedMM)[3], sep=""))

	x = which(threshold==orderedMM[,3])
	cutIndex = max(x)
	fileName = paste("CytoNet","_", metric, threshold,  ".txt", sep="")
	write.table(orderedMM[1:cutIndex,],file=fileName , sep="\t", quote=F, row.names=F, col.names=T)
}


