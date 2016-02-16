
#' Function creating adjacency matrices for each threshold
#'
#' For each threshold an binary adjacency matrix is constructed.
#'
#' @param orderedMMo measure matrix including all pairwise comparisions between vertex IDs based on a particular metric defining assosiction; association measured are sorted from best to worst
#' @param thresholds thresholds defining similarity, list of numeric values
#' @param metric charater (i.e. "SP", "PE", "EU", etc) indicating the metric used to define association between vertex IDs
#'
#' @return None, output files are generated
#'
#' @export

makeAdjMatrices <-function(orderedMMo, thresholds, metric){
	create_smMMs(orderedMMo, thresholds, metric)
	numEdges <- rep(0,length(thresholds))
	for(currentTIndex in 2:length(thresholds)){
		numEdges = makeCurrentSymAdjMatrix(metric,currentTIndex, thresholds,numEdges)	#KS function
	}
}
