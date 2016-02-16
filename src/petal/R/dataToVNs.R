
#' Builds a small-world scale-free network and extracts VNs based on genes of interest
#'
#' This function calls createSWSFnetFromFile() and downstreamAnalysis(). Its a one strep function to create a scale-free small-world network, determines the 'best' network model, from which then vicinity networks (VNs) are extracted based on vertex IDs the user can upload. For more specifics refer to the two functions createSWSFnetFromFile() and  downstreamAnalysis()    
#'
#' @param dataFile data file name (.txt) in tab-separated format, first column must be vertex ID (genes), first row are condition identifiers
#' @param metric charcter defining association between each vertex ID, correlation options: "SP" - Spearman Correlation, "PE" - Pearson Correlation, and "KE" - Kendall as define in cor{stats}; distances: "EU" - Euclidean, "MA" - Manhattan, and "CA" - canaberra as defined in dist{stats}
#' @param thresholds numeric vector representing a series of thresholds that the user can choose; this variable can only be assigned with the metric is one of the three correlation metrics
#' @param GoIFile a .txt file including vertex IDs of interest, IDs must be tab-separated
#' @param annoFile wheter there is an annoation file for the vertex IDs
#'
#' @return None, output files are generated
#'
#' @export

dataToVNs <- function(dataFile, GoIFile=NULL, annoFile=NULL, metric=NULL, thresholds=NULL){
	if(is.null(GoIFile)) stop("GoIFile must be specified, if no GoIFile please use createSWSFnetFromFile() instead.")

	createSWSFnetFromFile(dataFile, metric, thresholds)
	load("forDA.RData")
	if(is.na(winnerT)) stop("downstreamAnalysis() cannot be ran as no threshold was determined to construct a scale-free small-world network, please review NetworkStats.txt")
	downstreamAnalysis(winnerT, association, GoIFile, "myAnalysis.txt", dataFile, annoFile)
}
