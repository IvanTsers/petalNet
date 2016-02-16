
#' Builds a small-world scale-free network 
#'
#' This function creates a number of network models to determine the best threshold producing a small-world scale-free network. All network models are accessible as .RData files. The file winningThresh.RData includes the possible thresholds to construct a sw-sf network with the first entry representing the 'best' threshold for the model construction.  
#'
#' @param expMFile data file name (.txt) in tab-separated format, first column must be vertex ID (genes), first row are condition identifiers
#' @param metric charcter defining association between each vertex ID, correlation options: "SP" - Spearman Correlation, "PE" - Pearson Correlation, and "KE" - Kendall as define in cor{stats}; distances: "EU" - Euclidean, "MA" - Manhattan, and "CA" - canaberra as defined in dist{stats}
#' @param thresholds numeric vector representing a series of thresholds that the user can choose; this variable can only be assigned with the metric is one of the three correlation metrics
#' @param MMoExists wheter the measure matrix exists (from a previous run), set to FALSE, if set to TRUE it will load the measure matrix (file must be in the same directory as the run)
#'
#' @return None, output files are generated
#'
#' @export

createSWSFnetFromFile <- function(expMFile=NULL, metric=c("SP", "PE", "KE", "EU", "MA", "CA"), thresholds=NULL, MMoExists=FALSE){
	options(warn=0)
	if(is.null(expMFile)) stop("Data matrix file must be provided")

	#checking for appropriate metric, if none entered Spearman is default
	metric = match.arg(metric)


	# validating thresholds if entered
	# thresholds can only be indicated for correlation measures, for all other metrics indicated thresholds are ignored
	if(metric=="SP" || metric=="PE" || metric=="KE"){

		if(is.null(thresholds)){
			CalculateThres = TRUE
		}else{
			if(min(thresholds)<0.7) stop("Threshold(s) must be greater or equal to 0.7")
			if(max(thresholds)>=1) stop("Threshold(s) must be smaller than 1")
			if(length(thresholds)>6) stop("The maximum number of thresholds is 5.")	
			CalculateThres = FALSE
		}
	}else{  
		print("A distance metric was chosen to define association, as a result the entered thresholds are ignored and will be calculated instead.")
		CalculateThres = TRUE
	}
	
	expM = readinExpM(expMFile)

	# if the number of rows is smaller than the number of columns 
	# this means there are less genes than measures, this is a problem, since the algorithm is written for large datasets
	# assuming the user uploaded the data in a transposed format
	if(dim(expM)[1]<dim(expM)[2]){
		expM = t(expM)
		warning("The input file has more columns than rows, as a result the dataset is transposed such that genes are in the rows and measures in the column.")
		warning("In case you are considered only a few genes over more measurements, this is not the program to use and we recommend to use another algorithm.")
	}

	if(MMoExists){
		fileToLoad = paste("MMo_", metric, ".RData", sep="")
		load(fileToLoad)
	}else{
		MMo        = create_MMo(expM, metric)
	}
	rm(expM)
	
	sortedV = as.numeric(MMo[,3])

	# checking if user specified thresholds or if they need to be calculated
	if(CalculateThres==FALSE){
		thresholds = c(sortedV[1], thresholds)
		thresholds = sort(thresholds, decreasing=TRUE)

		# making sure the thresholds are values within the MM
		for(i in 2:(length(thresholds))){
			tmp = which(sortedV>=thresholds[i])
			thresholds[i] = sortedV[max(tmp)]
		}
	}else{
		if(isTRUE(metric=="PE" || metric=="SP" || metric=="KE")){
			thresholds = generateThresholds(sortedV, 7)
		}else{thresholds = generateThresholds(sortedV, 7, FALSE)}
	}
	rm(sortedV)

	save(thresholds, file="thresholds.RData") 

	
	makeAdjMatrices(MMo, thresholds, metric)

	# cleaning directory by deleting unnecessary files
	system("rm AdjM*")
	system("rm Neigh*")
	system("rm ogenes.RData")
	system("rm edgeNum.RData")
	system("rm smMM*")
	rm(MMo)

	load("numGenes.RData")
	makeThresholdTable(metric, thresholds, numGenes)
	
	# cleaning directory by deleting unnecessary files
	system("rm Sym.*")

	# determining best threshold and alternatives
	winningThresh = findWinningThresh("NetworkStats.txt")
	save(winningThresh, file="winningThresh.RData")

	# .RData file made for downstreamAnalysis()
	winnerT     = winningThresh[1]
	association = metric
	save(winnerT, association, file = "forDA.RData")
	
	# writing general information about the data to file
	dataAssessToFile(expMFile, metric, CalculateThres, thresholds, winningThresh)
}
