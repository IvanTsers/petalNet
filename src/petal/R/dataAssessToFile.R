
#' Creates an output file providing general information about the data file and its analysis
#'
#' This function creates an output file including multiple infromation about the data at hand as well as steps taken in the analysis. This is an internal function of createSWSFnet
#'
#' @param expMFile filename containing data matrix 
#' @param metricUsed metric being used to define association
#' @param CalculateThres whether user specified thresholds or not, boolean: TRUE or FALSE
#' @param thresholds list of thresholds considered for the network model generation
#' @param winningThresh 'best' threshold that generates a scale-free, small-world network
#'
#' @return None; output file is generated
#'
#' @export

dataAssessToFile <-function(expMFile, metricUsed, CalculateThres, thresholds, winningThresh){
	
	outFile = "GeneralInfo.txt"
	load("expM.RData")

	# checking how many values are missing within the expression matrix
	x = which(is.na(expM))
	# checking how many zero entries within the expression matrix
	y = which(expM==0)
	
# writing general information to file
	toFile = paste("General information provided about the analyis of : ", expMFile, sep="")
	writeToFile(toFile, outFile, FALSE)
	rm(toFile)

# dimension of data matrix
	toFile = paste("Dimension: ", dim(expM)[1], " (number of genes/rows) across ", dim(expM)[2], " (measures/columns)", sep="")
	writeToFile(toFile, outFile)
	rm(toFile)

# info about missing values
	if(isTRUE(length(x)==0)){
		toFile = paste("There are no missing values.")
	}else{  toFile = paste("There are ", length(x), " missing values out of ", length(expM), " possible expression measures (", round((length(x)/length(expM))*100,digits=3), "% of the entries are missing within the data matrix)", sep="")}
	writeToFile(toFile, outFile)
	rm(toFile)

# info about zero entries
	if(isTRUE(length(y)!=0)){
		toFile = paste("There are ", length(y), " entries with a value of zero.", sep="")
		writeToFile(toFile, outFile)
		rm(toFile)
	}

# info about range of data
	toFile = paste("Expression measures range between ", round(min(expM, na.rm=TRUE), digits=3), " and ", round(max(expM, na.rm=TRUE), digits=3), sep="")
	writeToFile(toFile, outFile, TRUE)
	rm(toFile)
	cat("\n", file=outFile, append=TRUE)

# specifying which metric is used
	if(metricUsed=="SP"){toFile = "Spearman is used to calculate the pairwise correaltion between all gene pairs."}
	if(metricUsed=="PE"){toFile = "Pearson is used to calculate the pairwise correaltion between all gene pairs."}
	if(metricUsed=="KE"){toFile = "Kendall is used to calculate the pairwise correaltion between all gene pairs."}
	if(metricUsed=="EU"){toFile = "Euclidean is used to calculate the pairwise distance between all gene pairs."}
	if(metricUsed=="MA"){toFile = "Manhattan is used to calculate the pairwise distance between all gene pairs."}
	if(metricUsed=="CA"){toFile = "Canberra is used to calculate the pairwise distance between all gene pairs."}
	writeToFile(toFile, outFile)
	rm(toFile)
		
# threshold information
	if(CalculateThres){writeToFile("No thresholds were entered, hence the algorithm will choose a threshold which produces a scale-free and small-world network", outFile)}

	writeToFile("Considered thresholds are:", outFile)
	writeToFile(t(thresholds), outFile)
	cat("\n", file=outFile, append=TRUE)

	writeToFile("The winning threshold is:", outFile)
	writeToFile(winningThresh[1], outFile)
	cat("\n", file=outFile, append=TRUE)
}

