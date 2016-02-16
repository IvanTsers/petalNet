
#' Calculatiing a varity of network statistics for a number of network models
#'
#' This function calls all adjancency matrices which were build by makeAdjMatrices(), removed unconnected vertices from each network model and calculates network statistics for each model. Network statistics include: "R^2" and "slope/power" of the log-transformed connectrivity distribution to determine scale-freeness, "meanCC" and "meanPath" standing for the mean cluster coefficient and mean pathlength within the model to classify the network as small-world, "%used" indicating the percentage vertex IDs used from the original dataset, "%bigComp" the percentage of vertices included in the biggest component,"numComp" the number of components, "kmax" the maximal number of connections, "kmean" the average connectity of all verticies in the network, "kmedian" the median connectity of all verticies in the network, "diameter" of the network which is equal to the shortest longest path within the network model. 
#'
#' @param metricChoice chracter indicating metric choice, must be one of the following: "PE", "SP", "KE", "EU", "MA", "CA"
#' @param thresholds threholds considered to create the different network models
#' @param numGenes number of rows in original data matrix
#' @param outputFile file name set to "NetworkStats.txt", but can be adjusted
#'
#' @return none, output file is generated
#'
#' @export

makeThresholdTable <-function(metricChoice, thresholds, numGenes, outputFile = "NetworkStats.txt"){
	
	outputFile = match.arg(outputFile)
	library("igraph")

	toFile=c("threshold","R^2","slope/power", "meanCC", "meanPath","%used", "%bigComp","numComp", "kmax", "kmean", "kmedian", "diameter")
	writeToFile(t(toFile),outputFile, FALSE)
	rm(toFile)
	
	for(currentThreshIndex in 2:(length(thresholds))){
		currentThresh = thresholds[currentThreshIndex]
		prevThresh    = thresholds[(currentThreshIndex -1)]

		txtSuffix         = paste(metricChoice, "_",prevThresh, "to", currentThresh,".txt", sep="")
		currentSymAdjfile = paste("Sym.",metricChoice, ".AdjM",txtSuffix, sep="")

		x = scan(currentSymAdjfile, what="", nlines=1, sep="\t")
		adjMatrix = matrix(scan(currentSymAdjfile, what=0, skip=1, sep="\t"), ncol=length(x), byrow=TRUE)
		colnames(adjMatrix) = rownames(adjMatrix) = x

		tmpV        = rowSums(adjMatrix)
		removeIndex = as.vector(which(tmpV==1))
		if(isTRUE(length(removeIndex)==0)){ adjM = adjMatrix
		}else{adjM = adjMatrix[-removeIndex,-removeIndex]}
		rm(adjMatrix, removeIndex,tmpV, x)

		save(adjM, file=paste("adjM", metricChoice, "_", currentThresh, ".RData", sep=""))

		#converting adjacency matrix to a graph format for igraph to use it
		GraphAdjM = graph.adjacency(adjM, mode="undirected", weighted=NULL, diag=FALSE, add.colnames=NULL)

		#diameter of network, longest shortes path
		Dmeter      = diameter(GraphAdjM, directed = FALSE)
		PathAverage = average.path.length(GraphAdjM, directed=FALSE)
		AverageCC   = transitivity(GraphAdjM, type="average")
		
		GraphComponents = clusters(GraphAdjM)
		numComponents   = GraphComponents$no
		PercentBigComp  = (max(GraphComponents$csize)/dim(adjM)[1])*100

		connectivity   = rev(sort(rowSums(adjM)-1)) 		# number of connections per row for adjM with diag=1
		PercentageUsed = (dim(adjM)[1]/numGenes)*100
		ConnectMedian  = median(connectivity)
		ConnectMean    = mean(connectivity)
		ConnectMax     = max(connectivity)

		h<-hist(connectivity, breaks=seq(min(connectivity)-1, max(connectivity)+1, by=1), plot=FALSE)

		#Removing nodes with no connectifity since log(0)=undf
		k  = h$breaks[-1]
		pk = h$counts

		#log-scale
		#plot(log2(k),log2(pk))

		#remove unexisting connectivity degrees
		if(!is.na(match(-Inf,log2(pk)))){
			removeIndex=grep(-Inf,log2(pk))
			newk  = log2(k)[-removeIndex]
			newpk = log2(pk)[-removeIndex]
		}else{  newk  = log2(k)
			newpk = log2(pk)
		}

		linearFit = lm(newpk~newk)
		#par(mfrow=c(2,2))
		#plot(linearFit)
		linearFitSummary = summary(linearFit)
		Intercept        = linearFitSummary$coefficients[1]
		slope            = linearFitSummary$coefficients[2]
		R2               = linearFitSummary$r.squared

		tmpT = c(currentThresh, R2, slope, AverageCC, PathAverage, PercentageUsed, PercentBigComp, numComponents, ConnectMax, ConnectMean, ConnectMedian,  Dmeter)
		writeToFile(round(t(tmpT), digits=4),outputFile)
	}
}
