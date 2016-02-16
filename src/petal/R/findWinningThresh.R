
#' Functions calculates best and alternative thresholds for a scale-free small-world network
#'
#' This function takes the network statistics file created by makeThresholdTable() and weights the parameters against each other to determine the best network that is scale-free and small-world and follows other important properties such as the consideration of the biggest component of a network model. 
#'
#' @param threshFile filename of file containing network statistics
#' @param minR2 maximum acceptable slope of the linear regression of the log transformed degree distribution; default set to 0.8
#' @param maxSlope maximum acceptable slope of the linear regression of the log transformed degree distribution; default set to -1; slope should be -3, -1
#' @param minUsed minimum acceptable percentage of the original dataset used to build the network model; default is 40
#' @param minBigComp minimum acceptable size of the biggest component in percentage; default set to 95 
#'
#' @return possibleThresh numeric vector presenting possible thresholds which construct scale-free small-world network models, the vector is sorted from best to worse 
#'
#' @export

findWinningThresh <-function(threshFile, minR2=0.80, maxSlope=-1, minUsed=40, minBigComp=95){

	m1 <- scan(threshFile,what="",sep="",nlines=1)
	topoT <- matrix(as.numeric(scan(threshFile,what="",sep="", skip=1)),ncol=length(m1),byrow=T)
	colnames(topoT) = m1

	ThresholdTable = topoT[order(as.numeric(topoT[,2]), decreasing=TRUE),]	

	# 1 - threshold
	# 2 - R^2
	# 3 - slope/power
	# 4 - meanCC
	# 5 - meanPath
	# 6 - %used
	# 7 - %bigComp

	r2 = which(ThresholdTable[,2]<minR2)
	r3 = which(ThresholdTable[,3]>maxSlope)
	r6 = which(ThresholdTable[,6]<minUsed)
	r7 = which(ThresholdTable[,7]<minBigComp)
	r = unique(c(r2,r3,r6,r7))
	if(length(r)>0){RankTT = ThresholdTable[-r,1:7]
	}else{RankTT = ThresholdTable[,1:7]}


	if(is.null(dim(RankTT))){ 				# this is true if there is only one threshold left
		possibleThresh = as.vector(RankTT[1])
	}else if(dim(RankTT)[1]==0){
		warning("None of the considered thresholds produce a scale-free and small-world network")
		possibleThresh = NA
	}else{ 
		RankTT[,2] = dim(RankTT)[1]-rank(round(RankTT[,2], digits=2), ties.method="max")+1
		RankTT[,3] = rank(round(RankTT[,3], digits=2),ties.method="min")
		RankTT[,4] = dim(RankTT)[1]-rank(round(RankTT[,4], digits=2), ties.method="max")+1
		RankTT[,5] = rank(round(RankTT[,5], digits=1), ties.method="min")
		RankTT[,6] = dim(RankTT)[1]-rank(round(RankTT[,6], digits=0), ties.method="max")+1
		RankTT[,7] = dim(RankTT)[1]-rank(round(RankTT[,7],digits=0), ties.method="max")+1

		numBest   = count(1, RankTT)
		SumRankTT = cbind(RankTT[,1],rowSums(RankTT[,2:7]),numBest)
		colnames(SumRankTT) = c("threshold", "rankSum", "numBest")
		rownames(SumRankTT) = seq(1, dim(SumRankTT)[1], by = 1)
		print(SumRankTT)
		SumRankTTo = SumRankTT[order(as.numeric(SumRankTT[,2]), decreasing=FALSE),]
		possibleThresh = as.numeric(SumRankTTo[,1])
	}

	return(possibleThresh)

}

