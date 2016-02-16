#' Calculating all pairwise assocations based on indicated metric
#'
#' This function takes a data matrix and calculates all pairwise associations of the rows within the matrix. Association measures are sorted from best to worst, the ordered measure matrix is returned and written to a .RData file 
#'
#' @param expMatrix numeric data matrix
#' @param metricChoice choice of metric, options include: Pearson Correlation "PE", Spearman Correlation "SP", Kendall Correlation "KE", Euclidean Distance "EU", Manhattan Distance "MA", and Canberra Distance "CA"
#'
#' @return ordered measure matrix, sorted from best to worst association measure
#'
#' @export


create_MMo <- function(expMatrix, metricChoice){

	load("numGenes.RData")
	MMoutputFile="MM.txt"

	# Calculating association matrix according to metric choice
	if(!is.na(match("PE", metricChoice))){
		corM = round(cor(t(expMatrix), method="pearson", use="pairwise.complete.obs"), digits=3)
	}

	if(!is.na(match("SP", metricChoice))){
		corM = round(cor(t(expMatrix), method="spearman", use="pairwise.complete.obs"), digits=3)
	}
	
	if(!is.na(match("KE", metricChoice))){
		corM = round(cor(t(expMatrix), method="kendal", use="pairwise.complete.obs"), digits=3)
	}

	if(!is.na(match("EU", metricChoice))){
		corM = round(as.matrix(dist(expMatrix,method="euclidean", diag=FALSE, upper=TRUE)), digits=2)
	}

	if(!is.na(match("MA", metricChoice))){
		corM = round(as.matrix(dist(expMatrix,method="manhattan", diag=FALSE, upper=TRUE)), digits=2)
	}

	if(!is.na(match("CA", metricChoice))){
		corM = round(as.matrix(dist(expMatrix,method="canberra", diag=FALSE, upper=TRUE)), digits=2)
	}

	# Writing correlation values to file: MM
	toFile = c("gene1", "gene2", metricChoice)
	writeToFile(t(toFile),MMoutputFile, FALSE)
	rm(toFile)
	tmp3 = NULL
	for(i in 1:(numGenes-1)){
		tmp1 = vector(mode="logical", length=numGenes-i)
		tmp1[1:length(tmp1)] = rownames(expMatrix)[i]
		tmp2 = rownames(expMatrix)[(i+1):numGenes]
		tmp3 = corM[i,(i+1):numGenes]
		v = cbind(tmp1, tmp2, tmp3)
		writeToFile(v,MMoutputFile)
	}
	rm(tmp1,tmp2, tmp3, v)

	# reading in MM and 
	MM = readinTable(MMoutputFile)
	if(isTRUE(metricChoice=="PE" || metricChoice=="SP" || metricChoice=="KE")){
		MMo = MM[order(as.numeric(MM[,3]),decreasing=TRUE),]	#ordering values in decreasing order, stongest correlation in first entry
	}else{  MMo = MM[order(as.numeric(MM[,3]),decreasing=FALSE),]}  #ordering values in increasing order, best distance in first entry

	save(MMo, file=paste("MMo_", metricChoice, ".RData", sep=""))
	rm(MM)
	system(paste("rm ",MMoutputFile, sep=""))

	return(MMo)
}

