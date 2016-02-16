
#' Generates a series of thresholds based on the distribution of the numeric vector 
#'
#' This function generates a series of thresholds based on the passed numeric vector. These thresholds are used to generate a number of visible network models.
#'
#' @param orderedV sorted numeric vector, for correlation and mutual infromation the vector should be in decending order, for distance the vector must be in increasing order (i.e. first entry should represent best association value)
#' @param adjMNum number indicating how many thresholds (i.e. network models) should be calculated 
#' @param decreasingOrder wheter the orderedV is in decending order (TRUE) or increasing order (FALSE)
#'
#' @return numeric vector
#'
#' @export

## AS IS THIS ONLY WORKS FOR LARGE NUMBER OF GENES NEED TO ADD SOMETHING TO INCLUDE SMALLER EXPERIMENTS
## CUTOFF DOES NOT WORK WITH n<301 AND IN GENERAL DOES NOT MAKE MUCH SENSE FOR ANY n<500 MAYBE EVEN MORE
## GO BACK TO PREVIOUS IDEA USING TOP 5%/10%

# orderedV should be decending for Pearson, Spearman, MI
# orderedV should be increasing for pCos, Euclidean

generateThresholds <-function(orderedV,adjMNum, decreasingOrder=TRUE){

	load("numGenes.RData")	

	cutoff = orderedV[150*numGenes]
	firstT = orderedV[0.5*numGenes]

	delta = abs(firstT-cutoff)/(adjMNum-1)

	if(decreasingOrder){
		tvalue = c(orderedV[1], firstT, firstT-delta)
		for(i in 2:(adjMNum-2)){tvalue = c(tvalue, firstT-i*delta)}
		tvalue = c(tvalue, cutoff)

		for(i in 2:(length(tvalue)-1)){
			tmp = which(orderedV>=tvalue[i])
			tvalue[i] = orderedV[max(tmp)]
		}
	}else{
		tvalue = c(orderedV[1], firstT, firstT+delta)
		for(i in 2:(adjMNum-2)){tvalue = c(tvalue, firstT+i*delta)}
		tvalue = c(tvalue, cutoff)

		for(i in 2:(length(tvalue)-1)){
			tmp = which(orderedV<=tvalue[i])
			tvalue[i] = orderedV[max(tmp)]
		}
	}
	return(tvalue)
}

