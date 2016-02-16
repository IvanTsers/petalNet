#' Custom file read function for a numeric data matrix
#'
#' Reads in the infromation from a file, file must be tab delimited, first line must be table header, first column must be identifiers, inner part of the matrix must be numeric
#'
#' @param fileName file name which to read in, character sting with .txt
#'
#' @return readin numeric data matrix, rownames correspond to first column, colnames correspond to first row
#'
#' @export

readinExpM <- function(fileName){

	m1 = scan(fileName, sep="\t", what="", nlines=1)
	m2 = scan(fileName, sep="\t", what="", nlines=1, skip=1)

	if(isTRUE(length(m1)==length(m2))){
		expression     = matrix(scan(fileName,what="",sep="\t", skip=1),ncol=length(m1), byrow=TRUE)
		expM           = matrix(as.numeric(expression[,-1]), ncol = length(m1)-1, byrow=FALSE)
		colnames(expM) = m1[-1]
		rownames(expM) = expression[,1]
	}else{
		expression     = matrix(scan(fileName, sep="\t", what="", skip=1), ncol=length(m2), byrow=TRUE)
		expM           = matrix(as.numeric(expression[,-1]), ncol=length(m2)-1, byrow=FALSE)
		colnames(expM) = m1
		rownames(expM) = expression[,1]
	}
	rm(m1, m2,expression)

	# number of genes
	numGenes = dim(expM)[1]
	# saving number of genes  
	save(numGenes, file="numGenes.RData")
	save(expM, file="expM.RData")

	return(expM)
}

