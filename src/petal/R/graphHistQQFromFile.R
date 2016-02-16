
#' Creates a tiff image of data histogram and the corresponding Q-Q plot from a file
#'
#' Creates a tiff files including the histogram of the data as well as the corresponding q-qplot see qqnorm{stats} for more information. Input file is read in and the data matrix is stored as expM.RData. 
#'
#' @param inputFile file name of a data matrix; file must be tab delimited, first line must be table header, first column must be identifiers, inner part of the matrix must be numeric
#' @param fileName prefix of the desired output file, defaul is "Data_Hist-QQ"
#'
#' @return None; output file is generated
#'
#' @export

graphHistQQFromFile <- function(inputFile, fileName="Data_Hist-QQ"){

	m1 = scan(inputFile, sep="\t", what="", nlines=1)
	m2 = scan(inputFile, sep="\t", what="", nlines=1, skip=1)

	if(isTRUE(length(m1)==length(m2))){
		expression     = matrix(scan(inputFile,what="",sep="\t", skip=1),ncol=length(m1), byrow=TRUE)
		expM           = matrix(as.numeric(expression[,-1]), ncol = length(m1)-1, byrow=FALSE)
		colnames(expM) = m1[-1]
		rownames(expM) = expression[,1]
	}else{
		expression     = matrix(scan(inputFile, sep="\t", what="", skip=1), ncol=length(m2), byrow=TRUE)
		expM           = matrix(as.numeric(expression[,-1]), ncol=length(m2)-1, byrow=FALSE)
		colnames(expM) = m1
		rownames(expM) = expression[,1]
	}
	rm(m1, m2,expression)

	save(expM, file="expM.RData")
	
	measures = as.vector(expM)
	minExp   = round(min(measures, na.rm=TRUE)-0.5,digits=0)
	maxExp   = round(max(measures, na.rm=TRUE)+0.5, digits=0)
	delta    = (maxExp -minExp)/100

	tiff(filename=paste(fileName,".tiff", sep=""), width=7, height=8.5, units="in", pointsize=12, compression="lzw", bg="white", res=600)
	par(mfrow = c(2,1))
	hist(measures, breaks=seq(minExp,maxExp , by=delta),main="Histogram of Data", col="blue")
	qqnorm(expM)
	qqline(expM, col="red")
	dev.off()
}	


