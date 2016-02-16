
#' Creates a tiff image of data histogram and the corresponding Q-Q plot
#'
#' Creates a tiff files including the histogram of the data as well as the corresponding q-qplot see qqnorm{stats} for more information based on a numeric data matrix
#'
#' @param dataMatrix numeric data matrix
#' @param fileName prefix of the desired output file, defaul is "Data_Hist-QQ"
#'
#' @return None; output file is generated
#'
#' @export

graphHistQQ <- function(dataMatrix, fileName="Data_Hist-QQ"){
	measures = as.vector(as.numeric(dataMatrix))
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


