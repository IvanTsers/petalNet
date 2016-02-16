
#' Graphing multiple vectors in one plot
#'
#' Plots multiple numeric vectors in rainbow color in not differently specified and saves the graph as a high resolution tiff
#'    
#' @param currentSubgroup names of identifiers to be plotted
#' @param expMatrix data matrix which includes identifiers in its rownames
#' @param fileName file name to which the plot to save to, no file extension needed
#' @param main overall title for the plot
#' @param cl specification for the default plotting color, default set to rainbow (each vector has its own color, note differentiation between colors might be quite small)
#'
#' @return None
#'
#' @export

plotExpProfiles <- function(currentSubgroup, expMatrix, fileName, main=NULL, cl=NULL){
 	index = match(currentSubgroup, rownames(expMatrix))
	n = length(index)
	if(is.null(cl)){farbe = rainbow(n)
		}else{  farbe = vector(mode="logical",length= n)
		        farbe[1:n] = cl
	}

	tiff(filename=paste(fileName,".tiff", sep=""), width=7, height=8.5, units="in", pointsize=12, compression="lzw", bg="white", res=600)
	plotProfiles(expMatrix[index,], colnames(expMatrix), main, farbe)
 	dev.off()
 }


