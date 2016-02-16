
#' Internal function of makeAdjMatrices
#'
#' Splits the large ordered measure matrix into smaller matrices according to the indicated thresholds 
#'    
#' @param orderedMM measure matrix including all pairwise comparisions between vertex IDs based on a particular metric defining assosiction; association measured are sorted from best to worst
#' @param thresholds vector of numeric values indicating the thresholds on which to build the network models
#' @param metric charater (i.e. "SP", "PE", "EU", etc) indicating the metric used to define association between vertex IDs
#'
#' @return None

create_smMMs <- function(orderedMM, thresholds, metric){
	startIndex = 1
        for(i in 2:length(thresholds)){

                x = which(thresholds[i]==orderedMM[,3])
                cutIndex = max(x)
			 
		prefix = paste("smMM.", metric, "_", thresholds[i-1], "to",thresholds[i], sep="")
		
		filetxt   = paste( prefix, ".txt", sep="")
		fileRData = paste( prefix, ".RData", sep="")

                currentMM = orderedMM[startIndex:cutIndex,]

                write.table(currentMM,file=filetxt , sep="\t", quote=F, row.names=F, col.names=T)
                save(currentMM, file = fileRData)
                startIndex = cutIndex+1
        }
        ogenes = sort(unique(c(orderedMM[1:cutIndex, 1], orderedMM[1:cutIndex, 2])))
        save(ogenes, file="ogenes.RData")
}
	


