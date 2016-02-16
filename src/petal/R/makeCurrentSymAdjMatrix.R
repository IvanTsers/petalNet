#' Internal function of makeAdjMatrices
#'
#' Adjacency matricies are produced and checked for correctness
#'
#' @param shm metric prefix
#' @param currentThreshIndex threshold index that is currently considered
#' @param thresholds vector of numeric values
#' @param edgeNum numer of edges for the considered threshold
#'
#' @return edgeNum, the number of edges for the network model under considered threshold

makeCurrentSymAdjMatrix <- function(shm,currentThreshIndex,thresholds,edgeNum){
	load("ogenes.RData")

	currentThresh <- thresholds[currentThreshIndex]
	prevThresh    <- thresholds[(currentThreshIndex -1)]
	print(currentThresh)
	cur.suffix <- paste(shm, "_",prevThresh, "to", currentThresh,".RData", sep="")
	cur.suffix.txt <- paste(shm, "_",prevThresh, "to", currentThresh,".txt", sep="")
	currentMMfile = paste("smMM.", cur.suffix, sep="")

	load(currentMMfile)

	edgeNum[currentThreshIndex] <- dim(currentMM)[1]

	currentadjfile = paste("AdjM", cur.suffix.txt, sep="")								# first line:list of genes (ogenes)
	currentneighbrfile = paste("Neighbors.", cur.suffix.txt, sep="")					# first column gene names (ogenes)

	write.table(NULL, file=currentneighbrfile, sep="\t", append=FALSE,row.names=FALSE, col.names=FALSE)
	write.table(t(ogenes), file=currentadjfile, sep="\t", append=FALSE,row.names=FALSE, col.names=FALSE)

	for(i in 1:length(ogenes)){
		neighbrvec <- rep(0, length(ogenes))
		ind1 = as.vector(which(ogenes[i]==currentMM[,1]))								# the first gene shows up in the MM here

		if(length(ind1)==0 && currentThreshIndex==2){
			write.table(ogenes[i],  file=currentneighbrfile, sep="\t", append=TRUE,row.names=FALSE, col.names=FALSE)
			write.table(t(neighbrvec),file=currentadjfile, sep="\t", append=TRUE,row.names=FALSE, col.names=FALSE)
		}else{
			neighbrs = match(currentMM[ind1,2],ogenes)									# the neighbors of this gene, at this threshold
			neighbrvec[neighbrs] <- 1
			if(currentThreshIndex == 2){ 												# add the current gene to the neighbors as a row name
				neighbrs <- c(ogenes[i],neighbrs[order(neighbrs)])
			}
			if(currentThreshIndex > 2){													# need to include the thresholds of earlier thresholds
				prev.suffix <- paste(shm, "_", thresholds[(currentThreshIndex-2)],"to",prevThresh,".txt", sep="")
				prev.neighbrfile = paste("Neighbors.", prev.suffix, sep="")
				prev.neighbrs <- as.numeric(scan(prev.neighbrfile, what="", nlines=1, skip = i-1, quiet=TRUE)[-1])
				neighbrvec[prev.neighbrs] <- 1
				tmp <- c(neighbrs, prev.neighbrs)
				neighbrs <- c(ogenes[i],tmp[order(tmp)])
			}

			write.table(t(neighbrs), file=currentneighbrfile, sep="\t", append=TRUE,row.names=FALSE, col.names=FALSE)
			write.table(t(neighbrvec),file=currentadjfile, sep="\t", append=TRUE,row.names=FALSE, col.names=FALSE)
		}
	}

	# CHECK THE NUMBER OF EDGES IN THE ADJ MATRIX BEFORE MAKING IT SYMMETRIC AND ADDING DIAGONAL
	# CURRENT ADJ MATRIX SHOULD THE NUMBER OF EDGES = THE SUM OF DIM smMM's FROM THRESH INDEX = 2 TO CURRENT THRESH INDEX
	currentAdjM <- matrix(scan(currentadjfile, what=0,sep="\t", skip=1, quiet=TRUE),ncol=length(ogenes),byrow=T)
	#print(edgeNum)
	if(sum(currentAdjM) != sum(edgeNum[1:currentThreshIndex])) stop(paste("problem at ", currentThresh))

	## MAKE THE MATRIX SYMMETRIC, ADD DIAGONAL, WRITE IT BACK OUT
	colnames(currentAdjM) = rownames(currentAdjM) = ogenes
	SymM = currentAdjM + t(currentAdjM)
	diag(SymM) <- 1
	write.table(SymM, file = paste("Sym.",shm, ".",currentadjfile, sep=""), sep="\t",row.names=FALSE, col.names=TRUE)
	save(edgeNum, file="edgeNum.RData")
	return(edgeNum)
}