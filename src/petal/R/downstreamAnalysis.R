
#' Automated subnetwork extraction
#'
#' This function finds all maximal cliques of a subgroup of vertex ids, then the vicinity network (1-neighbor network) is identified for each maximal clique.  
#'
#' @param winnerT threshold used to construct the network model
#' @param metric metrix used to construct the network model
#' @param GoIFile gene of interest (GoI) file, a txt file inlucuding vertex ids which are of particular interest
#' @param outFile indicate a filename to which to write the information to
#' @param expMatrixFile filename of original datafile (expression matrix)
#' @param annoTFile annotation file if there is one, first column must correspond to vertex IDs, first row must be column headers
#' @param adjMatrix custom adjacency matrix can be loaded if petal network construction was not used
#'
#' @return None, but multiple files are created 
#'
#' @export

downstreamAnalysis <-function(winnerT, metric, GoIFile, outFile, expMatrixFile, annoTFile=NULL, adjMatrix=NULL){

# if no adjMatrix passed to the function it will load it from within the directory
	if(is.null(adjMatrix)){
		adjF    = paste("adjM", metric, "_", winnerT, ".RData", sep="")
		load(adjF)
	}else{
		adjM = adjMatrix
		rm(adjMatrix)
	}

#  reading in annotation file
	if(is.null(annoTFile)){
		noAnnoT = TRUE
	}else{
		annoT = readinTable(annoTFile)
		noAnnoT = FALSE
	}

	expM = readinExpM(expMatrixFile)

# reading in node stats table
	nodeStats = getVertexStats(adjM, TRUE)
									  
# load genes and checking if GoIs are duplicated, removing duplicates
	genes = scan(GoIFile, sep="\t", what="")
	genes = unique(genes) 

# removing possible NAs within GoI list
	removeIndex = as.vector(which(is.na(genes)))
	if(length(removeIndex)!=0){genes = genes[-removeIndex]}
	rm(removeIndex)

# to file: the number of unique gene identifiers identified from the input file
	toFile = paste(length(genes), " unique gene identifiers are determined from the input file (",GoIFile ,").", sep="")
	writeToFile(toFile, outFile, FALSE)
	rm(toFile)

# checking if GoIs are in the network, if not they are removed
	tmpIndex = match(genes, rownames(adjM))
	tmpX     = as.vector(which(is.na(tmpIndex)))

	if(length(tmpX)!=0){
		if(length(tmpX)==1){
			toFile = paste("Gene ", genes[tmpX], " is not in the ", metric, "-", winnerT, " network model.", sep="")
			writeToFile(toFile, outFile, TRUE)
			rm(toFile)
		}else{
			toFile = paste(length(tmpX), " genes are not in the ",  metric, "-", winnerT, " network model, see file 'GenesNotInNet.txt' for the list of genes.", sep="")
			writeToFile(toFile, outFile, TRUE)
			writeToFile(genes[tmpX], "GenesNotInNet.txt", FALSE)
			rm(toFile)
		}
		genes = genes[-tmpX]
		cat("\n", file=outFile, append=TRUE)
		toFile = paste(length(genes), " genes remain for furhter downstream analysis.", sep="")
		writeToFile(toFile, outFile, TRUE)
		cat("\n", file=outFile, append=TRUE)
		rm(toFile)
	}

	geneIndexAdjM = match(genes, rownames(adjM))
	g_tmp = NULL

	if(length(genes)==1){
		writeToFile("One gene of interest was entered.", outFile, TRUE)
		neighborM = makeVNM(genes, adjM)
		toFile = paste(genes, " has ", dim(neighborM)[1]-1, " neighbors:",sep="")
		writeToFile(toFile, outFile, TRUE)
		rm(toFile)

		poop = match(genes,rownames(neighborM))
		neighbors = rownames(neighborM)[-poop]

		writeToFile(t(neighbors), outFile, TRUE)
		densityM(neighborM, outFile,FALSE)
	}else{
		nodeIndex = match(genes, rownames(adjM))
		GoIM      = adjM[nodeIndex,nodeIndex]
		MCs       = findMaximalCliques(GoIM)
		writeToFile("The genes being investigated are: ", outFile, TRUE)

# sorting node stats table and adding colnames into the datafram to avoid warning messeges 
		tmpx = colnames(nodeStats)
		tmpNodeStats = nodeStats[geneIndexAdjM,]
		orderNStats = tmpNodeStats[order(tmpNodeStats[,1]),]
		orderNStats = rbind(tmpx,orderNStats)
		write.table(orderNStats, outFile, TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
		if(length(MCs)==1){
			toFile = paste("Within their sub-network there is ", length(MCs), " maximal cliques", sep="")
		}else{
			toFile = paste("Within their sub-network there are ", length(MCs), " maximal cliques", sep="")
		}
		writeToFile(toFile, outFile, TRUE)
		rm(toFile)
		cat("\n", file=outFile, append=TRUE)
		cat("\n", file=outFile, append=TRUE)

		for(i in 1:length(MCs)){
			currentFile = paste("VN", i, ".txt", sep="")
			group  = MCs[[i]]
			gNames = rownames(GoIM)[group]

			toFile = paste("Vicinity Network ",i, " is based on:", sep="")
		# writing information to general output file
			writeToFile(toFile, outFile, TRUE)
			writeToFile(t(gNames), outFile, TRUE)
			cat("\n", file=outFile, append=TRUE)

		# writing information to specific vicinity network file
			writeToFile(toFile, currentFile, FALSE)
			writeToFile(t(gNames), currentFile, TRUE)
			cat("\n", file=currentFile, append=TRUE)

			neighborM = makeVNM(gNames, adjM)
			toFile = paste(dim(neighborM)[1]-length(gNames), " neighbors", sep="")
			writeToFile(toFile, currentFile, TRUE)
			neighDense = densityM(neighborM, currentFile)

		# if an annotation file was provided, it is writing to file, otherwise just the genes are written to file
			if(noAnnoT){
				writeToFile("Genes:", currentFile, TRUE)
				writeToFile(rownames(neighborM),currentFile,TRUE)
			}else{
				cat("\n", file=currentFile, append=TRUE)
				writeToFile("Annotation:", currentFile, TRUE)
				indexA = match(rownames(neighborM), annoT[,1])
				writeToFile(t(colnames(annoT)),currentFile,TRUE)
				writeToFile(annoT[indexA,],currentFile,TRUE)
			}

			cat("\n", file=currentFile, append=TRUE)
			cat("\n", file=currentFile, append=TRUE)

			filePrefix = paste("VN", i, "_pofiles", sep="")
			plotExpProfiles(rownames(neighborM), expM, filePrefix)

			l_tmp = c(i, dim(neighborM)[1], length(gNames), round(neighDense, digits=2))
			g_tmp = rbind(g_tmp, l_tmp)
		}

		writeToFile("Vicinity Network Summary:", outFile, TRUE)
		t_header = c("VN", "VNsize", "GoInum", "density")
		writeToFile(t(t_header), outFile, TRUE)
		writeToFile(g_tmp, outFile, TRUE)
	}
}