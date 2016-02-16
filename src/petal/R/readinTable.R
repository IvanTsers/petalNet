
#' Custom file read function for a table format
#'
#' Reads in the information from a file in table format, file must be tab delimited, first line must be table header
#'
#' @param fileName file name which to read in, character sting with .txt
#'
#' @return matrix with no rownames
#'
#' @export

readinTable <- function(fileName){
	m1 = scan(fileName, sep="\t", what="", nlines=1)
	MM = matrix(scan(fileName, sep="\t", what="", skip=1), ncol=length(m1), byrow=TRUE)
	colnames(MM) = m1
	return(MM)
}

