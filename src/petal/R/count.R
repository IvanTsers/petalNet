
#' counting function
#'
#' This function counts the appearance of a particular variable within each row of a data matrix
#'
#' @param variable character or numeric, variable that should be counted
#' @param M data matrix in which to find the variable in
#'
#' @return numeric vector indicating the number of appearances of the variable of interest for each row of the provided data matrix
#'
#' @export

count <-function(variable, M){
	c = vector(mode="numeric", length = dim(M)[1])
	for(i in 1:dim(M)[1]){
		c[i] = length(which(M[i,]==variable))
	}
	return(c)
}

