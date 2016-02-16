
#' Custom writing to file function 
#'
#' This is a custom function to write information to file without row or column names
#'
#' @param information material to be written to file, matrix or vector
#' @param fileName file name which to write to, charcter sting with .txt
#' @param appending whether to append to the specified file or overwrite, boolean variable TRUE(default) or FALSE
#'
#' @return None
#'
#' @export

writeToFile <- function(information,fileName, appending=TRUE){
	write.table(information,file=fileName, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=appending)
}



# @examples
# Appends the vector to myFile.txt
# Each numeric value is in a separate line
# x = c(0,1,2,3,4)	
# writeToFile(x, "myFile.txt", TRUE)

# Creates a new myFile.txt
# Vector x is writting to file within one line, indices are separated by tab
# writeToFile(t(x), "myFile.txt", FALSE)