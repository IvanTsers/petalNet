
#' Internal function of plotExpProfiles
#'
#' Creates canvas and scales the plot for plotExpProfiles().
#'    
#' @param expressionM datamatrix which to plot all indicies 
#' @param cond colnames of the data matrix
#' @param header overall title for the plot
#' @param farbe specification for the default plotting color
#'
#' @return None

plotProfiles <- function(expressionM, cond, header, farbe){
  plot.new()						
  par(ps=12,mar=c(4,4,2,2),font=2)			
 

  v = as.vector(expressionM)
  MAX<-max(v, na.rm=TRUE)
  MIN<-min(v, na.rm=TRUE)
 
  X<-c(1:length(cond))

  plot(c(1,X[length(X)]),c(MIN-1,MAX),pch=16,cex=0, font.lab=1,cex.lab=1.2,cex.axis=1.2,main=header, xlab="",ylab="Expression",xaxt="n")
  axis(1, at = 1:length(cond), labels =substr(colnames(expressionM), 1,6))
  if(length(expressionM) > length(cond)){
        for(i in 1:dim(expressionM)[1]){
              Y<-expressionM[i,]
              points(X,Y,pch=16,cex=.5,col=farbe[i])
              lines(X,Y,lwd=1,col=farbe[i])
        }
  }
}

