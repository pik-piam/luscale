#' Spam2Mapping converter
#' 
#' Converts a spam object into a mapping and in addition returns statistics about cells
#' and clusters per region
#' 
#' @param spam spam object
#' @param cellregions a region vector with names cell of the spam object
#' @param clusterregions a region vector with names cluster of the spam object
#' @return data.frame containing the mapping and weights
#' @author Jan Philipp Dietrich
#' @export
#' @examples
#' spam <- spam::as.spam(matrix(c(0.2,0,0.8,0,0,1),2,3))
#' spam2mapping(spam,paste0("CE",1:3), paste0("CL",1:2))
#' @seealso \code{\link{clusterspam}}


spam2mapping <- function(spam, cellregions=NULL, clusterregions=NULL) {
  .cleanregions <- function(x,len) {
    if(length(x)==0) return(paste("GLO",1:len,sep="."))
    if(length(x)==1) return(paste(x,1:len,sep="."))
    if(length(x)!=len) stop("wrong number of regions supplied (!=",len,")!")
    x <- sub("\\..*$","",x)
    return(paste(x,1:len,sep="."))
  }

  out <- data.frame(cell=spam@colindices,cluster=NA,weight=spam@entries,row.names = NULL)
  for(c in 1:(length(spam@rowpointers)-1)) {
   out$cluster[spam@rowpointers[c]:min(length(out$cluster),spam@rowpointers[c+1])] <- c 
  }
  if(out$weight[1]==0 & dim(out)[1]>dim(spam)[2]) {
    out <- tail(out,-1)
  }
  out <- out[order(out$cell),]
  out$cell <- .cleanregions(cellregions,max(dim(spam)))[out$cell]
  out$cluster <- .cleanregions(clusterregions,min(dim(spam)))[out$cluster]
  rownames(out) <- NULL
  
  cluster <- table(sub("\\..*$","",unique(out$cluster)))
  cells <- table(sub("\\..*$","",out$cluster))
  stats <- rbind(cells,cluster)
  stats <- cbind(stats,Total=rowSums(stats))
  print(stats)
  attr(out,"stats") <- stats
  return(out)
}
