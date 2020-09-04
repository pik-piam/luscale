#' Spam2Mapping converter
#' 
#' Converts a spam object into a mapping and in addition returns statistics about cells
#' and clusters per region
#' 
#' @param spam spam object
#' @param regions a region vector with names for each column of the spam object
#' @return data.frame containing the mapping and weights
#' @author Jan Philipp Dietrich
#' @export
#' @examples
#' spam <- spam::as.spam(matrix(c(0.2,0,0.8,0,0,1),2,3))
#' spam2mapping(spam,c("A","A","B"))
#' @seealso \code{\link{clusterspam}}


spam2mapping <- function(spam,regions="GLO") {
  regions <- sub("\\..*$","",regions)
  out <- data.frame(cell=spam@colindices,cluster=NA,weight=spam@entries,row.names = NULL)
  for(c in 1:(length(spam@rowpointers)-1)) {
   out$cluster[spam@rowpointers[c]:min(length(out$cluster),spam@rowpointers[c+1])] <- c 
  }
  if(out$weight[1]==0 & dim(out)[1]>dim(spam)[2]) {
    out <- tail(out,-1)
  }
  out <- out[order(out$cell),]
  out$cell <- paste(regions,out$cell,sep=".")
  out$cluster <- paste(regions,out$cluster,sep=".")
  rownames(out) <- NULL
  
  cluster <- table(sub("\\..*$","",unique(out$cluster)))
  cells <- table(sub("\\..*$","",out$cluster))
  stats <- rbind(cells,cluster)
  stats <- cbind(stats,Total=rowSums(stats))
  print(stats)
  attr(out,"stats") <- stats
  return(out)
}
