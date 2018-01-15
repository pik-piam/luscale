#' MAgPIE Kmeans clustering
#' 
#' Performs MAgPIE kmeans clustering and calculates corresponding spam relation
#' matrix
#' 
#' 
#' @usage mag_kmeans(cdata,ncluster)
#' @param cdata a cluster data file as produced by cluster_base
#' @param ncluster The desired total number of clusters.
#' @return A spam relation matrix
#' @author Jan Philipp Dietrich
#' @importFrom stats kmeans
#' @seealso \code{\link{cluster_per_region}}, \code{\link{mag_hierarchical}},
#' \code{\link{clusterspam}}
mag_kmeans <- function(cdata,ncluster) {
  cpr <- cluster_per_region(cdata,ncluster)
  spam <- spam::spam(0,nrow=ncluster,ncol=dim(cdata)[1])
  ccount <- 0
  for(r in dimnames(cpr)[[1]]) {
    cells <- grep(r,dimnames(cdata)[[1]])
    fit <- kmeans(cdata[cells,],cpr[r,"clusters"],iter.max=10000)
    spam[cbind(fit$cluster+ccount,cells)] <- rep(1,length(fit$cluster))
    ccount <- ccount + cpr[r,"clusters"]
  }
  return(spam)
}
