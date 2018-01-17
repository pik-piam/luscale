#' MAgPIE Kmeans clustering
#' 
#' Performs MAgPIE kmeans clustering and calculates corresponding spam relation
#' matrix
#' 
#' 
#' @param cdata a cluster data file as produced by cluster_base
#' @param ncluster The desired total number of clusters.
#' @param weight named vector with weighting factors for each region for the cluster distribution 
#' ,e.g. weight=c(AFR=3,EUR=0.5). weight > 1 will grant more cluster to a region and
#' weight < 1 less cluster than by default. 
#' @return A spam relation matrix
#' @author Jan Philipp Dietrich
#' @importFrom stats kmeans
#' @seealso \code{\link{cluster_per_region}}, \code{\link{mag_hierarchical}},
#' \code{\link{clusterspam}}
mag_kmeans <- function(cdata,ncluster,weight=NULL) {
  cpr <- cluster_per_region(cdata,ncluster,weight)
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
