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
#' @param cpr cells-per-region information as returned by cluster_per_region. Weight and ncluster are
#' ignored in case that cpr is provided!
#' @param seed a single value, interpreted as an integer, or NULL, to define seed for random calculations
#' @return A spam relation matrix
#' @author Jan Philipp Dietrich
#' @importFrom stats kmeans
#' @seealso \code{\link{cluster_per_region}}, \code{\link{mag_hierarchical}},
#' \code{\link{clusterspam}}
mag_kmeans <- function(cdata,ncluster=NULL,weight=NULL,cpr=NULL,seed=42) {
  if(is.null(cpr)) cpr <- cluster_per_region(cdata,ncluster,weight)
  spam <- spam::spam(0,nrow=sum(cpr[,"clusters"]),ncol=dim(cdata)[1])
  ccount <- 0
  set.seed(seed)
  for(r in dimnames(cpr)[[1]]) {
    cells <- grep(r,dimnames(cdata)[[1]])
    fit <- kmeans(cdata[cells,],cpr[r,"clusters"],iter.max=10000)
    spam[cbind(fit$cluster+ccount,cells)] <- rep(1,length(fit$cluster))
    ccount <- ccount + cpr[r,"clusters"]
  }
  set.seed(NULL)
  return(spam)
}
