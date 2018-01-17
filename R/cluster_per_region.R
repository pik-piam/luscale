#' Cluster per Region
#' 
#' This function calculates an appropriate number of clusters per region as it
#' is needed for mag_kmeans
#' 
#' @param cdata a cluster data file as produced by cluster_base
#' @param ncluster The desired total number of clusters.
#' @param weight named vector with weighting factors for each region for the cluster distribution 
#' ,e.g. weight=c(AFR=3,EUR=0.5). weight > 1 will grant more cluster to a region and
#' weight < 1 less cluster than by default. 
#' @return A matrix with regions in rows and number of cells and clusters in
#' columns
#' @author Jan Philipp Dietrich
#' @seealso \code{\link{mag_kmeans}}, \code{\link{cluster_base}},
#' \code{\link{clusterspam}}
cluster_per_region <- function(cdata,ncluster,weight=NULL){
  tmp <- dimnames(cdata)[[1]]
  regions <- unique(sub("\\..*$","",tmp))
  if(length(regions) > ncluster) stop("More regions than cluster. Clustering stopped!")
  calcw <- function(weight, regions) {
    w <- rep(1,length(regions))
    names(w) <- regions
    if(!is.null(weight)) w[names(weight)] <- weight
    return(w)
  }
  cpr <- rep(NA,length(regions))
  names(cpr) <- regions
  for(r in regions) cpr[r] <- length(grep(r,tmp))
  cpr <- cbind(cpr,cpr,calcw(weight,regions))
  dimnames(cpr)[[2]] <- c("cells","clusters","weight")
  cpr[,"clusters"] <- round(cpr[,"cells"]*cpr[,"weight"]/sum(cpr[,"cells"]*cpr[,"weight"])*ncluster)
  cpr[cpr[,"clusters"]==0,"clusters"] <- 1
  while(sum(cpr[,"clusters"])!=ncluster) {
    correct <- ncluster-sum(cpr[,"clusters"])
    m <- which(cpr[,"clusters"]==max(cpr[,"clusters"]))[1]
    cpr[m,"clusters"] <- cpr[m,"clusters"] + sign(correct)
  }
  return(cpr)
}

