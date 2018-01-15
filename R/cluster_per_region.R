#' Cluster per Region
#' 
#' This function calculates an appropriate number of clusters per region as it
#' is needed for mag_kmeans
#' 
#' 
#' @usage cluster_per_region(cdata,ncluster)
#' @param cdata a cluster data file as produced by cluster_base
#' @param ncluster The desired total number of clusters.
#' @return A matrix with regions in rows and number of cells and clusters in
#' columns
#' @author Jan Philipp Dietrich
#' @seealso \code{\link{mag_kmeans}}, \code{\link{cluster_base}},
#' \code{\link{clusterspam}}
cluster_per_region <- function(cdata,ncluster){
  tmp <- dimnames(cdata)[[1]]
  regions <- unique(sub("\\..*$","",tmp))
  if(length(regions) > ncluster) stop("More regions than cluster. Clustering stopped!")
  cpr <- rep(NA,length(regions))
  names(cpr) <- regions
  for(r in regions) cpr[r] <- length(grep(r,tmp))
  cpr <- cbind(cpr,cpr)
  dimnames(cpr)[[2]] <- c("cells","clusters")
  cpr[,"clusters"] <- round(cpr[,"cells"]/sum(cpr[,"cells"])*ncluster)
  cpr[cpr[,"clusters"]==0,"clusters"] <- 1
  while(sum(cpr[,"clusters"])!=ncluster) {
    correct <- ncluster-sum(cpr[,"clusters"])
    m <- which(cpr[,"clusters"]==max(cpr[,"clusters"]))[1]
    cpr[m,"clusters"] <- cpr[m,"clusters"] + sign(correct)
  }
  return(cpr)
}

