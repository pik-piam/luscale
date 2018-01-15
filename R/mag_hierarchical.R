#' MAgPIE Hierarchical clustering
#' 
#' Performs MAgPIE hierarchical clustering and calculates corresponding spam
#' relation matrix
#' 
#' As the creation of a clustering tree is very time consuming the function
#' checks first in the input folder if the corresponding data already exists
#' and if not it stores the tree information in the input folder for later use
#' in the next execution of this function.
#' 
#' @usage mag_hierarchical(cdata,ncluster,ifolder,mode="h")
#' @param cdata a cluster data file as produced by cluster_base
#' @param ncluster The desired total number of clusters.
#' @param ifolder input folder where the MAgPIE input files are located
#' @param mode Clustering type. At the moment you can choose between complete
#' linkage clustering (h), single linkage clustering (s) and Ward clustering
#' (w).
#' @return A spam relation matrix
#' @author Jan Philipp Dietrich
#' @importFrom magclass getCells ncells getRegions 
#' @importFrom stats hclust cutree
#' @seealso \code{\link{cluster_per_region}}, \code{\link{mag_kmeans}},
#' \code{\link{clusterspam}}
mag_hierarchical <- function(cdata,ncluster,ifolder,mode="h") {
  tdata_file <- path(ifolder,paste(mode,digest::digest(dimnames(cdata)[[1]],"md5"),"tree.Rdata",sep="_"))
  if(file.exists(tdata_file)) {
    load(tdata_file)
  } else {
    spam <- spam::spam(0,nrow=dim(cdata)[1],ncol=ncluster)
    fullfit <- list()
    fullfit$labels <- getCells(cdata)
    fullfit$order <- rep(NA,ncells(cdata))
    names(fullfit$order) <- fullfit$labels
    nrow <- 0
    nrows <- NULL
    for(r in getRegions(cdata)) {
      cells <- grep(r,dimnames(cdata)[[1]])
      dist <- dist(cdata[cells,], method = "euclidean")
      if(mode=="h") {
        fit <- hclust(dist, method="complete")
      } else if(mode=="w") {
        fit <- hclust(dist^2, method="ward")
      } else if(mode=="s") {
        fit <- hclust(dist, method="single")
      }
      cellnumber <- as.numeric(sub("^.*\\.","",fit$labels))
      fit$merge[fit$merge<0] <- cellnumber[fit$merge[fit$merge<0]*-1]*-1
      fit$merge[fit$merge>0] <- fit$merge[fit$merge>0]+nrow
      fit$order <- cellnumber[fit$order]
      
      fullfit$merge <- rbind(fullfit$merge,fit$merge)
      fullfit$height <- c(fullfit$height,fit$height)
      
      fullfit$order[fit$labels] <- fit$order
      
      nrow <- nrow(fullfit$merge)
      nrows <- c(nrows,nrow)   
    }
    #add links between regions with huge heights
    fullfit$height <- c(fullfit$height,rep(2*max(fullfit$height),length(nrows)-1))
    for(n in nrows[-length(nrows)]) {
      fullfit$merge <- rbind(fullfit$merge,c(n,dim(fullfit$merge)[1]))  
    }  
    #bring data in the right order, adapt row numbers accordingly
    o <- order(fullfit$height)
    fullfit$merge <- fullfit$merge[o,]
    b <- order(fullfit$merge)
    b <- b[min(which(fullfit$merge[b]>0)):length(b)] 
    fullfit$merge[b] <- order(o)[-length(o)]
    fullfit$height <- fullfit$height[o]
    attr(fullfit,"class") <- "hclust"
    save(fullfit,file=tdata_file)
  }
  clusters <- cutree(fullfit,k=ncluster)
  #sort clusters by regions
  cl <- NULL
  for(r in getRegions(cdata)) {
    cl <- c(cl,unique(clusters[grep(paste0(r,"."),names(clusters),fixed=TRUE)]))
  }
  if(length(cl)!=ncluster) stop("Something went wrong during the clustering. Some cluster seem to exist across region borders or there are less clusters than demanded!")
  tmp <- order(cl)[clusters]
  names(tmp) <- names(clusters)
  clusters <- tmp
  spam <- spam::spam(0,nrow=ncluster,ncol=dim(cdata)[1])
  spam[cbind(clusters,1:dim(cdata)[1])] <- rep(1,dim(cdata)[1])
  return(spam)
}
