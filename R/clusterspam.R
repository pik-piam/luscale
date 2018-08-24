#' Main MAgPIE clustering function
#' 
#' This is the main MAgPIE clustering function which you should use if you want
#' to create a spam relation matrix based on clustering methods
#' 
#' The Spam relation matrix is written at the same time in the output folder.
#' 
#' @param lr low resolution containing the clustering mode as first letter in
#' the name, e.g. h100. Available modes are at the moment n: normed kmeans
#' clustering, h: hierarchical complete linkage clustering, s: hierarchical
#' single linkage clustering, w: hierarchical Ward clustering and c: combined
#' hierarchical/k-means-clustering. In the latter hierarchical complete-linkage
#' clustering is used to determine the distribution of clusters among regions
#' and normed-k-means clustering is used to determine clusters within regions
#' @param hr high resolution, e.g. 0.5
#' @param ifolder input folder where the MAgPIE input files are located
#' @param ofolder output folder where spam relation matrix should be written
#' to.
#' @param cfiles a vector containin the names of the MAgPIE input files
#' (beginning of the name is enough)
#' @param years2use A vector with years with should be taken into account for
#' the clustering
#' @param spatial_header A vector of the form c("REG.1","REG.2") (region name,
#' cell number) with entries for each spatial entity of the MAgPIE input files
#' which should be used to replace the names given in the inputs (required for
#' flexible region aggregation as here region names might change.). If set to
#' NULL the original information is used.
#' @param use_cache Read data from cache file if available (dangerous as changes
#' in settings will not be considered if an existing cache file is found).
#' @param weight named vector with weighting factors for each region for the cluster distribution 
#' ,e.g. weight=c(AFR=3,EUR=0.5). weight > 1 will grant more cluster to a region and
#' weight < 1 less cluster than by default. 
#' @return A spam relation matrix
#' @author Jan Philipp Dietrich
#' @export
#' @importFrom magclass read.magpie nyears wrap
#' @importFrom spam rowSums
#' @seealso \code{\link{cluster_per_region}}, \code{\link{mag_kmeans}},
#' \code{\link{mag_hierarchical}}
clusterspam <- function(lr,hr="0.5", ifolder=".", ofolder=".", cfiles=c("lpj_yields", "lpj_airrig", "transport_distance"), years2use="y1995", spatial_header=NULL, use_cache=TRUE, weight=NULL) {
  mode <- substr(lr,0,1)
  ncluster <- as.integer(substring(lr,2))
  cdata <- cluster_base(ifolder,cfiles,years2use,spatial_header,use_cache)

  if(mode=="n") {
    spam <- mag_kmeans(cdata,ncluster,weight)
  } else if(mode=="h" | mode=="w" | mode=="s") {
    spam <- mag_hierarchical(cdata,ncluster,ifolder,mode,weight)
  } else if(mode=="c"){
    calcCPR <- function(spam, cell2reg) {
      reg <- unique(cell2reg)
      cluster2reg <- as.factor(spam%*%as.numeric(cell2reg)/rowSums(spam))
      levels(cluster2reg) <- levels(cell2reg)
      cpr <- t(rbind(table(cell2reg),table(cluster2reg)))
      dimnames(cpr)[[2]] <- c("cells","clusters")
      return(cpr)  
    }
    tmpspam <- mag_hierarchical(cdata,ncluster,ifolder,mode="h",weight)
    cell2reg <- as.factor(sub("\\..*$","",dimnames(cdata)[[1]]))
    spam <- mag_kmeans(cdata,cpr=calcCPR(tmpspam,cell2reg))
  } else {
    stop("Unkown clustering mode ",mode,"!")
  }
  wkey <- ifelse(is.null(weight), "", gsub(".","",paste0("_",names(weight),weight,collapse=""),fixed=TRUE))
  write.spam(spam,path(ofolder,paste(hr,"-to-",lr,wkey,"_sum.spam",sep="")))
  saveRDS(spam2mapping(spam,rownames(cdata)), path(ofolder,paste(hr,"-to-",lr,wkey,"_mapping.rds",sep="")))
  return(spam)
}
