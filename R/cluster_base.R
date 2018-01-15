#' Cluster Base
#' 
#' Reads a series of MAgPIE files and combines them to a matrix which is then
#' used by mag_kmeans or mag_hierarchical for calculating a clustering
#' 
#' As this procedure is typically quite time consuming the function is a the
#' begining looking for a file cluster.Rdata in the input folder and is loading
#' the data from this file instead of cfiles if it exists. If it does not exist
#' data is read from cfiles but written to cluster.Rdata so that it is already
#' available in a second attempt.
#' 
#' @param ifolder input folder where the MAgPIE input files are located
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
#' @return A matrix containing the data
#' @author Jan Philipp Dietrich
#' @seealso \code{\link{mag_kmeans}}, \code{\link{mag_hierarchical}},
#' \code{\link{clusterspam}}
cluster_base <- function(ifolder=".",cfiles=c("lpj_yields_rf", "lpj_yields_ir", "lpj_airrig", "transport_distance"),years2use="y1995", spatial_header=NULL, use_cache=TRUE) {
  rdata_file <- path(ifolder,"cluster.Rdata")
  if(file.exists(rdata_file) & use_cache) {
    load(rdata_file)
  } else {
    cdata <- NULL    
    for(f in cfiles) {
      tmp <- read.magpie(paste(path(ifolder,f),"*",sep=""))
      if(nyears(tmp)>1) tmp <- tmp[,years2use,]
      tmp <- wrap(tmp,list(1,c(2,3)))
      cdata <- cbind(cdata,tmp)
    }
    unlink(tmp)
    cdata <- scale(cdata)
    try(save(cdata,file=rdata_file))
    
  }
  if(!is.null(spatial_header))  dimnames(cdata)[[1]] <- spatial_header
  return(cdata)  
}
