#' getAggregationMatrix
#' 
#' Function which converts the supplied regionmapping file to a transformation
#' matrix which then can be used for aggregation with
#' \code{\link{speed_aggregate}}.
#' 
#' 
#' @usage getAggregationMatrix(rel, from=NULL, to=NULL)
#' @param rel file name of a region mapping (".csv" as file ending, ";" used as
#' separator and with columns containing the region names) or a mapping matrix
#' (matrix with region names in columns)
#' @param from Name of the first column to be used (if not set the first or
#' second column will be used)
#' @param to Name of the second column to be used (if not set the second or
#' third column will be used)
#' @return A matrix nregions1 x nregions2 with 1s and 0s showing the mapping of
#' countries to regions.
#' @author Jan Philipp Dietrich
#' @export
#' @examples
#' \dontrun{
#'   x <- cbind(reg=c("REG1","REG2","REG2","REG3"),country=c("C1","C2","C3","C4")) 
#' 
#'   getAggregationMatrix(x)
#'   
#'   getAggregationMatrix(x,from="reg",to="reg")
#' }
getAggregationMatrix <- function(rel,from=NULL,to=NULL) {
  
  if(!(is.matrix(rel) | is.data.frame(rel))) {
    if(!file.exists(rel)) stop("Cannot find given region mapping file!")
    rel <- read.csv(rel, as.is = TRUE, sep = ";")     
  }
    
  if(is.null(from)) {
    from <- ifelse(dim(rel)[2]==3,2,1)
  }
  if(is.null(to)) {
    to <- ifelse(dim(rel)[2]==3,3,2)
  }
  
  regions <- unique(rel[,to])
  countries <- unique(rel[,from])
  m <- matrix(data=0, nrow=length(regions),ncol=length(countries),dimnames=list(regions=regions,countries=countries))
  m[cbind(match(rel[,to],rownames(m)),match(rel[,from],colnames(m)))] <- 1
  if(is.numeric(to)) to <- dimnames(rel)[[2]][to]
  if(is.numeric(from)) from <- dimnames(rel)[[2]][from]
  names(dimnames(m)) <- c(to,from)
  return(m)
}