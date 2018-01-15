#' weighted_mean.groupAggregate
#' 
#' Function which aggregates subsets of an array via an weighted mean function.
#' 
#' 
#' @usage weighted_mean.groupAggregate(data, weight, dim, na.rm=TRUE, ...)
#' @param data A MAgPIE object or array
#' @param weight An object with the same dimensions and dimnames as data.
#' @param dim The dimension over which is aggregated. Can be 1,3 or higher.
#' This dimension will be replaced by the larger categories specified in the
#' query.
#' @param na.rm removes NAs in data and weight
#' @param ... additional arguments that will be handed on to groupAggregate().
#' @return Returns a magpie object. The name dimension "dim" is grouped
#' according to the query, and a weighted mean is calculated for each
#' subset(group) of "dim".
#' @author Benjamin Bodirsky
#' @export
#' @importFrom magclass as.magpie
#' @seealso \code{\link{groupAggregate}},\code{\link{colSums}},
#' \code{\link{magpieSums}}, \code{\link{superAggregate}}
#' @examples
#' 
#' a<-new.magpie(cells_and_regions = c("ZAF.1","ZAR.2","VEN.3"), 
#'               years = c("y1995","y2005"), 
#'               names = c("Wheat","Barley","Sugar_cane","Sugar_beet"), 
#'               fill=1:30)
#' #> a
#' #An object of class "magpie"
#' #, , Wheat
#' #
#' #      y1995 y2005
#' #ZAF.1     1     4
#' #ZAR.2     2     5
#' #VEN.3     3     6
#' #
#' #, , Barley
#' #
#' #      y1995 y2005
#' #ZAF.1     7    10
#' #ZAR.2     8    11
#' #VEN.3     9    12
#' #
#' #, , Sugar_cane
#' #
#' #      y1995 y2005
#' #ZAF.1    13    16
#' #ZAR.2    14    17
#' #VEN.3    15    18
#' #
#' #, , Sugar_beet
#' #
#' #      y1995 y2005
#' #ZAF.1    19    22
#' #ZAR.2    20    23
#' #VEN.3    21    24
#' 
#' b<-weighted_mean.groupAggregate(a,weight=a,dim=3)
#' 
#' #> b
#' #An object of class "magpie"
#' #, , tece
#' #
#' #      y1995     y2005
#' #ZAF.1  6.25  8.285714
#' #ZAR.2  6.80  9.125000
#' #VEN.3  7.50 10.000000
#' #
#' #, , sugr_cane
#' #
#' #      y1995 y2005
#' #ZAF.1    13    16
#' #ZAR.2    14    17
#' #VEN.3    15    18
#' #
#' #, , sugr_beet
#' # 
#' #      y1995 y2005
#' #ZAF.1    19    22
#' #ZAR.2    20    23
#' #VEN.3    21    24
#' 
#' 
#' 
weighted_mean.groupAggregate<-function(data,weight,dim,na.rm=TRUE,...) {
  data<-as.magpie(data)
  weight<-as.magpie(weight)
  if (!is.null(weight)){
    if (!identical(c(dimnames(weight)[[1]],dimnames(weight)[[2]]),c(dimnames(data)[[1]],dimnames(data)[[2]]))) {
      stop("Weight and data need to have the same dimnames!")
    }
  }
  if(any(is.na(data))) {
    weight=as.magpie(weight*(!is.na(data)))
  }
  vectorfunction=function(x){sum(x,na.rm=na.rm)}
  product<-data*weight
  out<-groupAggregate(product,vectorfunction=vectorfunction,dim=dim,...)/groupAggregate(weight,vectorfunction=vectorfunction,dim=dim,...)
  return(out)
}