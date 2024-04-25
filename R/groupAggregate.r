#' groupAggregate
#' 
#' Function which applies aggregation functions over a subset of an array or
#' magpie object.
#' 
#' 
#' @usage groupAggregate(data, vectorfunction=function(x){sum(x,na.rm=TRUE)},
#' dim=3, query=NULL, from=NULL, to=NULL, ...)
#' @param data A MAgPIE object or array
#' @param vectorfunction Aggregation Type. Can be any vector function, and will
#' be applied over the sub-vectors of "dim". E.g. the function is applied to
#' the vector of temperate cereals in each region and each timestep.
#' @param dim The dimension over which is aggregated. Can be 1,3 or higher.
#' This dimension will be replaced by the larger categories specified in the
#' query.
#' @param query Query assigns the name dimension into categories. If Null,
#' query is automatically searched for. If no query can be found, you have to
#' manually select a query (csv or array), as well as specifiy "from" and "to"
#' @param from Only needed if query is not null. From indicates the column
#' within the query where the dimnames of data can be found.
#' @param to Only needed if query is not null. To indicates the column within
#' the query into which the dimnames shall be grouped.
#' @param ... additional arguments for the vectorfunction.
#' @return Returns a magpie object. The name dimension "dim" is grouped
#' according to the query, and the vectorfunction is applied on this groups to
#' aggregate them to each on column.
#' @author Benjamin Bodirsky
#' @export
#' @importFrom utils read.csv
#' @importFrom magclass is.magpie unwrap as.magpie 
#' @seealso \code{\link{colSums}}, \code{\link{superAggregate}}
#' @examples
#' \dontrun{
#' a <- new.magpie(cells_and_regions = c("ZAF.1","ZAR.2","VEN.3"),
#'                 years = c("y1995","y2005"), 
#'                 names = c("Wheat","Barley","Sugar_cane","Sugar_beet"), 
#'                 fill=1:30)
#' 
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
#' b<-groupAggregate(a,dim=3)
#' #> b
#' #An object of class "magpie"
#' #, , tece
#' #
#' #      y1995 y2005
#' #ZAF.1     8    14
#' #ZAR.2    10    16
#' #VEN.3    12    18
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
#' c<-groupAggregate(a,dim=1)
#' #> c
#' #An object of class "magpie"
#' #, , Wheat
#' #
#' #      y1995 y2005
#' #AFR.1     3     9
#' #LAM.2     3     6
#' #
#' #, , Barley
#' #
#' #      y1995 y2005
#' #AFR.1    15    21
#' #LAM.2     9    12
#' #
#' #, , Sugar_cane
#' #
#' #      y1995 y2005
#' #AFR.1    27    33
#' #LAM.2    15    18
#' #
#' #, , Sugar_beet
#' #
#' #      y1995 y2005
#' #AFR.1    39    45
#' #LAM.2    21    24 
#' }

groupAggregate<-function(data, vectorfunction=function(x){sum(x,na.rm=TRUE)}, dim=3, query=NULL, from=NULL, to=NULL, ...) {
  #vectorfunction=function(x){mean(x,na.rm=TRUE)}
  rename<-TRUE
  if (dim>=3) {
    if (is.magpie(data)){data<-unwrap(data)}
    if ((is.null(query))&(length(unique(dimnames(data)[[dim]]))<dim(data)[[dim]])) {
       data<-data
    } else {
      data<-rename_dimnames(data=data,query=query,dim=dim,from=from,to=to)
    }
    data<-as.magpie(data,spatial=1,temporal=2)
    
    groups<-unique(dimnames(data)[[3]])
    out<-array(NA, dim=c(dim(data)[[1]],dim(data)[[2]],length(groups)), dimnames=list(dimnames(data)[[1]],dimnames(data)[[2]],groups))
    
    for (element_x in groups) {
      part_x<-data[,,which(dimnames(data)[[3]]==element_x),drop=FALSE]
      out[,,element_x]<-apply(part_x, MARGIN=c(1,2), FUN=(vectorfunction),...)
    }

    return(as.magpie(out))

  } else if (dim==1){
    # no query exists,
    # a) first check if there are several entries per region --> then simply aggregate within the existing regions
    # b) otherwhise check for a standard query
    
    data<-as.array(suppressMessages(rename_dimnames(data=data,query=query,dim=dim,from=from,to=to)))
    #remove cell numbers if there are any
    dimnames(data)[[1]] <- gsub("\\.[0-9]*$","",dimnames(data)[[1]])
    groups <- sort((unique(dimnames(data)[[1]])))
    
    out<-array(NA, dim=c(length(groups),dim(data)[[2]],dim(data)[[3]]), dimnames=list(groups,dimnames(data)[[2]],dimnames(data)[[3]]))
    
    for (element_x in groups) {
      part_x<-data[dimnames(data)[[1]]==element_x,,,drop=FALSE]
      out[element_x,,]<-apply(part_x, MARGIN=c(2,3), FUN=(vectorfunction),...)
    }
    
    return(as.magpie(out))
  } else {stop("This dim is not yet supported by this function")}
}

