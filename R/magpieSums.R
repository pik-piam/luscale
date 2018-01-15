#' Different summations of MAgPIE-objects
#' 
#' Returns the sum over all cells belonging to a specific aggregation unit
#' (countries,regions,world) of a 0.5 degree MAgPIE data set.
#' 
#' 
#' @usage magpieSums(x, level="country", year="y1995")
#' @param x A 0.5 degree MAgPIE-object
#' @param level Summation level. Currently available: "country", "region" and
#' "world"
#' @param year The year for which data should be summed up
#' @return A matrix whith aggregation units as rows and data sums in columns
#' @author Jan Philipp Dietrich
#' @export
#' @importFrom magclass is.magpie ncells ndata getNames nregions getRegions
#' @seealso \code{\link{superAggregate}}
#' @examples
#' 
#'   #magpieSums(x)
#' 
magpieSums <- function(x,level="country",year="y1995") {
  warning("magpieSums is deprecated! Please use alternative functions such as dimSums instead!") 
  if(!is.magpie(x)) stop("Wrong data format! Input is not an MAgPIE-Object!")
  if(ncells(x)!=59199) stop("Wrong number of cells. Data must be provided as 0.5 degree data set!")
  if(level=="country") {
    stop("Country support has been dropped in version 2.9 of the package. Please use dimSums instead!")
  } else if(level=="region") {
    y <- matrix(0,nregions(x),ndata(x))
    rownames(y) <- getRegions(x)
    colnames(y) <- getNames(x)
    for(c in getRegions(x)) {
      y[c,] <- colSums(x[c,year,])
    }
  } else if(level=="world"){
    y <- colSums(x[,year,])
    dimnames(y)[[1]] <- list("GLO")
  } else {
    stop("Unknown aggregation level!")
  }       
  return(y)
}
