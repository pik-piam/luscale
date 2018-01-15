#' Create SPAM
#' 
#' Creates a spam relation matrix based on another relation matrix and a 1D
#' MAgPIE object (only 1 year and 1 data column)
#' 
#' 
#' @usage create_spam(x,rel,year=NULL,elem=NULL,fname=NULL)
#' @param x Vector, MAgPIE object or file name of a MAgPIE object. If MAgPIE
#' object it must either contain only a single year and a single data column or
#' year and elem have to be defined.
#' @param rel Input spam relation matrix or file name of a input spam relation
#' matrix
#' @param year Year which should be used for calculation (has to be specified
#' if the MAgPIE object contains more than 1 year)
#' @param elem Data column which should be used for calculation (has to be
#' specified if the MAgPIE object contains more than 1 data column)
#' @param fname file name of a file the output should be written to. If
#' fname=NULL the spam matrix is returned by the function
#' @return If fname=NULL, the spam matrix, otherwise nothing is returned.
#' @author Jan Philipp Dietrich
#' @export
#' @seealso \code{\link{speed_aggregate}}, \code{\link{read.spam}},
#' \code{\link{write.spam}}
#' @examples
#' 
#'  \dontrun{rel2 <- create_spam(x,rel)}
#' 

create_spam <- function(x,rel,year=NULL,elem=NULL,fname=NULL) {
  if(is.character(x)) x <- read.magpie(x)  
  if(is.character(rel)) rel <- read.spam(rel)
  if(!is.null(year)) x <- x[,year,]
  if(!is.null(elem)) x <- x[,,elem] 

  if(is.magpie(x) | is.array(x)) {
    if(dim(x)[2]!=1) stop("Too many years, function can only handle a single year!")
    if(dim(x)[3]!=1) stop("Too many data sets, function can only handle a single data set!")
  }
  x <- as.numeric(x)
  y <- as.numeric(rel %*% x)
  
  if(any(y==0)) {
    x <- x + 10^-100 
    y <- as.numeric(rel %*% x)
  }
  out <- spam::diag.spam(1/y) %*% rel %*% spam::diag.spam(x)
  if(!is.null(fname)) {
    write.spam(out,fname)
  } else {
    return(spam::diag.spam(1/y) %*% rel %*% spam::diag.spam(x))
  }
}