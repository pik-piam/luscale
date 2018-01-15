#' Reshape File
#' 
#' Disaggregates a cellular MAgPIE output to 0.5 degree based on spam matrices.
#' 
#' 
#' @usage reshape_file(filename,spam=NULL,search_folders=NULL)
#' @param filename file name of a file that should be disaggregated or an R
#' object to be disaggregated.
#' @param spam a vector of spam files and the corresponding needle. For more
#' details please check the sourcecode of this function with the predefined
#' spam vector. In case of an R object to be disagregated, no needle is
#' required.
#' @param search_folders Folders in which the script should search for spam
#' files. IMPORTANT: the paths are meant relative to the location where the
#' file given with fielname is stored! Not used if an R object is to be
#' disaggregated.
#' @return NULL, if no disaggregation rule is found, otherwise the
#' disaggregated MAgPIE object
#' @author Jan Philipp Dietrich,Markus Bonsch
#' @export
#' @importFrom magclass read.magpie getYears mbind
#' @seealso \code{\link{reshape_folder}}, \code{\link{speed_aggregate}}
#' @examples
#' 
#'  \dontrun{a <- reshape_file("cell.land.m.csv")}
#' \dontrun{a <- read.magpie("cell.land.m.csv")}
#' \dontrun{b <- reshape_file(a,spam="0.5_to_h100_area_weighted_mean.spam")}
#' 
reshape_file <- function(filename,spam=NULL,search_folders=NULL) {
  folder <- gsub("[^/]+$","",filename)

  # Check, if filename is a file or an R object.
  isfile<-FALSE
  if(is.character(filename)) isfile<-TRUE 
  if(!isfile & is.null(spam)) stop("Input is an R object. A spam file has to be supplied")
  
  if(isfile){
    # Determine spam file locations relative to the location of the data file
    # that should be reshaped
    if(is.null(spam)) {
      spam["m"]   <- "*-to-*_sum.spam"
      spam["aw"]  <- "*-to-*_area_weighted_mean.spam"
      spam["caw"] <- "*-to-*_crop_weighted_mean_%year.spam"
      spam["iaw"] <- "*-to-*_irrig_area_weighted_mean.spam"
    }
  
    if(is.null(search_folders)){
      search_folders <- c(".", "../../../input/cellular", "../../input/cellular",
                          "../input/cellular", "input/cellular")
    }
    
    #search for needle in filename
    found <- sapply(paste(".",names(spam),".",sep=""),grepl,filename,fixed=TRUE)
    if(sum(found)>1) stop("More than one type indicator found in the file name!")
    if(sum(found)==0) return(NULL)
  
    #create path to spam file based on search folders
    for(s in search_folders) {
      if(length(Sys.glob(gsub("%year","*",path(folder,s,spam[found]),fixed=TRUE)))>0) {
        spam <- path(folder,s,spam[found])
        break
      }
    }
    if(length(spam)>1) stop("Required spam file ",spam[found]," not found!")
  }

  if(grepl("%year",spam,fixed=TRUE)) {
    if(isfile){
      m <- read.magpie(filename)
    } else{
      m<-filename
    }
    years <- getYears(m)
    tmp <- NULL
    for(y in years) {
      s <- read.spam(gsub("%year",y,spam,fixed=TRUE))
      tmp <- mbind(tmp,speed_aggregate(m[,y,],t(s)))
    }
    return(tmp)
  } else {
    s <- read.spam(spam)
    return(speed_aggregate(filename,t(s)))
  } 
} 
