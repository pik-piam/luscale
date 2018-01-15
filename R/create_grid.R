#' Create Grid
#' 
#' Create a spam relation matrix between two regular grid resolutions (measured
#' in degree)
#' 
#' 
#' @usage create_grid(res_in,res_out,fname=NULL,folder=NULL)
#' @param res_in input resolution in degree (e.g. 0.5)
#' @param res_out output resolution in degree (e.g. 12)
#' @param fname file name of a file the output should be written to. If
#' fname=NULL the relation matrix is returned by the function, if
#' fname="default" the default name format res_in-to-res_out_sum.spam is used
#' @param folder folder the file should be written to
#' @return If fname=NULL, the relation matrix, otherwise nothing is returned.
#' @author Jan Philipp Dietrich
#' @export
#' @importFrom magclass getCPR read.magpie is.magpie 
#' @seealso \code{\link{create_spam}}, \code{\link{read.spam}},
#' \code{\link{write.spam}}
#' @examples
#' 
#'  \dontrun{rel <- create_grid(0.5,12.0)}
#' 
create_grid <- function(res_in,res_out,fname=NULL,folder=NULL) {
  coordinates_int <- getCoordinates(degree=FALSE)
  cpr <- getCPR(res_in)

  max_lon <- c(360/res_in,360/res_out)
  max_lat <- c(180/res_in,180/res_out)

  output <- NULL
  map <- array(" ",c(length(cpr),max_lon[2],max_lat[2]))
  dimnames(map)[[1]] <- names(cpr)
  cell_count <- 0
  for(region in names(cpr)) {
    for(i in (cell_count + (1:cpr[region]))) {
      x <- ceiling(coordinates_int$lon[i]/res_out*res_in)
      y <- ceiling(coordinates_int$lat[i]/res_out*res_in)
      map[region,x,y] <- ifelse((map[region,x,y] == " "),i,paste(map[region,x,y],i,sep=";"))  
    }
    for(x in 1:max_lon[2]) {
      for(y in 1:max_lat[2]) {
      if(map[region,x,y]!= " ") {
        output <- c(output, x, y,map[region,x,y])
        }
      }
    }
    cell_count <- cell_count + cpr[region]
  }
  
  output <- t(array(output,c(3,length(output)/3)))
  
  relation_list <- list()
  relation_list$indices <- matrix(0,sum(cpr),2)
  relation_list$values <- rep(1,sum(cpr))
  m <- 1
  for(i in 1:dim(output)[1]) {
    highres_cells <- as.integer(strsplit(as.character(output[i,3]),";")[[1]])
    for(j in 1:length(highres_cells)) {
      relation_list$indices[m,1] <- i
      relation_list$indices[m,2] <- highres_cells[j]
      m <- m + 1 
    }
  }
  if(!is.null(fname)) {
    if(fname=="default") fname <- paste(format(as.numeric(res_in),nsmall=1),"-to-",format(as.numeric(res_out),nsmall=1),"_sum.spam",sep="") 
    write.spam(spam::spam(relation_list),path(folder,fname))
  } else {
    return(spam::spam(relation_list))
  }
}