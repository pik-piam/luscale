getCoordinates <- function(degree=FALSE) {
  res <- 0.5
  out <- ludata_dat
  if(degree) {
    out$lon <- (out$lon-0.5)*res-180
    out$lat <- (out$lat-0.5)*res-90  
  }
  return(as.matrix(out))
}