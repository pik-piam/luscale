path <- function(...,ftype=NULL) {
  if(!is.null(ftype)) if(substring(ftype,1,1)!='.') ftype <- paste('.',ftype,sep='')
  out <- paste(...,sep="/")
  out <- paste(gsub("//","/",out),ftype,sep="")
  first <- list(...)[[1]]
  .tmp <- function(first,out) {
    manipulate <- FALSE
    if(is.null(first)) manipulate <- TRUE
    else if(first=="") manipulate <- TRUE
    if(manipulate) out <- gsub("^/+","",out)
    return(out)
  }
  if(length(first)>1) {
    for(i in 1:length(first)) out[i]<- .tmp(first[i],out[i])
  } else {
    out <- .tmp(first,out)
  }
  return(out)
}