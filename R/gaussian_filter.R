#' gaussian_filter
#' 
#' Smears cellular magpie objects with a gaussian filter.
#' 
#' The calculation is based on the \code{\link[raster]{focal}} function in the
#' \code{raster} library.
#' 
#' @usage gaussian_filter(x,sigma=2,matrix_size=NULL)
#' @param x The object to be smeared.
#' @param sigma The width of the gaussian filter in degree.)
#' @param matrix_size The size of the filter matrix in degree. All cells within
#' the range of the matrix size are taken into account for the smearing. For
#' NULL, the matrix size is 3 times sigma.
#' @return The smeared MAgPIE object.
#' @author Markus Bonsch
#' @export
#' @importFrom magclass is.magpie ncells nyears ndata as.magpie setCells getYears getYears<- getCells setYears setNames mbind as.data.frame
#' @importFrom methods as
#' @seealso
#' \code{\link[raster]{focal}},\code{\link[raster]{focalWeight}},\code{\link{interpolate}}
#' @examples
#' 
#'  \dontrun{a <- gaussian_filter(x=land,sigma=3,matrix_size=4)}
#' 
gaussian_filter<-function(x,sigma=2,matrix_size=NULL){
  
  ###############################
  #Perform some consistency checks on the input
  ################################
  
  if(!is.magpie(x))stop("x has to be a magpie object")
  if(is.null(attr(x,"coordinates"))){
    if(ncells(x)==59199){
      warning("Missing coordinates for the input. Assuming that it is in the standard MAgPIE cellorder.")
    } else {
      stop("Cannot determine coordinates of input object")
    }
  }
  
  #####################################
  #Define necessary subfunctions
  #####################################
  
  #Function to convert MAgPIE objects to rasterLayer objects
  magpie2raster<-function(x,na2zero=TRUE,coverGlobe=TRUE){
    #perform some consistency checks
    if(nyears(x)!=1) stop("Only supported for magpie objects containing a single year")
    if(ndata(x)!=1) stop("Only supported for magpie objects containing a single data column")
    #get the coordinates of the magpie cells
    if(is.null(attr(x,"coordinates"))){
      coord <- getCoordinates(degree=TRUE)
    } else {
      coord<-attr(x,"coordinates")
    }
    
    if(na2zero){
      x[is.na(x)]<-NaN
    }
    #get the data as a vector
    data<-as.vector(x)
    
    if(coverGlobe){
      #get the coordinates of all cells covering the globe
      globeLon<-seq(-189.75,189.75,by=0.5)
      globeLat<-seq(-99.75,99.75,by=0.5)
      globeCoord<-as.matrix(expand.grid(globeLat,globeLon))
      dimnames(globeCoord)[[2]]<-c("lat","lon")
      #check, which cells are missing in the input data
      missing<-setdiff(paste(globeCoord[,"lon"],globeCoord[,"lat"]),paste(coord[,"lon"],coord[,"lat"]))
      missing<-match(missing,paste(globeCoord[,"lon"],globeCoord[,"lat"]))
      #add those coordinates at the end of the coordinates list
      coord<-rbind(coord,globeCoord[missing,c("lon","lat")])
      #add zeros to the data for all missing cells
      data<-c(data,rep(0,length(missing)))    
    }
    
    data <- as.data.frame(data)
    coord<-as.data.frame(coord)
    grid <- suppressWarnings(sp::SpatialPixelsDataFrame(points = coord[c("lon", "lat")], data = data))
    raster<-raster(grid)
    if(na2zero){
      #replace NAs (but not NaNs) by 0
      tmp<-raster@data@values
      tmp[is.na(tmp)&!is.nan(tmp)]<-0
      raster@data@values<-tmp
    }
    return(raster)
  }
  
  #Function to convert rasterLayer objects to MAgPIE objects
  raster2magpie<-function(x,soilcells=TRUE){
    grid<-as(x,"SpatialPixelsDataFrame")
    sp::fullgrid(grid)<-TRUE
    out<-grid@data[[1]]
    names(out)<-paste("GLO",1:length(grid@data[[1]]),sep=".")
    out<-as.magpie(out)
    tmp<-coordinates(grid)
    dimnames(tmp)[[2]]<-c("lon","lat")
    attr(out,"coordinates")<-tmp
    if(soilcells){
      #get the coordinates of the standard 59199 magpie cells
      ref<-getCoordinates(degree=TRUE)
      #match the reference coordinates with the actual ones
      order<-match(paste(ref[,1],ref[,2]),paste(tmp[,1],tmp[,2]))
      out<-setCells(out[order,,],paste("GLO",1:59199,sep="."))
      attr(out,"coordinates")<-ref
    }
    return(out)
  }
  
  ############################################
  #Do the aggregation
  ############################################
  
  #get the gaussian weight matrix
  tmp<-magpie2raster(x[,1,1])
  if(is.null(matrix_size)){
    filter<-raster::focalWeight(x=tmp,d = sigma,type = "Gauss")
  } else {
    filter<-raster::focalWeight(x=tmp,d = c(sigma,matrix_size),type = "Gauss")
  }  
                      
  #Creatre output object
  outtmp<-list()
  out<-NULL
  revert_years<-FALSE
  if(is.null(getYears(x))){
    revert_years<-TRUE
    getYears(x)<-"y9999"
  }
  
  #Perform the filtering
  for(t in getYears(x)){
    for(d in 1:ndata(x)){
      tmp<-magpie2raster(x[,t,d])
      tmp<-raster::focal(tmp,w=filter,na.rm=FALSE)
      tmp<-raster2magpie(tmp,soilcells=TRUE)
      tmp<-setCells(tmp,getCells(x))
      tmp<-setYears(tmp,t)
      tmp<-setNames(tmp,magclass::getNames(x)[d])
      coordinates<-attr(tmp,"coordinates")
      outtmp[[t]]<-mbind(outtmp[[t]],tmp)
    }
    out<-mbind(out,outtmp[[t]])
  }
  attr(out,"coordinates")<-coordinates
  if(revert_years) getYears(out)<-NULL
  return(out)
}