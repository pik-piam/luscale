#' superAggregate
#' 
#' Function which applies aggregation functions over a subset of an magpie or
#' lpj object.
#' 
#' 
#' @aliases superAggregate superAggregate,lpj-method
#' @usage superAggregate(data, aggr_type, level="reg", weight=NULL, na.rm=TRUE,
#' crop_aggr=FALSE, ...)
#' @param data An MAgPIE or LPJ object
#' @param aggr_type Aggregation Type. Can be any function for one or two
#' objects (data and weight) of the same size. Currently pre-supported
#' functions: "sum","mean","weighted_mean".
#' @param level Aggregation level: Either a level type (as a name) or a vector
#' of 3-Character region names with the same length as the data has cells.
#' Allowed level types are global "glo", regional "reg", per Country "country"
#' and per REMIND regions "remind_reg".  "country" and "remind_reg" are only
#' supported for 0.5 grid datacells with 59199 cells and are always returned as
#' arrays. If you use a vector of regions the aggregation will take place
#' according to your regions.
#' @param weight Currently only used for weighted_mean (see
#' \code{\link{weighted.mean}}, yet also applicable for individualized
#' functions. Has to be of the same size as data.
#' @param na.rm If TRUE, NAs are ignored both in data and weight.
#' @param crop_aggr determines whether output should be crop-specific (FALSE)
#' or aggregated over all crops (TRUE). The method used for aggregation is set
#' by aggr_type (Currently works only for levels "reg" and "glo")
#' @param ... additional arguments for the aggregation method for the standard
#' functions (not for self-created ones)
#' @return In the case of level="glo" or "reg", the function returns a MAgPIE
#' object. In the case of level="country" an array is returned. In the case of
#' an LPJ object being aggregated, a list of MAgPIE objects is returned, each
#' entry being one of the 4th dimension slices.
#' @author Benjamin Bodirsky, Jan Philipp Dietrich, Florian Humpenoeder
#' @export
#' @importFrom magclass as.magpie
#' @seealso \code{\link{colSums}}, \code{\link{magpieSums}}
#' @examples
#' 
#' data(population_magpie)
#' superAggregate(population_magpie,"sum",level="glo")
#' superAggregate(population_magpie,"mean",level="glo")
#' superAggregate(population_magpie,"weighted_mean",level="glo",weight=population_magpie)
#' aggregation_function<-function(func_data,func_weight) {
#'   colMeans(func_data)
#' }
#' superAggregate(population_magpie,aggregation_function,level="glo",weight=population_magpie)
#' 
superAggregate<-function(data, aggr_type, level="reg", weight=NULL, na.rm=TRUE, crop_aggr=FALSE, ...) {

  aggr_vector <- NULL

  if(length(level)>1) {  #vector instead of single element
    aggr_vector <- level
    level <- "reg"
  }
  
  if (is.null(weight)==FALSE) {
    if (identical(dim(weight),dim(data))==FALSE) {
       stop("data and weight have to have the same dimensions!")
    } else {
       weight<-as.magpie(weight)
       weight <- magpiesort(weight)
    } 
  }
  
  data <- magpiesort(data)
  data("luqueries", envir=environment() ,package = "luscale")
  
  if ((level=="reg")&(all(getRegions(data)%in%luqueries$spatial$iso_reg[,"iso"]))){
    if (length(getRegions(data))!=dim(data)[1]) {stop("some countries are represented more than once")}
    data<-as.array(data)
    dimnames(data)[[1]]<-substr(dimnames(data)[[1]],1,3)
    data<-rename_dimnames(data,query=NULL,dim=1,from="iso",to="reg")
    data<-data[order(dimnames(data)[[1]]),,,drop=FALSE]
    data<-as.magpie(data)
    if (is.null(weight)==FALSE) {
      weight<-as.array(weight)    
      dimnames(weight)[[1]]<-substr(dimnames(weight)[[1]],1,3)
      weight<-rename_dimnames(weight,query=NULL,dim=1,from="iso",to="reg")
      weight<-weight[order(dimnames(weight)[[1]]),,,drop=FALSE]
      weight<-as.magpie(weight)    
    }
    level<-"reg"
  }

  if(level == "remind_reg") {
    if(dim(data)[1]!=59199) stop("Level 'remind_reg' is only allowed for 0.5deg data sets! Dimension of input data shows a wrong number of elements in first dimension!")
    aggr_vector <- as.character(landusedata$cellbelongings$remind.regions)
    level <- "reg"
  }

  ### Emulation of array method  
  if (is.magpie(data)==FALSE) { 
    data<-as.magpie(data)
#    retour<-callGeneric(data=data, aggr_type=aggr_type, level = level, weight = weight, aggregation_function = aggregation_function)
    if((level=="reg")|(level=="glo")){warning("array transformed into a MAgPIE object")}
  }
  
  #replace regions with aggr_vector
  if(!is.null(aggr_vector)) {
    if(dim(data)[1]!=length(aggr_vector)) stop("Level vector has wrong number of elements! Level and Data do not fit!")
    if(any(nchar(aggr_vector)!=3)) stop("Level vector contains elements with more or less than 3 characters! This is not allowed at the moment!")
    dimnames(data)[[1]] <- paste(aggr_vector,1:ncells(data),sep=".")
  }

  ### Checks

  if (is.null(weight)==FALSE) {
    if (!is.magpie(weight)) {stop("Data Format of weight not supported")}
    if (identical(dim(weight),dim(data))==FALSE) {
      stop("data and weight have to have the same dimensions!")
    }  
  }   


  if (!is.magpie(data)) {stop("Data Format of data not supported")}


  if ((level=="country")&(ncells(data)!=59199)) {
      stop("Wrong number of cells. Data must be provided as 0.5 degree data set!")
    }
    
  if(!is.null(aggr_vector)) {
  
  }


  ### Define Aggregation Function
  if (is.function(aggr_type)==TRUE) {
    aggregation_function<-aggr_type
  } else {
    if (aggr_type=="sum") {
      if (crop_aggr) {
        aggregation_function<-function(func_data,func_weight) {
          rowSums(colSums(func_data,na.rm=na.rm),na.rm=na.rm)
        }
      } else {
        aggregation_function<-function(func_data,func_weight) {
          colSums(func_data,na.rm=na.rm)
        }  
      }
    } else if (aggr_type=="mean"){
      if (crop_aggr) {
        aggregation_function<-function(func_data,func_weight) {
          rowMeans(colMeans(func_data,na.rm=na.rm),na.rm=na.rm)
        }
      } else {
        aggregation_function<-function(func_data,func_weight) {
          colMeans(func_data,na.rm=na.rm)
        }
      }
    } else if (aggr_type=="weighted_mean"){
      aggregation_function<-function(func_data,func_weight) {
        if (crop_aggr) {
          if(any(is.na(func_weight))&(na.rm==TRUE)){ 
            func_weight[,,]<-apply(func_weight,c(2,3),function(x){replace(x, is.na(x), 0)})
          }
          # this is programmed very inefficient! sorry. BB
          retour<-array(NA,dim=c(dim(func_data)[2],1))
          for (column1 in 1:dim(func_data)[2]) {
            #for (column2 in 1:dim(func_data)[3]) {
              retour[column1,]<-weighted.mean(as.vector(func_data[,column1,]),as.vector(func_weight[,column1,]),na.rm=na.rm)
            #}
          }
          return(retour)
        } else {
          if(any(is.na(func_weight))&(na.rm==TRUE)){ 
            func_weight[,,]<-apply(func_weight,c(2,3),function(x){replace(x, is.na(x), 0)})
          }
          # this is programmed very inefficient! sorry. BB
          retour<-array(NA,dim=c(dim(func_data)[2],dim(func_data)[3]))
          for (column1 in 1:dim(func_data)[2]) {
            for (column2 in 1:dim(func_data)[3]) {
              retour[column1,column2]<-weighted.mean(as.vector(func_data[,column1,column2]),as.vector(func_weight[,column1,column2]),na.rm=na.rm)
            }
          }
          return(retour)
        }
      }  
    }else {
      stop("This aggregation type does not exist yet")
    }
  }

### Calculate on the aggregation of level

  if (level=="glo") {
    if (crop_aggr) {
      outdata<-array(NA, dim=c(1, nyears(data), 1),dimnames=list("GLO",getYears(data),NULL))
      if(!is.null(getSets(data))) names(dimnames(outdata))[2]<-getSets(data,fulldim=FALSE)[2]
      outdata[,,]<-aggregation_function(func_data=data[,,],func_weight=weight[,,])
      outdata<-as.magpie(outdata)
    } else {
      outdata<-array(NA, dim=c(1, nyears(data), ndata(data)),dimnames=list("GLO",getYears(data),getNames(data)))
      if(!is.null(getSets(data))) {names(dimnames(outdata))[2]<-getSets(data,fulldim=FALSE)[2] ; names(dimnames(outdata))[3]<-getSets(data,fulldim=FALSE)[3]}
      outdata[,,]<-aggregation_function(func_data=data[,,],func_weight=weight[,,])
      outdata<-as.magpie(outdata)
    }
  } else if (level=="reg") {
    if (crop_aggr) {
      outdata<-array(NA, dim=c(nregions(data), nyears(data), 1),dimnames=list(getRegions(data),getYears(data),NULL))
      if(!is.null(getSets(data))){ names(dimnames(outdata))[1] <- "i" ; names(dimnames(outdata))[2]<-getSets(data,fulldim=FALSE)[2]}
      for (region_x in getRegions(data)) {
        cells<-data[region_x,,]
        if(is.null(weight)==FALSE) {
          weight_cells<-weight[region_x,,]
        } else { weight_cells <- NULL }
        outdata[region_x,,]<-aggregation_function(func_data=cells,func_weight=weight_cells)
      }
      outdata<-outdata[order(dimnames(outdata)[[1]]),,,drop=FALSE]
      outdata<-as.magpie(outdata)
    } else {
      outdata<-array(NA, dim=c(nregions(data), nyears(data), ndata(data)),dimnames=list(getRegions(data),getYears(data),getNames(data)))
      if(!is.null(getSets(data))){ names(dimnames(outdata))[1] <- "i" ; names(dimnames(outdata))[2]<-getSets(data,fulldim=FALSE)[2] ; names(dimnames(outdata))[3] <-getSets(data,fulldim=FALSE)[3]}
      for (region_x in getRegions(data)) {
        cells<-data[region_x,,]
        if(is.null(weight)==FALSE) {
          weight_cells<-weight[region_x,,]
        } else { weight_cells <- NULL }
        outdata[region_x,,]<-aggregation_function(func_data=cells,func_weight=weight_cells)
      }
      outdata<-outdata[order(dimnames(outdata)[[1]]),,,drop=FALSE]
      outdata<-as.magpie(outdata)
    }
  } else if (level=="regglo") {
    if (crop_aggr) {
      #reg
      outdata<-array(NA, dim=c(nregions(data), nyears(data), 1),dimnames=list(getRegions(data),getYears(data),NULL))
      if(!is.null(getSets(data))){ names(dimnames(outdata))[1] <- "i" ; names(dimnames(outdata))[2]<-getSets(data,fulldim=FALSE)[2]}
      for (region_x in getRegions(data)) {
        cells<-data[region_x,,]
        if(is.null(weight)==FALSE) {
          weight_cells<-weight[region_x,,]
        } else { weight_cells <- NULL }
        outdata[region_x,,]<-aggregation_function(func_data=cells,func_weight=weight_cells)
      }
      outdata<-outdata[order(dimnames(outdata)[[1]]),,,drop=FALSE]
      reg<-as.magpie(outdata)
      #glo
      outdata<-array(NA, dim=c(1, nyears(data), 1),dimnames=list("GLO",getYears(data),NULL))
      if(!is.null(getSets(data))) names(dimnames(outdata))[2]<-getSets(data,fulldim=FALSE)[2]
      outdata[,,]<-aggregation_function(func_data=data[,,],func_weight=weight[,,])
      glo<-as.magpie(outdata)
      outdata <- mbind(reg,glo)
    } else {
      #reg
      outdata<-array(NA, dim=c(nregions(data), nyears(data), ndata(data)),dimnames=list(getRegions(data),getYears(data),getNames(data)))
      if(!is.null(getSets(data))){ names(dimnames(outdata))[1] <- "i" ; names(dimnames(outdata))[2]<-getSets(data,fulldim=FALSE)[2] ; names(dimnames(outdata))[3] <-getSets(data,fulldim=FALSE)[3]}
      for (region_x in getRegions(data)) {
        cells<-data[region_x,,]
        if(is.null(weight)==FALSE) {
          weight_cells<-weight[region_x,,]
        } else { weight_cells <- NULL }
        outdata[region_x,,]<-aggregation_function(func_data=cells,func_weight=weight_cells)
      }
      outdata<-outdata[order(dimnames(outdata)[[1]]),,,drop=FALSE]
      reg<-as.magpie(outdata)
      #glo
      outdata<-array(NA, dim=c(1, nyears(data), ndata(data)),dimnames=list("GLO",getYears(data),getNames(data)))
      if(!is.null(getSets(data))) {names(dimnames(outdata))[2]<-getSets(data,fulldim=FALSE)[2] ; names(dimnames(outdata))[3]<-getSets(data,fulldim=FALSE)[3]}
      outdata[,,]<-aggregation_function(func_data=data[,,],func_weight=weight[,,])
      glo<-as.magpie(outdata)
      outdata <- mbind(reg,glo)
    }
  } else if (level=="5regions") {
    five_regions<-c("O90","AFM","LAM","CPA","XAS")
    outdata<-new.magpie(five_regions,dimnames(data)[[2]],dimnames(data)[[3]])
    if(!is.null(getSets(data))){names(dimnames(outdata))[2]<-getSets(data,fulldim=FALSE)[2] ; names(dimnames(outdata))[3] <-getSets(data,fulldim=FALSE)[3]}
    alloc<-list(c("EUR","FSU","NAM","PAO"),c("AFR","MEA"),c("LAM"),"CPA",c("PAS","SAS"))
    
    for (region_x in 1:5) {
      cells<-data[alloc[region_x][[1]],,]
      if(is.null(weight)==FALSE) {
        weight_cells<-weight[alloc[region_x][[1]],,]
      } else { weight_cells <- NULL }
      outdata[region_x,,]<-aggregation_function(func_data=cells,func_weight=weight_cells)
    }
    outdata<-outdata[order(dimnames(outdata)[[1]]),,,drop=FALSE]
    outdata<-as.magpie(outdata)
    
  } else if (level=="SSP") {
    five_regions<-c("ASIA","LAM","MAF","OECD","REF")
    outdata<-new.magpie(five_regions,dimnames(data)[[2]],dimnames(data)[[3]])
    if(!is.null(getSets(data))){names(dimnames(outdata))[1]<-"" ; names(dimnames(outdata))[2]<-getSets(data,fulldim=FALSE)[2] ; names(dimnames(outdata))[3] <-getSets(data,fulldim=FALSE)[3]}
    alloc<-list(c("CPA","PAS","SAS"),c("LAM"),c("AFR","MEA"),c("EUR","NAM","PAO"),c("FSU"))
    
    for (region_x in 1:5) {
      cells<-data[alloc[region_x][[1]],,]
      if(is.null(weight)==FALSE) {
        weight_cells<-weight[alloc[region_x][[1]],,]
      } else { weight_cells <- NULL }
      outdata[region_x,,]<-aggregation_function(func_data=cells,func_weight=weight_cells)
    }
    outdata<-outdata[order(dimnames(outdata)[[1]]),,,drop=FALSE]
    outdata<-as.magpie(outdata)
    
  } else if (level=="country") {
    outdata<-array(NA, dim=c(length(getCountries()), nyears(data), ndata(data)),dimnames=list(getCountries(),getYears(data),getNames(data)))
    if(!is.null(getSets(data))){names(dimnames(outdata))[2]<-getSets(data,fulldim=FALSE)[2] ; names(dimnames(outdata))[3] <-getSets(data,fulldim=FALSE)[3]}
    names(dimnames(outdata))<-getSets(data)
    for (country in getCountries()) {
      cells<-data[countrycells(country),,]
      if(length(cells)>0){
        if(is.null(weight)==FALSE) {
          weight_cells<-weight[countrycells(country),,]
        } else { weight_cells <- NULL }
          outdata[country,,]<-aggregation_function(func_data=cells,func_weight=weight_cells)    
      } else {
        if(aggr_type=="sum") outdata[country,,] <- 0
      }
    }
    message("Object transformed to array")
  } else  { stop("level does not exist")}

return(outdata)
}


setMethod("superAggregate",
    signature(data = "lpj"),
    function (data, aggr_type, level = "reg", weight = NULL) 
    {
        if (is.null(weight)==FALSE) {
          if (identical(dim(weight),dim(data))==FALSE) {
            stop("data and weight have to have the same dimensions!")
          }  
        }
        retour<-list()
        for (four_d in 1:dim(data)[4]) {
          data_x<-as.magpie(data[,,,four_d])[[1]]
          if (is.null(weight)==FALSE) {
            weight_x<-as.magpie(weight[,,,four_d])[[1]]
          } else {weight_x<-NULL}
          retour[[four_d]]<-superAggregate(data_x, aggr_type, level = level, weight = weight_x)

        } 
        if (is.null(dimnames(data)[[4]])==FALSE) {names(retour)[[1]]<-dimnames(data)[[4]]}
        if(level=="country"){warning("LPJ object transformed into list of arrays")}
        if((level=="reg")|(level=="glo")){warning("LPJ object transformed into list of MAgPIE objects")}
        return(retour)
    }
)
