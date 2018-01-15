#' Interpolate
#' 
#' Disaggregates a cellular MAgPIE output to 0.5 degree based on a spam matrix
#' and information about the initial 0.5 degree distribution.
#' 
#' The function is based on the following assumption: \code{x} is an object
#' with more than one data dimension and the sum over the (normally two) data
#' dimensions is constant over time. Example: One column cropland, the other
#' one (cell size - cropland). \code{x_ini_hr} and \code{x_ini_lr} have to be
#' of the same structure (except for the time dimension of course). The
#' function calculates the amount by which the individual data columns of x
#' change in each timestep. The output is based on \code{x_ini_hr} and only the
#' differences in later timesteps to ths starting point are disaggregated by
#' the spam matrix. This assures that as little information as possible is lost
#' from the original dataset \code{x_ini_hr}.
#' 
#' The disaggregation procedure itself works as follows:
#' 
#' 1. Differences in distribution between years are derived for the low
#' resolution data set.
#' 
#' 2. Based on these differences extension and reduction shares are calculated
#' for the different pools. Reduction shares are calculated relative to the
#' pool itself (e.g. a reduction in a cropland pool from 10ha to 6ha leads to a
#' reduction share of (10ha-6ha)/10ha = 40\%). At the same time extension
#' shares are calculated relative to the pool which was made available by
#' reductions of the other pools (e.g. cropland is reduced from 10ha to 6ha,
#' forest area is reduced from 2ha to 1ha, but pastureland increases from 20ha
#' to 25ha. In this case the extension share of pasture will be
#' (25ha-20ha)/(10ha-6ha+2ha-1ha)=5ha/5ha=100\%). This difference in
#' calculation of reduction and extension share is crucial for the application
#' at the high resolution level because otherwise the calculation will not add
#' up.
#' 
#' 3. Reduction and extension shares are disaggregated to the high resolution
#' level by just assigning the same low resolution shares to all belonging
#' cells at the higher resolution.
#' 
#' 4.Starting with the provided high resolution pool data set for the initial
#' year reduction shares are applied on all pools in all cells. The pool which
#' is made available for expansions is calculated by summing up all values
#' which were released by the pool reductions.
#' 
#' 5. Pool expansions are calculated based on the pool made available in 4 for
#' the first time step.
#' 
#' 6. Steps 4 and 5 are repeated for all the following years based on the newly
#' created high resolution data.
#' 
#' Applying this procedure makes sure that relative pool reductions are
#' identical for the low resolution cell and all belonging high resolutions
#' cells whereas the extension shares relative to the areas made available per
#' cell are identical between low resolution cell and belonging high resolution
#' cells.
#' 
#' @usage interpolate(x, x_ini_lr, x_ini_hr, spam, add_avail_hr=NULL,
#' prev_year="y1985")
#' @param x The object to be disaggregated. See details for further important
#' information
#' @param x_ini_lr The low resolution distribution of x before optimization
#' (e.g. obtainable from input/cellular.)
#' @param x_ini_hr The initial 0.5 degree distribution of x before optimization
#' (e.g. obtainable in the raw_data folder in the 0.5set folder of the
#' inputdata).
#' @param spam The spam matrix that is to be used for disaggregation.
#' @param add_avail_hr This argument is deprecated and can't be used anymore.
#' @param prev_year Timestep that is assumed for the initial distributions
#' x_ini_hr and x_ini_lr.
#' @return The disaggregated MAgPIE object containing x_ini_hr as first
#' timestep
#' @export
#' @author Jan Philipp Dietrich,Markus Bonsch
#' @importFrom magclass is.magpie nregions nyears getNames getYears<- mbind dimSums setYears getYears as.magpie new.magpie as.array
#' @seealso \code{\link{reshape_folder}},\code{\link{reshape_file}},
#' \code{\link{speed_aggregate}}
#' @examples
#' 
#'  \dontrun{a <- interpolate(x = land, 
#'                            x_ini_lr = land_ini_lr,
#'                            x_ini_hr = land_ini_hr,
#'                            spam = "0.5-to-n500.sum.spam")}
#' 
interpolate<-function(x,x_ini_lr,x_ini_hr,spam,add_avail_hr=NULL,prev_year="y1985"){
  if(!is.magpie(x) || !is.magpie(x_ini_lr)|| !is.magpie(x_ini_hr)) stop("x, x_ini_lr and x_ini_hr have to be magpie objects")
  if(nregions(x)!=nregions(x_ini_lr)) stop("x and x_ini_lr have to be of the same spatial aggregation")
  if(nyears(x_ini_lr)>1 || nyears(x_ini_hr)>1) stop("Initialization data must only have one timestep")
  if(!all(getNames(x)==getNames(x_ini_lr))||!all(getNames(x)==getNames(x_ini_hr))) stop("dimnames[[3]] of x, x_ini_lr and x_ini_hr have to be the same")
  if(!is.null(add_avail_hr)){ 
    stop("The add_avail functionality is deprecated and can't be used anymore")
  }
  
  getYears(x_ini_hr) <- prev_year
  getYears(x_ini_lr) <- prev_year
  lr<-mbind(x_ini_lr,x)
  #Test if the total sum is constant
  if(is.null(add_avail_hr)){
    test<-dimSums(lr,dim=c(1,3.1))[,2:nyears(lr),]- setYears(dimSums(lr,dim=c(1,3.1))[,1:nyears(lr)-1,],getYears(lr)[2:nyears(lr)])
    if(max(test)>0.1||min(test)< -0.1) warning("Total stock is not constant over time. See help for details")
  }
  #calculate reduction and extension shares which then can be disaggregated
  diff <- as.array(lr[,2:nyears(lr),]-setYears(lr[,1:(nyears(lr)-1),],getYears(lr)[2:nyears(lr)]))
  less <- diff; less[less>0] <- 0
  more <- diff; more[more<0] <- 0
  reduct <- -less/(lr[,1:(nyears(lr)-1),]+10^-100)
  avail <- rowSums(more,dims=2)
  extent <- as.magpie(more)
  for(e in getNames(extent)) extent[,,e] <- more[,,e]/(avail+10^-100)
  #disaggregate shares
  if(is.character(spam)){
    if(!file.exists(spam))stop("spam file ",spam," not found")
    rel <- read.spam(spam)
  }
  reduct_hr <- speed_aggregate(as.magpie(reduct),t(rel))
  extent_hr <- speed_aggregate(as.magpie(extent),t(rel))

  #calculate land pools in high res (hr)
  hr <- new.magpie(dimnames(reduct_hr)[[1]],c(prev_year,dimnames(reduct_hr)[[2]]),dimnames(reduct_hr)[[3]])
  if(is.null(add_avail_hr)){
    add_avail_hr<-array(0,dim=c(dim(reduct_hr)[1:2],1),dimnames=list(dimnames(reduct_hr)[[1]],dimnames(reduct_hr)[[2]],"add_avail_hr"))
  }
  add_avail_hr<-as.array(add_avail_hr)
  dimnames(x_ini_hr)[[1]] <- dimnames(reduct_hr)[[1]]
  hr[,prev_year,] <- x_ini_hr
  
  for(y in 2:nyears(hr)) hr[,y,] <- (1-reduct_hr[,y-1,])*setYears(hr[,y-1,],getYears(reduct_hr)[y-1]) + (rowSums(reduct_hr[,y-1,]*setYears(hr[,y-1,],getYears(reduct_hr)[y-1]))+add_avail_hr[,y-1,])*extent_hr[,y-1,]
  return(hr)
}