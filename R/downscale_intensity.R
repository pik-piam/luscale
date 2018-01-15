#' downscale intensity
#' 
#' Downscaling of intensities like GDP/cap, Energy/GDP, Emissions/Energy from
#' regional to country level. Intensities have the structure
#' dataNominator/dataDenominator.
#' 
#' downscale_intensity is implementing a convergence based downscaling
#' mechanism. This means that the intensity under investigation is assumed to
#' converge to a common value for all countries in a region at some point in
#' the future. Thus, the differences
#' 
#' @usage downscale_intensity(nomiR, denomR, nomibaseC, denomC, mapping,
#' convmethod="exp", convrate="standard", replacemissing=FALSE,
#' adjustnomiR=TRUE, adjustdenomR=TRUE, adjustshare="by nomiC")
#' @param nomiR the data in the nominator (e.g. GDP,Energy,Emissions) for all
#' regions for all relevant years.
#' @param denomR the data in the denominator (e.g. Population,GDP,Energy) for
#' all regions for all relevant years.
#' @param nomibaseC base year data of the nominator at country (or cell, if
#' applied as country to cell) level, e.g. GDP in 2010
#' @param denomC data for the denominator at country level for all relevant
#' years, e.g. population, GDP, Energy
#' @param mapping file for mapping countries to regions and vice versa.
#' @param convmethod the method of convergence of the intensities. Exponential
#' ("exp") or linear ("lin").
#' @param convrate A number between 0 and 1 which specifies the degree of
#' convergence in the last data year. 1 means total convergence. Standard is
#' 1-exp(-1) for exponential (i.e. the difference has decayed to e^(-1) of its
#' initial value) and 1 for linear.
#' @param replacemissing specifies whether unavailable data should be replaced
#' with region averages or not. If TRUE, downscale_intensity uses all available
#' information of the input data, i.e. it uses data from nomibaseC and denomC
#' even if the corresponding value of denomC resp. nomibaseC is not available
#' by using region averages. If FALSE, information is ignored when it is not
#' available for one of the two.Results will differ slightly depending on the
#' choice of replacemissing. This is because of the different information
#' content used.
#' @param adjustnomiR if replacemissing is FALSE, setting adjustnomiR to TRUE
#' will adjust the input data by substracting the data of those countries that
#' are ignored from their respective regions.
#' @param adjustdenomR same as adjustnomiR for the denominator.
#' @param adjustshare specifies how discrepancies between the downscaled data
#' and the regional Remind output are attributed to the countries. "by nomiC"
#' uses the denominator as weight (e.g. countries with higer GDP receive a
#' bigger share of the difference term). "by growth" uses the growth in every
#' timestep as weight, e.g. countries that have higher growth receive a bigger
#' share of the difference term.
#' @return \item{downscaled_intensity}{Magpie-object containing the downscaled
#' intensity for all available years and countries.}
#' @note The first year for which data is available is automatically taken as
#' the base year for the calculation of convergence. Thus, the data that is fed
#' in should always start at the intended base year.
#' @author Roman Julius Hennig
#' @export
#' @importFrom magclass getYears getYears<- getNames<- as.magpie
#' @seealso \code{\link{speed_aggregate}}
#' @examples
#' 
#' \dontrun{d1 <- downscale_intensity(gdpRegion,gdpCountry[,2010,],popCountry,mapping)}
#' 
#' \dontrun{d2 <- downscale_intensity(energyRegion,energyCountry[,2010,],gdpCountry,mapping)}
#' 
#convergence based downscaling function; only for weighted variables where convergence is assumed 
#(e.g. per capita GDP, energy per $GDP, emissions per unit energy)
#nomiR - data for regions that should be downscaled to country level
#nomibaseC - country level data in base year to initialize
#denomC - gdp projections for all years for all countries (e.g. from OECD input data)
#mapping - country to region mapping
downscale_intensity <- function(nomiR,denomR,nomibaseC,denomC,mapping,convmethod="exp",convrate="standard",replacemissing=FALSE,
                                adjustnomiR=TRUE,adjustdenomR=TRUE,adjustshare="by nomiC") {
  if (length(getYears(nomiR)) != length(getYears(denomC))) {warning("years for nomiR and denomC input data do not agree.
                                                                    Downscaling can only be done for years where both are available.")}
  yrs <- c()
  for(k in getYears(nomiR)) {
    if(sum(nomiR[,k,])>0 & sum(denomC[,k,])>0) {
      yrs <- c(yrs,as.numeric(sub("y","",k)))
    }
  }  
  l <- length(yrs) 
  ldy <- yrs[l]
  by <- yrs[1]
  t <- ldy - by
  map <- getAggregationMatrix(mapping)
  denomC <- speed_aggregate(denomR,mapping,weight=denomC)
  denomR <- speed_aggregate(denomC,mapping)
  nomibaseC <- speed_aggregate(nomiR[,2010,],mapping,weight=nomibaseC)
  if(!replacemissing) {
    if(adjustnomiR) {
      #summing up data where the respective data for denomC or nomiC is missing 
      temp1 <- nomibaseC
      temp1[denomC[,by,]!=0,,] <- 0
      temp2 <- speed_aggregate(temp1,mapping)/nomiR[,by,]
      getYears(temp2) <- NULL
      #adjusting the region level data (i.e. substracting the missing data from the region level data)
      nomiR <- (1-temp2)*nomiR
    }
    if(adjustdenomR) {
      temp3 <- denomC
      temp3[nomibaseC[,by,]!=0,,] <- 0   
      denomR <- denomR-speed_aggregate(temp3,mapping)
    }    
    #reconciling data by replacing missing values with 0 if not available
    nomibaseC[denomC[,by,]==0,,] <- 0
    denomC[nomibaseC[,by,]==0,,] <- 0 
  }  
  getNames(nomiR) <- NULL
  getNames(denomR) <- NULL
  intensityR <- nomiR/denomR
  dummyR <- denomC
  for(k in yrs) {dummyR[,k,] <- as.magpie(intensityR[,k,]%*%map)}
  initR <- dummyR[,by,]
  getYears(initR) <- NULL
  initC <- speed_aggregate(nomiR[,2010,],mapping,weight=nomibaseC)/denomC[,by,]
  #whenever an intensity is not available for a country in the base year for whatever reasons, it is replaced with the
  #respective intensity in the region where the country is located, in order to be able to calculate further results.
  if(replacemissing) {initC[is.nan(initC)|is.infinite(initC)|initC==0,,] <- initR[is.nan(initC)|is.infinite(initC)|initC==0,,]}
  getYears(initC) <- NULL
  initRel <- initR/initC
  intensityC <- denomC
  
  if(convmethod=="exp"){
    if(convrate=="standard"){convrate <- (1-exp(-1))}
    if(convrate>=1){stop("Full convergence not possible in exponential mode. Please specify a convrate < 1.")}
    x <- log(1-convrate)/t
    for(k in yrs) {intensityC[,k,] <- dummyR[,k,]/(1+exp(x*(k-by))*(initRel-1))} 
    }
  if(convmethod=="lin"){
    if(convrate=="standard"){convrate <- 1}
    if(convrate > 1){warning("Convrate > 1 means convergence will be reached before the last year of data availability. convrate should be </= 1.")}
    m <- convrate*(1-initRel)/t
    for(k in yrs) {intensityC[,k,] <- dummyR[,k,]/(initRel+(k-by)*m)}  
  }    
  #intensityC[is.nan(intensityC)|is.infinite(intensityC),,] <- dummyR[is.nan(intensityC)|is.infinite(intensityC),,]
  if(!replacemissing){for(k in yrs) {intensityC[nomibaseC==0,k,] <- 0}}

  #redistribute total data with new weights 
  dsnomiC <- intensityC*denomC
  sum <- speed_aggregate(dsnomiC,mapping)
  diff <- nomiR-sum
  
  if(adjustshare=="by growth") {
    #calculate share weights
    tempsum <- denomC      
    tempdiff <- denomC     #placeholders
    for(k in yrs) {
      tempsum[,k,] <- as.magpie(sum[,k,]%*%map)
      tempdiff[,k,] <- as.magpie(diff[,k,]%*%map)
    }
    share <- dsnomiC[,2:l,]
    for(k in 2:l) {
      temp2 <- dsnomiC[,yrs[k-1],]
      getYears(temp2) <- NULL
      temp3 <- tempsum[,yrs[k-1],]
      getYears(temp3) <- NULL
      share[,k-1,] <- (dsnomiC[,k,]-temp2)/(tempsum[,k,]-temp3)
      #share[,k-1,] <-       
    }    
    final <- intensityC
    final[,2:l,] <- intensityC[,2:l,] + tempdiff[,2:l,]*share/denomC[,2:l,]
  }
  
  if(adjustshare=="by nomiC") {
    final <- intensityC
    final[,2:l,] <- speed_aggregate(nomiR,mapping,weight=(final*denomC))[,2:l,]/denomC[,2:l,]  
  }
  if (replacemissing) {final[denomC[,by,]==0,,] <- dummyR[denomC[,by,]==0,,]}
  return(final)
}
