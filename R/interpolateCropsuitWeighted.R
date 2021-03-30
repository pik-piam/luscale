#' interpolateCropsuitWeighted
#'
#' Disaggregates a modelled time series of land pools after optimisation from the model resolution to the resolution
#' of the high resolution land initialisation data set, based on a spam matrix and cropland suitability information.
#'
#' The function requires the following input data:
#' \itemize{
#' \item \code{x} is an object containing a time series of land pools (model output). The sum over all land pools
#' is constant over time.
#' \item \code{x_ini_lr} and \code{x_ini_hr} provide the initial land pools (Mha) at high (hr) and low resolution (lr)
#' before the optimisation. They only contain the initial time step, but share the three-dimensional structure
#' with \code{x}.
#' \item \code{cropsuit_lr} and \code{cropsuit_hr} provide information about the amount (Mha) of suitable cropland
#' at high and low resolutions.
#' \item \code{spam} spam matrix containing information about cell belongings between the high and low resolution.
#' }
#'
#' The weighted disaggregation works as follows:
#'
#' 1. From the time series of land pools, the total area (Mha) af cropland expansion and reduction  at low resolution
#' is calculated for each time step.
#'
#' 2. In order to allocate the area of cropland expansion and reduction, for each grid cell at high resolution
#' expansion and reduction weights are calculated:
#' \itemize{
#' \item The expansion weight is calculated as the ratio between available cropland at the grid cell level
#' (high resolution) and total available cropland at the low resolution spatial unit. Available cropland is
#' the difference between suitable cropland and the amount of cropland at the previous time step minus urban land,
#' since it assumed that cropland cannot be allocated to urban land.
#' \item The reduction weight is given by the ratio between the amount of cropland per grid cell and the total area
#' of cropland at the low resolution spatial unit. This assumes that the cropland reduction is equally distributed
#' among all high resolution grid cells in terms of the amount of cropland within each cell.
#' }
#'
#' 3. The total area (Mha) of cropland expansion and reduction in the current time step at low resolution is
#' multiplied with the allocation weights calculated at the previous time step (see step 2.) to distribute cropland
#' at the high resolution. This procedure is repeated for all time steps.
#'
#' 4. Following the cropland allocation, the land area for the remaining non-cropland vegetation pools is calculated
#' by substracting the allocated cropland and urban land areas from the total land area in each grid cell.
#'
#' 5. The non-cropland vegetation pool at the high resolution, calculated in step 4., is then multiplied by the
#' respective shares of the remaining non-cropland vegetation pools at the previous time step (temporary allocation). 
#' This, however, is not sufficient to also account for changes within these land pools. Therefore, the residual
#' changes in these land poolsare calculated by comparing the sum of the temporarily allocated non-cropland pools with
#' the pools at low resolution. The residual area of land expansion and reduction is then allocated by based on 
#' reduction and expansion weights, similar as in 2.. The reduction weight is calculated as the ratio between the 
#' given temporary land pool at high resolution and total temporary land pool at low resolution. The expansion weight
#' is calculated as the ratio between the remaining land to be filled in each land pool and the total amount of 
#' residual land to be allocated in the current time step. Urban land is assumed to be constant over time.
#'
#' @param x Time series of land pools (model output) to be disaggregated.
#' @param x_ini_lr The low resolution distribution of x before model optimization.
#' @param x_ini_hr The initial high resolution distribution of x (land-use initialisation) before model optimization.
#' @param cropsuit_lr The area of land suitable for cropland at the low resolution.
#' @param cropsuit_hr The area of suitable cropland at the high resolution.
#' @param spam The spam file (sum) that is to be used for disaggregation.
#' @param cropland_scen different options are
#' \itemize{
#' \item \code{"allmarginal_0pNonCropVeg"}: All marginal land, 0 \% conservation of non-cropland vegetation
#' \item \code{"allmarginal_10pNonCropVeg"}: All marginal land, 10 \% conservation of non-cropland vegetation
#' \item \code{"allmarginal_20pNonCropVeg"}: All marginal land, 20 \% conservation of non-cropland vegetation
#' \item \code{"allmarginal_30pNonCropVeg"}: All marginal land, 30 \% conservation of non-cropland vegetation
#' \item \code{"halfmarginal_0pNonCropVeg"}: Half of the marginal land excluded, 0 \% conservation of non-cropland vegetation
#' \item \code{"halfmarginal_10pNonCropVeg"}: Half of the marginal land excluded, 10 \% conservation of non-cropland vegetation
#' \item \code{"halfmarginal_20pNonCropVeg"}: Half of the marginal land excluded, 20 \% conservation of non-cropland vegetation
#' \item \code{"halfmarginal_30pNonCropVeg"}: Half of the marginal land excluded, 30 \% conservation of non-cropland vegetation
#' \item \code{"nomarginal_0pNonCropVeg"}: No marginal land, 0 \% conservation of non-cropland vegetation
#' \item \code{"nomarginal_10pNonCropVeg"}: No marginal land, 10 \% conservation of non-cropland vegetation
#' \item \code{"nomarginal_20pNonCropVeg"}: No marginal land, 20 \% conservation of non-cropland vegetation
#' \item \code{"nomarginal_30pNonCropVeg"}: No marginal land, 30 \% conservation of non-cropland vegetation
#' }
#' @param year_ini Timestep that is assumed for the initial distributions \code{x_ini_hr} and \code{x_ini_lr}.
#' @param unit Unit of the output. "Mha" or "share"
#' @return The disaggregated MAgPIE object containing x_ini_hr as first
#' timestep
#' @export
#' @author Patrick von Jeetze
#' @importFrom magclass is.magpie nregions nyears getNames getYears mbind dimSums setYears getYears new.magpie
#' @seealso \code{\link{interpolate}}
#' \code{\link{speed_aggregate}}
#' @examples
#' \dontrun{
#' a <- interpolateCropsuitWeighted(
#'   x = land,
#'   x_ini_lr = land_ini_lr,
#'   x_ini_hr = land_ini_hr,
#'   cropsuit_lr = "avl_cropland.cs3",
#'   cropsuit_hr = "avl_cropland_0.5.mz",
#'   cropland_scen="allmarginal_0pNonCropVeg",
#'   spam = "0.5-to-c200_sum.spam"
#' )
#' }
#'
interpolateCropsuitWeighted <- function(x, x_ini_lr, x_ini_hr, cropsuit_lr, cropsuit_hr, cropland_scen="allmarginal_0pNonCropVeg", spam, year_ini = "y1985", unit = "Mha") {

  # test whether data can be handled by function
  if (!is.magpie(x) || !is.magpie(x_ini_lr) || !is.magpie(x_ini_hr)) stop("x, x_ini_lr and x_ini_hr have to be magpie objects")
  if (nregions(x) != nregions(x_ini_lr)) stop("x and x_ini_lr have to be of the same spatial aggregation")
  if (nyears(x_ini_lr) > 1 || nyears(x_ini_hr) > 1) stop("Initialization data must only have one timestep")
  if (!all(getNames(x) == getNames(x_ini_lr)) || !all(getNames(x) == getNames(x_ini_hr))) stop("dimnames[[3]] of x, x_ini_lr and x_ini_hr have to be the same")

  # get land-use intialisation year
  getYears(x_ini_hr) <- year_ini
  getYears(x_ini_lr) <- year_ini

  # create data set with magpie output data (x) and land-use initialisation (x_ini_lr)
  lr <- mbind(x_ini_lr, x)

  # Test if the total sum is constant
  test <- dimSums(lr, dim = c(1, 3))[, 2:nyears(lr), ] - setYears(dimSums(lr, dim = c(1, 3))[, 1:nyears(lr) - 1, ], getYears(lr)[2:nyears(lr)])
  if (max(test) > 0.1 || min(test) < -0.1) warning("Total stock is not constant over time. See help for details")

  ### calculate land expansion and reduction

  # difference in land pools over time (Mha)
  land_diff_lr <- lr[, 2:nyears(lr), ] - setYears(lr[, 1:(nyears(lr) - 1), ], getYears(lr)[2:nyears(lr)])
  # land reduction
  land_reduc_lr <- land_diff_lr
  land_reduc_lr[land_reduc_lr > 0] <- 0
  land_reduc_lr <- abs(land_reduc_lr)
  # land expansion
  land_expan_lr <- land_diff_lr
  land_expan_lr[land_expan_lr < 0] <- 0

  # calculate shares of non crop land pools in each time step
  nocrop_pools <- getNames(x_ini_lr)[-which(getNames(x_ini_lr) == "crop")]
  nocrop_veg <- nocrop_pools[-which(nocrop_pools == "urban")]

  # read crop suitabiliy data
  if (is.character(cropsuit_hr)) {
    if (!file.exists(cropsuit_hr)) stop("cropland suitability file not found")
    # low resolution
    cropland_suit_lr <- read.magpie(cropsuit_lr)
    # high resolution
    cropland_suit_hr <- read.magpie(cropsuit_hr)

    # correct for urban land because it is constant
    # where suitable cropland is larger than total non urban land chose the smaller value [pmin()]
    # high resolution
    land_non_urban_hr <- (dimSums(x_ini_hr, dim = 3) - x_ini_hr[, , "urban"])
    getCells(cropland_suit_hr) <- getCells(land_non_urban_hr)
    cropland_suit_hr <- pmin(cropland_suit_hr[, , cropland_scen], land_non_urban_hr)
    # low resolution
    land_non_urban_lr <- (dimSums(x_ini_lr, dim = 3) - x_ini_lr[, , "urban"])
    cropland_suit_lr <- pmin(cropland_suit_lr[, , cropland_scen], land_non_urban_lr)
  }

  # compute available cropland for each time step
  cropland_avail_lr <- cropland_suit_lr[, , cropland_scen] - lr[, , "crop"]
  getNames(cropland_avail_lr) <- "crop"

  # get spam file
  if (is.character(spam)) {
    if (!file.exists(spam)) stop("spam file ", spam, " not found")
    rel <- read.spam(spam)
  }

  ### disaggregate low resolution output data

  # total land expansion and reduction at low resolution
  land_reduc_lr_dagg <- speed_aggregate(land_reduc_lr, t(rel))
  land_expan_lr_dagg <- speed_aggregate(land_expan_lr, t(rel))
  # available cropland at low resolution
  cropland_avail_lr_dagg <- speed_aggregate(cropland_avail_lr, t(rel))
  # total amount of cropland at low resolution
  land_lr_dagg <- speed_aggregate(lr, t(rel))

  ########### allocate land pools in high res (hr) ###########################

  ### allocate cropland

  # create new magpie object for cropland allocation (high resolution)
  cropland_hr <- new.magpie(getCells(land_expan_lr_dagg), c(year_ini, getYears(land_expan_lr_dagg)), "crop")
  getCells(x_ini_hr) <- getCells(land_expan_lr_dagg)
  cropland_hr[, year_ini, ] <- x_ini_hr[, , "crop"]

  for (t in 2:nyears(cropland_hr)) {

    # cropland in previous time step
    cropland <- cropland_hr[, t - 1, ]
    # cropland expansion
    cropland_expan <- land_expan_lr_dagg[, t - 1, "crop"]
    # cropland reduction
    cropland_reduc <- land_reduc_lr_dagg[, t - 1, "crop"]
    # available cropland in each grid cell (high resolution)
    cropland_avail <- setNames(cropland_suit_hr[, , cropland_scen] - cropland_hr[, t - 1, ], "crop")
    cropland_avail[cropland_avail < 0, , ] <- 0
    # expansion weight: divide available cropland in each grid cell
    # by total available cropland of the respective cluster (low resolution)
    cropland_expan_weight <- cropland_avail / cropland_avail_lr_dagg[, t - 1, ]
    cropland_expan_weight[is.na(cropland_expan_weight), , ] <- 0
    # reduction weight: divide cropland in each grid cell by the
    # total amount of cropland in each cluster (low resolution)
    cropland_reduc_weight <- cropland_hr[, t - 1, ] / land_lr_dagg[, t - 1, "crop"]
    cropland_reduc_weight[is.na(cropland_reduc_weight) | is.infinite(cropland_reduc_weight), , ] <- 0

    # disaggregate cropland for each time step
    cropland_hr[, t, "crop"] <- cropland + cropland_expan * cropland_expan_weight - cropland_reduc * cropland_reduc_weight
  }

  ### allocate non-cropland pools

  # calculate total land area in each grid cell (high resolution)
  land_tot_hr <- dimSums(x_ini_hr, dim = 3)
  # calculate non-cropland vegetation pool after disaggregation
  land_tot_nocrop_veg_hr <- setNames(setYears(land_tot_hr - x_ini_hr[, , "urban"] - cropland_hr, getYears(cropland_hr)), NULL)
  land_tot_nocrop_veg_hr[land_tot_nocrop_veg_hr < 0] <- 0
  
  # create new magpie object
  land_nocrop_veg_hr <- new.magpie(getCells(cropland_hr), getYears(cropland_hr), nocrop_veg)
  land_nocrop_veg_hr[, year_ini, nocrop_veg] <- x_ini_hr[, , nocrop_veg]


  for (t in 2:nyears(land_nocrop_veg_hr)) {
    
    # calculate share of non-cropland vegetation pools in previous time step
    shr_prev_nocrop_veg_hr <- land_nocrop_veg_hr[, t - 1, ] / (land_tot_nocrop_veg_hr[, t - 1, ] + 1e-100)
    shr_prev_nocrop_veg_hr[is.na(shr_prev_nocrop_veg_hr) | is.infinite(shr_prev_nocrop_veg_hr)] <- 0
    # multiply shares of non-cropland pools in previous time step with available land in current time step
    land_nocrop_veg_hr[, t, ] <- shr_prev_nocrop_veg_hr * land_tot_nocrop_veg_hr[, t, ]

    # sum temporary non-cropland land pool at low resolution to compare them with the land pools in x 
    tmp_land_nocrop_lr <- speed_aggregate(land_nocrop_veg_hr[, t, ], rel)
    tmp_land_nocrop_lr_dagg <- speed_aggregate(tmp_land_nocrop_lr, t(rel))

    # calculate the residual difference that still needs to be allocated 
    residual_diff_lr <- land_lr_dagg[, t, nocrop_veg] - tmp_land_nocrop_lr_dagg
    # calculate the residual land reduction (absolute)
    residual_reduc_lr <- residual_diff_lr
    residual_reduc_lr[residual_reduc_lr > 0] <- 0
    residual_reduc_lr <- abs(residual_reduc_lr)
    
    # calculate the allocation weight for the residual land reduction:
    # divide temporay non-croplad pools at high resolution by sum of the temporary at low resolution
    nocrop_reduc_weight <- land_nocrop_veg_hr[, t, ] / (tmp_land_nocrop_lr_dagg )
    nocrop_reduc_weight[is.na(nocrop_reduc_weight) | is.infinite(nocrop_reduc_weight)] <- 0
    nocrop_reduc_weight[nocrop_reduc_weight<0] <- 0 ; nocrop_reduc_weight[nocrop_reduc_weight>1] <- 1
    
    # allocate the residual non-cropland pool reduction 
    land_nocrop_veg_hr[, t, ] <- land_nocrop_veg_hr[, t, ] - residual_reduc_lr * nocrop_reduc_weight
    
    # calculate the residual land expansion
    residual_expan_lr <- residual_diff_lr
    residual_expan_lr[residual_expan_lr < 0] <- 0 
    
    #calculate the remaining available land at high resolution
    nocrop_avail_hr <- (land_tot_nocrop_veg_hr[, t, ] - dimSums(land_nocrop_veg_hr[, t, ], dim=3))
    nocrop_avail_hr[nocrop_avail_hr < 0] <- 0
    
    # calculate the expansion weight for the residual expansion
    # divide the available land by the total area of the residual expansion at low resolution
    nocrop_expan_weight <- nocrop_avail_hr/(residual_expan_lr)
    nocrop_expan_weight[is.na(nocrop_expan_weight) | is.infinite(nocrop_expan_weight)] <- 0
    nocrop_expan_weight[nocrop_expan_weight<0] <- 0 ; nocrop_expan_weight[nocrop_reduc_weight>1] <- 1
    
    # allocate the residual non-cropland pool expansion 
    land_nocrop_veg_hr[, t, ] <- land_nocrop_veg_hr[, t, ] + residual_expan_lr*nocrop_expan_weight
  }
  
  # calculate disaggregated urban land pool
  land_urban_hr <- setYears(setNames(land_tot_hr - (cropland_hr + dimSums(land_nocrop_veg_hr, dim = 3)), "urban"), getYears(cropland_hr))

  # combine output
  hr <- mbind(cropland_hr, land_nocrop_veg_hr, land_urban_hr)

  if (unit == "share") {
    # divide land pools by total amount of land per grid cell
    out <- hr / land_tot_hr
  } else if (unit == "Mha") {
    out <- hr
  }

  return(out)
}


