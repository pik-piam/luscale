#' interpolateAvlCroplandWeighted
#'
#' Disaggregates a modelled time series of land pools after optimisation from the model resolution (low resolution)
#' to the resolution of the land initialisation data set (high resolution), based on a spam matrix and available
#' cropland.
#'
#' The function requires the following input data:
#' \itemize{
#' \item \code{x} is an object containing a time series of land pools (model output). The sum over all land pools
#' is constant over time.
#' \item \code{x_ini_lr} and \code{x_ini_hr} provide the initial land pools (Mha) at high (hr) and low resolution (lr)
#' before the optimisation. They only contain the initial time step, but share the three-dimensional structure
#' with \code{x}.
#' \item \code{avl_cropland_lr} and \code{avl_cropland_hr} provide information about the amount (Mha) of available cropland
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
#' \item The expansion weight is calculated as the ratio between the remaining cropland at the grid cell level
#' (high resolution) and overall remaining cropland at the low resolution spatial unit. Remaining cropland is
#' the difference between the available cropland and the amount of cropland at the previous time step minus urban land,
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
#' 5. The non-cropland vegetation pool at the high resolution (except of primary forest), calculated in step 4., is
#' then multiplied by the respective shares of the remaining non-cropland vegetation pools at the previous time step
#' (temporary allocation). This, however, is not sufficient to also account for changes within these land pools.
#' Therefore, the residual changes in these land pools are calculated by comparing the sum of the temporarily
#' allocated non-cropland pools with the pools at low resolution. The residual area of land expansion and reduction
#' is then allocated by based on reduction and expansion weights, similar as in 2.. The reduction weight is calculated
#' as the ratio between the given temporary land pool at high resolution and total temporary land pool at low
#' resolution. The expansion weight is calculated as the ratio between the remaining land to be filled in each land
#' pool and the total amount of residual land to be allocated in the current time step.
#'
#' 6. Primary forest is treated in a slightly different way, as primary forest cannot be expanded over time. In
#' cropland cells with no cropland expansion, primary forest is, at first, assumed to remain constant and transferred
#' from the previous time step to the current time step. Once again, the sum of the temporarary allocation is compared
#' to the sum of primary forest at low resolution to determine the residual primary forest land, which still needs to
#' be allocated. Where there is an surplus of primary forest, the reduction weight is calculated similarly as in 5.,
#' the land area is reduced accordingly. In areas where the temprorily allocated primary forest falls short, the
#' allocation weight is calculated as a function of the difference in primary land between the previous time step and
#' in the current time step. This makes sure that there is no expansion of primary forest.
#'
#' 7. Urban land is assumed to be constant over time.
#'
#' @param x Time series of land pools (model output) to be disaggregated.
#' @param x_ini_lr The low resolution distribution of x before model optimization.
#' @param x_ini_hr The initial high resolution distribution of x (land-use initialisation) before model optimization.
#' @param avl_cropland_lr The area of land available for crop cultivation at the low resolution.
#' @param avl_cropland_hr The area of available cropland at the high resolution.
#' @param spam The spam file (sum) that is to be used for disaggregation.
#' @param cropland_scen Depending on the cropland suitability data, standard options are
#' \itemize{
#' \item \code{"all_marginal"}: Cropland can be allocated to marginal land
#' \item \code{"half_marginal"}: Half of the marginal land is excluded from cropland allocation
#' \item \code{"no_marginal"}: Marginal land is fully excluded from cropland
#' }
#' @param set_aside_shr Share of available cropland that is witheld for other land cover types
#' @param set_aside_fader Fader for share of set aside policy
#' @param year_ini Timestep that is assumed for the initial distributions \code{x_ini_hr} and \code{x_ini_lr}.
#' @param unit Unit of the output. "Mha" or "share"
#' @return The disaggregated MAgPIE object containing x_ini_hr as first
#' timestep
#' @export
#' @author Patrick von Jeetze
#' @importFrom magclass is.magpie nregions nyears getNames getYears mbind dimSums setYears getYears new.magpie where
#' @seealso \code{\link{interpolate}}
#' \code{\link{speed_aggregate}}
#' @examples
#' \dontrun{
#' a <- interpolateAvlCroplandWeighted(
#'   x = land,
#'   x_ini_lr = land_ini_lr,
#'   x_ini_hr = land_ini_hr,
#'   avl_cropland_lr = "avl_cropland.cs3",
#'   avl_cropland_hr = "avl_cropland_0.5.mz",
#'   spam = "0.5-to-c200_sum.spam",
#'   cropland_scen = "all_marginal"
#' )
#' }
#'
interpolateAvlCroplandWeighted <- function(x, x_ini_lr, x_ini_hr, avl_cropland_lr, avl_cropland_hr, spam,
                                           cropland_scen = "all_marginal", set_aside_shr = 0,
                                           set_aside_fader = NULL, year_ini = "y1985", unit = "Mha") {

  # test whether data can be handled by function
  if (!is.magpie(x) || !is.magpie(x_ini_lr) || !is.magpie(x_ini_hr)) stop("x, x_ini_lr and x_ini_hr have to be magpie objects")
  if (nregions(x) != nregions(x_ini_lr)) stop("x and x_ini_lr have to be of the same spatial aggregation")
  if (nyears(x_ini_lr) > 1 || nyears(x_ini_hr) > 1) stop("Initialization data must only have one timestep")
  if (!all(getNames(x) == getNames(x_ini_lr)) || !all(getNames(x) == getNames(x_ini_hr))) stop("dimnames[[3]] of x, x_ini_lr and x_ini_hr have to be the same")

  # ========================================================================
  # prepare data for land allocation
  # ========================================================================

  # get land-use intialisation year
  getYears(x_ini_hr) <- year_ini
  getYears(x_ini_lr) <- year_ini

  # create data set with magpie output data (x) and land-use initialisation (x_ini_lr)
  lr <- mbind(x_ini_lr, x)

  # Test if the total sum is constant
  test <- dimSums(lr, dim = c(1, 3))[, 2:nyears(lr), ] - setYears(dimSums(lr, dim = c(1, 3))[, 1:nyears(lr) - 1, ], getYears(lr)[2:nyears(lr)])
  if (max(test) > 0.1 || min(test) < -0.1) warning("Total stock is not constant over time. See help for details")

  #------------------------------------------------------------------------
  # calculate land expansion and reduction
  #------------------------------------------------------------------------

  # difference in land pools over time (Mha)
  land_diff_lr <- lr[, 2:nyears(lr), ] - setYears(lr[, 1:(nyears(lr) - 1), ], getYears(lr)[2:nyears(lr)])
  # land reduction
  land_reduc_lr <- land_diff_lr
  land_reduc_lr[land_reduc_lr > 0] <- 0
  land_reduc_lr <- abs(land_reduc_lr)
  # land expansion
  land_expan_lr <- land_diff_lr
  land_expan_lr[land_expan_lr < 0] <- 0

  #------------------------------------------------------------------------
  # read crop suitabiliy data
  #------------------------------------------------------------------------

  if (is.character(avl_cropland_hr) | is.character(avl_cropland_lr)) {
    if (!file.exists(avl_cropland_hr)) stop("high resolution available cropland data not found")
    if (!file.exists(avl_cropland_lr)) stop("low resolution available cropland file not found")

    # low resolution
    avl_cropland_lr <- read.magpie(avl_cropland_lr)
    # high resolution
    avl_cropland_hr <- read.magpie(avl_cropland_hr)

    # correct for urban land because it is constant
    # where available cropland is larger than total non urban land chose the smaller value [pmin()]
    # high resolution
    land_non_urban_hr <- (dimSums(x_ini_hr, dim = 3) - x_ini_hr[, , "urban"])
    getCells(land_non_urban_hr) <- getCells(avl_cropland_hr)
    avl_cropland_hr <- pmin(avl_cropland_hr[, , cropland_scen], land_non_urban_hr)
    # low resolution
    land_non_urban_lr <- (dimSums(x_ini_lr, dim = 3) - x_ini_lr[, , "urban"])
    avl_cropland_lr <- pmin(avl_cropland_lr[, , cropland_scen], land_non_urban_lr)
  }

  # expand available cropland data over time
  # high resolution
  avl_cropland_hr_tmp <- new.magpie(getCells(avl_cropland_hr), getYears(lr), cropland_scen)
  avl_cropland_hr_tmp[, getYears(lr), ] <- avl_cropland_hr
  avl_cropland_hr <- avl_cropland_hr_tmp
  # low resolution
  avl_cropland_lr_tmp <- new.magpie(getCells(avl_cropland_lr), getYears(lr), cropland_scen)
  avl_cropland_lr_tmp[, getYears(lr), ] <- avl_cropland_lr
  avl_cropland_lr <- avl_cropland_lr_tmp


  if (set_aside_shr != 0 & is.null(set_aside_fader)) stop("Share of withheld cropland given, but no policy fader for target year provided")
  if (set_aside_shr != 0 & !is.null(set_aside_fader)) {
    # correct available cropland with policy restriction
    for (t in 1:nyears(lr)) {
      # high resolution
      avl_cropland_hr[, t, ] <- avl_cropland_hr[, t, ] * (1 - set_aside_shr * set_aside_fader[getYears(lr)[t], ])
      # low resolution
      avl_cropland_lr[, t, ] <- avl_cropland_lr[, t, ] * (1 - set_aside_shr * set_aside_fader[getYears(lr)[t], ])
    }
  }

  # compute available cropland for expansion in each time step
  cropland_remain_lr <- avl_cropland_lr[, , cropland_scen] - lr[, , "crop"]
  getNames(cropland_remain_lr) <- "crop"
  getYears(cropland_remain_lr) <- getYears(lr)

  # get spam file
  if (is.character(spam)) {
    if (!file.exists(spam)) stop("spam file ", spam, " not found")
    rel <- read.spam(spam)
  }

  #------------------------------------------------------------------------
  # disaggregate low resolution output data to be used in calculations
  #------------------------------------------------------------------------

  # total land expansion and reduction at low resolution
  land_reduc_lr_dagg <- speed_aggregate(land_reduc_lr, t(rel))
  land_expan_lr_dagg <- speed_aggregate(land_expan_lr, t(rel))
  # available cropland at low resolution
  cropland_remain_lr_dagg <- speed_aggregate(cropland_remain_lr, t(rel))
  # total amount of cropland at low resolution
  land_lr_dagg <- speed_aggregate(lr, t(rel))

  # ========================================================================
  # allocate cropland at high res (hr)
  # ========================================================================

  # create new magpie object for cropland allocation (high resolution)
  cropland_hr <- new.magpie(getCells(land_expan_lr_dagg), getYears(lr), "crop")
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
    cropland_remain <- setNames(avl_cropland_hr[, t, cropland_scen] - cropland_hr[, t - 1, ], "crop")
    cropland_remain[cropland_remain < 0, , ] <- 0
    # expansion weight: divide available cropland in each grid cell
    # by total available cropland of the respective cluster (low resolution)
    cropland_expan_weight <- cropland_remain / cropland_remain_lr_dagg[, t - 1, ]
    cropland_expan_weight[is.na(cropland_expan_weight), , ] <- 0
    getYears(cropland_expan_weight) <- getYears(cropland_hr)[t]
    # reduction weight: divide cropland in each grid cell by the
    # total amount of cropland in each cluster (low resolution)
    cropland_reduc_weight <- cropland_hr[, t - 1, ] / land_lr_dagg[, t - 1, "crop"]
    cropland_reduc_weight[is.na(cropland_reduc_weight) | is.infinite(cropland_reduc_weight), , ] <- 0

    # disaggregate cropland for each time step
    cropland_hr[, t, "crop"] <- cropland + cropland_expan * cropland_expan_weight - cropland_reduc * cropland_reduc_weight
  }

  # ========================================================================
  # allocate non-cropland pools
  # ========================================================================

  # calculate total land area in each grid cell (high resolution)
  land_tot_hr <- setYears(dimSums(x_ini_hr, dim = 3), NULL)
  # calculate non-cropland vegetation pool after disaggregation
  land_tot_nocrop_veg_hr <- setNames(setYears(land_tot_hr - x_ini_hr[, , "urban"] - cropland_hr, getYears(cropland_hr)), NULL)
  land_tot_nocrop_veg_hr[land_tot_nocrop_veg_hr < 0] <- 0

  # define sets used during land allocation
  nocrop_pools <- getNames(x_ini_lr)[-which(getNames(x_ini_lr) == "crop")]
  secd_veg <- nocrop_pools[-which(nocrop_pools == "urban" | nocrop_pools == "primforest")]

  # create new magpie object
  land_nocrop_hr <- new.magpie(getCells(cropland_hr), getYears(cropland_hr), nocrop_pools)
  land_nocrop_hr[, year_ini, nocrop_pools] <- x_ini_hr[, , nocrop_pools]

  #------------------------------------------------------------------------
  # allocate urban (constant)
  #------------------------------------------------------------------------

  land_nocrop_hr[, , "urban"] <- x_ini_hr[, , "urban"]

  #------------------------------------------------------------------------
  # allocate primary forest
  #------------------------------------------------------------------------

  for (t in 2:nyears(land_nocrop_hr)) {

    # check where there was cropland expansion
    where_crop_expan <- where(cropland_hr[, t, ] > cropland_hr[, t - 1, ])[["true"]][["regions"]]
    # calculate share of primary forest in previous time step
    shr_prev_primforest_hr <- land_nocrop_hr[, t - 1, "primforest"] / (land_tot_nocrop_veg_hr[, t - 1, ])
    shr_prev_primforest_hr[is.na(shr_prev_primforest_hr) | is.infinite(shr_prev_primforest_hr)] <- 0
    # multiply share of primary forest in previous time step with available land in current time step
    # to allocate primary forest in cells with cropland expansion
    land_nocrop_hr[where_crop_expan, t, "primforest"] <- shr_prev_primforest_hr[where_crop_expan, , ] * land_tot_nocrop_veg_hr[where_crop_expan, t, ]

    # check in which cells there was no cropland expansion
    where_crop_noexpan <- where(cropland_hr[, t, ] > cropland_hr[, t - 1, ])[["false"]][["regions"]]
    # allocate the area of primary forest in the previous time step to those cells
    land_nocrop_hr[where_crop_noexpan, t, "primforest"] <- land_nocrop_hr[where_crop_noexpan, t - 1, "primforest"]

    # sum temporary primary forest pool at low resolution to compare it with the pool in x
    tmp_land_primforest_lr <- speed_aggregate(land_nocrop_hr[, t, "primforest"] * 1e+10, rel)
    tmp_land_primforest_lr_dagg <- speed_aggregate(tmp_land_primforest_lr, t(rel))
    tmp_land_primforest_lr_dagg <- tmp_land_primforest_lr_dagg / 1e+10

    # calculate the residual difference that still needs to be allocated
    primf_residual_diff_lr <- land_lr_dagg[, t, "primforest"] - tmp_land_primforest_lr_dagg
    primf_residual_reduc_lr <- primf_residual_diff_lr
    primf_residual_reduc_lr[primf_residual_reduc_lr > 0] <- 0
    primf_residual_reduc_lr <- abs(primf_residual_reduc_lr)

    # calculate the allocation weight for the residual land reduction:
    # divide temporay primary forest pool at high resolution by sum of the temporary pool at low resolution
    primf_reduc_weight <- land_nocrop_hr[, t, "primforest"] / (tmp_land_primforest_lr_dagg)
    primf_reduc_weight[is.na(primf_reduc_weight) | is.infinite(primf_reduc_weight)] <- 0
    primf_reduc_weight[primf_reduc_weight < 0] <- 0

    # allocate the residual primary forest reduction
    land_nocrop_hr[, t, "primforest"] <- land_nocrop_hr[, t, "primforest"] - primf_residual_reduc_lr * primf_reduc_weight

    # calculate the residual land expansion
    primf_residual_alloc_lr <- primf_residual_diff_lr
    primf_residual_alloc_lr[primf_residual_alloc_lr < 0] <- 0

    # calculate the remaining available for primary land at high resolution
    # primary forest cannot be bigger in the current time step than in the previous time step
    primf_avail_hr <- land_nocrop_hr[, t - 1, "primforest"] - land_nocrop_hr[, t, "primforest"]

    # calculate the allocatiob weight for the residual primary forest land
    # divide the available land by the total area of the residual allocation at low resolution
    primf_alloc_weight <- primf_avail_hr / (primf_residual_alloc_lr)
    primf_alloc_weight[is.na(primf_alloc_weight) | is.infinite(primf_alloc_weight)] <- 0

    # allocate residual primary forest
    land_nocrop_hr[, t, "primforest"] <- land_nocrop_hr[, t, "primforest"] + primf_residual_alloc_lr * primf_alloc_weight
  }

  # calculate remaining secondary vegetation pool after allocation of cropland and primary forest
  land_tot_secd_veg_hr <- land_tot_nocrop_veg_hr - land_nocrop_hr[, , "primforest"]
  getNames(land_tot_secd_veg_hr) <- NULL

  #------------------------------------------------------------------------
  # allocate secondary vegetation land pools
  #------------------------------------------------------------------------

  for (t in 2:nyears(land_nocrop_hr)) {

    # calculate share of non-cropland vegetation pools in previous time step
    shr_prev_nocrop_veg_hr <- land_nocrop_hr[, t - 1, secd_veg] / (land_tot_secd_veg_hr[, t - 1, ])
    shr_prev_nocrop_veg_hr[is.na(shr_prev_nocrop_veg_hr) | is.infinite(shr_prev_nocrop_veg_hr)] <- 0
    # multiply shares of non-cropland pools in previous time step with available land in current time step
    land_nocrop_hr[, t, secd_veg] <- shr_prev_nocrop_veg_hr * land_tot_secd_veg_hr[, t, ]

    # sum temporary non-cropland land pool at low resolution to compare them with the land pools in x
    tmp_land_nocrop_lr <- speed_aggregate(land_nocrop_hr[, t, secd_veg] * 1e+10, rel)
    tmp_land_nocrop_lr_dagg <- speed_aggregate(tmp_land_nocrop_lr, t(rel))
    tmp_land_nocrop_lr_dagg <- tmp_land_nocrop_lr_dagg / 1e+10

    # calculate the residual difference that still needs to be allocated
    residual_diff_lr <- land_lr_dagg[, t, secd_veg] - tmp_land_nocrop_lr_dagg
    # calculate the residual land reduction (absolute)
    residual_reduc_lr <- residual_diff_lr
    residual_reduc_lr[residual_reduc_lr > 0] <- 0
    residual_reduc_lr <- abs(residual_reduc_lr)

    # calculate the allocation weight for the residual land reduction:
    # divide temporay non-croplad pools at high resolution by sum of the temporary at low resolution
    nocrop_reduc_weight <- land_nocrop_hr[, t, secd_veg] / (tmp_land_nocrop_lr_dagg)
    nocrop_reduc_weight[is.na(nocrop_reduc_weight) | is.infinite(nocrop_reduc_weight)] <- 0

    # allocate the residual non-cropland pool reduction
    land_nocrop_hr[, t, secd_veg] <- land_nocrop_hr[, t, secd_veg] - residual_reduc_lr * nocrop_reduc_weight

    # calculate the residual land expansion
    residual_expan_lr <- residual_diff_lr
    residual_expan_lr[residual_expan_lr < 0] <- 0

    # calculate the remaining available land at high resolution
    nocrop_avail_hr <- (land_tot_secd_veg_hr[, t, ] - dimSums(land_nocrop_hr[, t, secd_veg], dim = 3)) * (residual_expan_lr / dimSums(residual_expan_lr, dim = 3))

    # calculate the expansion weight for the residual expansion
    # divide the available land by the total area of the residual expansion at low resolution
    nocrop_expan_weight <- nocrop_avail_hr / (residual_expan_lr)
    nocrop_expan_weight[is.na(nocrop_expan_weight) | is.infinite(nocrop_expan_weight)] <- 0

    # allocate the residual non-cropland pool expansion
    land_nocrop_hr[, t, secd_veg] <- land_nocrop_hr[, t, secd_veg] + residual_expan_lr * nocrop_expan_weight
  }

  # ========================================================================
  # return output
  # ========================================================================

  # combine output
  hr <- mbind(cropland_hr, land_nocrop_hr)
  # correct for small negative deviations
  hr[hr < 0] <- 0

  if (unit == "share") {
    # divide land pools by total amount of land per grid cell
    out <- hr / dimSums(hr, dim = 3)
  } else if (unit == "Mha") {
    out <- hr
  }

  return(out)
}
