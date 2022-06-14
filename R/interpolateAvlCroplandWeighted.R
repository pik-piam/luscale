#' interpolateAvlCroplandWeighted
#'
#' Disaggregates a modelled time series of land pools after optimisation from the model resolution (low resolution)
#' to the resolution of the land initialisation data set (high resolution), based on a relation map and available
#' cropland.
#'
#' The function requires the following input data:
#' \itemize{
#' \item \code{x} is an object containing a time series of land pools (model output). The sum over all land pools
#' is constant over time.
#' \item \code{x_ini_lr} and \code{x_ini_hr} provide the initial land pools (Mha) at high (hr) and low resolution (lr)
#' before the optimisation. They only contain the initial time step, but share the three-dimensional structure
#' with \code{x}.
#' \item  \code{avl_cropland_hr} provides information about the amount (Mha) of available cropland at high resolution.
#' \item \code{map} relation map containing information about cell belongings between the high and low resolution.
#' }
#'
#' The weighted disaggregation works as follows:
#'
#' 1. The share of cropland in terms of total available cropland is calculated at the previous time step and
#' then multiplied by the available cropland at the current time step (as available cropland can change over time
#' - e.g. by policy restriction as can be specified in \code{snv_pol_shr}). This temporary cropland pool is then
#' compared to the low resolution cropland pool and the residual area of cropland expansion and reduction is
#' determined.
#'
#' 2. In order to allocate residual area of cropland expansion and reduction, for each grid cell at high resolution
#' expansion and reduction weights are calculated and multiplied by the residual area:
#' \itemize{
#' \item The reduction weight is given by the ratio between the amount of cropland per grid cell and the total area
#' of the temporary cropland at the low resolution spatial unit. This assumes that the cropland reduction is equally
#' distributed among all high resolution grid cells.
#' \item The expansion weight is calculated as the ratio between the remaining cropland at the grid cell level
#' (high resolution) and the overall remaining cropland at the low resolution spatial unit in the current time step.
#' The remaining cropland given by the difference between the available cropland and the temporaryl cropland pool
#' minus urban land, since it assumed that cropland cannot be allocated to urban land.
#' }
#'
#' 3. Following the cropland allocation, the land area for the remaining non-cropland vegetation pools is calculated
#' by substracting the allocated cropland and urban land areas from the total land area in each grid cell.
#'
#' 4. The non-cropland vegetation pool at the high resolution (except of primary forest), calculated in step 3., is
#' then multiplied by the respective shares of the remaining non-cropland vegetation pools at the previous time step
#' (temporary allocation). Similar to the cropland allocation, is not sufficient to also account for changes within
#' these land pools. Therefore, the temporarily allocated non-cropland pools are, once again, compared with the pools
#' at low resolution. The residual area of land expansion and reduction is then allocated by based on reduction and
#' expansion weights, similar as in 2.. The reduction weight is calculated as the ratio between the given temporary
#' land pool at high resolution and total temporary land pool at low resolution. The expansion weight is calculated
#' as the ratio between the remaining land to be filled in each land pool and the total amount of residual land to
#' be allocated in the current time step.
#'
#' 5. Primary forest is treated in a slightly different way, as primary forest cannot be expanded over time. In
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
#' @param avl_cropland_hr The area of available cropland at the high resolution.
#' @param urban_land_hr Either a magpie object of the cellular urban input data, or "static" string
#' @param map A relation map between low and high resolution
#' @param marginal_land Depending on the cropland suitability data, standard options are
#' \itemize{
#' \item \code{"all_marginal"}: Cropland can be allocated to marginal land
#' \item \code{"q33_marginal"}: The bottom tertile of the marginal land area is excluded
#' \item \code{"no_marginal"}: Marginal land is fully excluded from cropland
#' }
#' @param snv_pol_shr Share of available cropland that is witheld for other land cover types. Can be supplied as a
#' single value or as a magpie object containing different values in each iso country.
#' @param snv_pol_fader Fader for share of set aside policy.
#' @param year_ini Timestep that is assumed for the initial distributions \code{x_ini_hr} and \code{x_ini_lr}.
#' @param unit Unit of the output. "Mha" or "share"
#' @param set_aside_shr depreciated, stop using.
#' @param set_aside_fader depreciated, stop using.
#' @return The disaggregated MAgPIE object containing x_ini_hr as first
#' timestep
#' @export
#' @author Patrick von Jeetze, David Chen
#' @importFrom magclass is.magpie nregions nyears getItems add_columns getNames getYears mbind dimSums setYears getYears new.magpie where
#' @importFrom madrat toolAggregate toolGetMapping toolConditionalReplace
#' @seealso \code{\link{interpolate2}}
#' \code{\link{toolAggregate}}
#' @examples
#' \dontrun{
#' a <- interpolateAvlCroplandWeighted(
#'   x = land,
#'   x_ini_lr = land_ini_lr,
#'   x_ini_hr = land_ini_hr,
#'   avl_cropland_hr = "avl_cropland_0.5.mz",
#'   map = "clustermap_rev4.59_c200_h12.rds",
#'   marginal_land = "all_marginal"
#' )
#'
#' sf <- read.magpie("f30_scenario_fader.csv")[, , "by2030"]
#'
#' b <- interpolateAvlCroplandWeighted(
#'   x = land,
#'   x_ini_lr = land_ini_lr,
#'   x_ini_hr = land_ini_hr,
#'   avl_cropland_hr = "avl_cropland_0.5.mz",
#'   map = "clustermap_rev4.59_c200_h12.rds",
#'   marginal_land = "all_marginal",
#'   snv_pol_shr = 0.2,
#'   snv_pol_fader = sf
#' )
#'
#' iso <- readGDX(gdx, "iso")
#' set_aside_iso <- readGDX(gdx, "policy_countries30")
#' set_aside_select <- readGDX(gdx, "s30_snv_shr")
#' set_aside_noselect <- readGDX(gdx, "s30_snv_shr_noselect")
#' snv_pol_shr <- new.magpie(iso, fill = snv_noselect)
#' snv_pol_shr[set_aside_iso, , ] <- set_aside_select
#'
#' c <- interpolateAvlCroplandWeighted(
#'   x = land,
#'   x_ini_lr = land_ini_lr,
#'   x_ini_hr = land_ini_hr,
#'   avl_cropland_hr = "avl_cropland_0.5.mz",
#'   map = "clustermap_rev4.59_c200_h12.rds",
#'   marginal_land = "all_marginal",
#'   snv_pol_shr = snv_pol_shr,
#'   snv_pol_fader = saf
#' )
#' }
#'
interpolateAvlCroplandWeighted <- function(x, x_ini_lr, x_ini_hr, avl_cropland_hr, map, urban_land_hr = "static",
                                           marginal_land = "all_marginal", snv_pol_shr = 0,
                                           snv_pol_fader = NULL, year_ini = "y1985", unit = "Mha",set_aside_shr=NULL,set_aside_fader=NULL) {
  if(!is.null(set_aside_shr)|!is.null(set_aside_fader)){
    warning("arguments set_aside_shr and set_aside_fader are depreciated, stop using them")
    snv_pol_fader <- set_aside_fader
    snv_pol_shr <- set_aside_shr
  }

  # test whether data can be handled by function
  if (!is.magpie(x) || !is.magpie(x_ini_lr) || !is.magpie(x_ini_hr) || (!is.magpie(urban_land_hr) & urban_land_hr != "static")) stop("x, x_ini_lr and x_ini_hr, urban_land_hr have to be magpie objects")
  if (nregions(x) != nregions(x_ini_lr)) stop("x and x_ini_lr have to be of the same spatial aggregation")
  if (nyears(x_ini_lr) > 1 || nyears(x_ini_hr) > 1) stop("Initialization data must only have one timestep")
  if (!all(getNames(x) == getNames(x_ini_lr)) || !all(getNames(x) == getNames(x_ini_hr))) stop("dimnames[[3]] of x, x_ini_lr and x_ini_hr have to be the same")
  if (!file.exists(map)) stop("relation map file ", map, " not found")

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
  # read available cropland data
  #------------------------------------------------------------------------

  if (is.character(avl_cropland_hr)) {
    if (!file.exists(avl_cropland_hr)) stop("high resolution available cropland data not found")
    # high resolution
    avl_cropland_hr <- read.magpie(avl_cropland_hr)
  }

  # expand available cropland data over time
  # high resolution
  avl_cropland_hr_tmp <- new.magpie(getCells(avl_cropland_hr), getYears(lr), marginal_land)
  avl_cropland_hr_tmp[, getYears(lr), ] <- avl_cropland_hr[, , marginal_land]
  avl_cropland_hr <- avl_cropland_hr_tmp

  if (any(snv_pol_shr != 0) & is.null(snv_pol_fader)) stop("Share of withheld cropland given, but no policy fader for target year provided")

  if (length(map) == 1) {
    map <- toolGetMapping(map, where = "local")
  }
  if (length(snv_pol_shr == 1)) {
    snv_pol_shr <- new.magpie(map[, "cell"], fill = snv_pol_shr)
  } else {
    snv_pol_shr <- toolAggregate(snv_pol_shr[unique(map[, "country"]), , ], map, from = "country", to = "cell")
  }

  if (!is.null(snv_pol_fader) & is.magpie(snv_pol_fader)) {
    if (ndata(snv_pol_fader) != 1) stop("snv_pol_fader has too many data dimensions. Please select one target year only for this disaggregation.")
    # correct available cropland with policy restriction
    for (t in 1:nyears(lr)) {
      avl_cropland_hr[, t, ] <- avl_cropland_hr[, t, ] * (1 - snv_pol_shr * snv_pol_fader[, getYears(lr)[t], ])
    }
  } else if (!is.null(snv_pol_fader) & !is.magpie(snv_pol_fader)) {
    if (ncol(snv_pol_fader) != 1) stop("snv_pol_fader has too many columns. Please select one target year only for this disaggregation.")
    # correct available cropland with policy restriction
    for (t in 1:nyears(lr)) {
      avl_cropland_hr[, t, ] <- avl_cropland_hr[, t, ] * (1 - snv_pol_shr * snv_pol_fader[getYears(lr)[t], ])
    }
  }

  # correct for urban land to constrain the calculation of the allocation weight
  # This is due to inconsistency in input dataset where in some regions urban land is already subtracted from avl Cropland
  # In other regions it is not subtracted even when urban centres there exist
  # high resolution
  land_non_urban_hr <- dimSums(x_ini_hr, dim = 3) - x_ini_hr[, , "urban"]
  getCells(avl_cropland_hr) <- getCells(avl_cropland_hr_tmp) <- getCells(land_non_urban_hr)
  # where available cropland is larger than total non urban land chose the smaller value [pmin()]
  for (t in 1:nyears(lr)) {
    avl_cropland_hr_tmp[, t, ] <- pmin(avl_cropland_hr[, t, ], land_non_urban_hr)
  }
  avl_cropland_hr <- avl_cropland_hr_tmp

  #
  ##### remove urban cells from potentially available cropland, if not static ###
  #
   if (is.magpie(urban_land_hr)){

    #check if urban land output and input correspond at low res level
    urban_input_lr <-  toolAggregate(urban_land_hr, map, from = "cell", to = "cluster")
    testurb    <- x[,,"urban"]  - urban_input_lr[getItems(x,dim=1),getYears(x),]

     if (max(testurb) > 0.05*diff(range(urban_input_lr)) || min(testurb) < (-0.05*diff(range(urban_input_lr)))){
        warning("Cluster level differences of urban land between magpie and input data are
          greater than 5% of the total range of urban land values in a cluster.
          Consider the urban implementation, maybe something's wrong.")}

    # remove urban land expansion from future available cropland
    urb_ex <- urban_land_hr - setYears(x_ini_hr[,,"urban"], NULL)
    urb_ex[urb_ex < 0 ] <- 0 #negatives from rounding

    avl_cropland_hr[,2:nyears(avl_cropland_hr),] <-  avl_cropland_hr[,2:nyears(avl_cropland_hr),] - setNames(urb_ex[, getYears(avl_cropland_hr)[-1],],NULL)
    avl_cropland_hr[avl_cropland_hr < 0] <- 0
  }



  #------------------------------------------------------------------------
  # disaggregate low resolution output data to be used in calculations
  #------------------------------------------------------------------------

  # total land expansion and reduction at low resolution
  land_reduc_lr_dagg <- toolAggregate(land_reduc_lr, map, from = "cluster", to = "cell")
  land_expan_lr_dagg <- toolAggregate(land_expan_lr, map, from = "cluster", to = "cell")
  # total amount of cropland at low resolution
  land_lr_dagg <- toolAggregate(lr, map, from = "cluster", to = "cell")

  # ========================================================================
  # allocate cropland at high res (hr)
  # ========================================================================

  # create new magpie object for cropland allocation (high resolution)
  cropland_hr <- new.magpie(getCells(avl_cropland_hr), getYears(lr), "crop")
  getCells(x_ini_hr) <- getCells(avl_cropland_hr)
  cropland_hr[, year_ini, ] <- x_ini_hr[, , "crop"]

  for (t in 2:nyears(cropland_hr)) {

    # calculate share of non-cropland vegetation pools in previous time step
    shr_prev_cropland_hr <- cropland_hr[, t - 1, "crop"] / (avl_cropland_hr[, t - 1, ])
    shr_prev_cropland_hr[is.na(shr_prev_cropland_hr) | is.infinite(shr_prev_cropland_hr)] <- 0
    # multiply shares of non-cropland pools in previous time step with available land in current time step
    cropland_hr[, t, "crop"] <- shr_prev_cropland_hr * avl_cropland_hr[, t, ]

    # sum temporary non-cropland land pool at low resolution to compare them with the land pools in x
    tmp_cropland_lr <- toolAggregate(cropland_hr[, t, "crop"] * 1e+10, map, from = "cell", to = "cluster")
    tmp_cropland_lr_dagg <- toolAggregate(tmp_cropland_lr, map, from = "cluster", to = "cell")
    tmp_cropland_lr_dagg <- tmp_cropland_lr_dagg / 1e+10

    # calculate the residual difference that still needs to be allocated
    residual_diff_lr <- land_lr_dagg[, t, "crop"] - tmp_cropland_lr_dagg
    # calculate the residual land reduction (absolute)
    residual_reduc_lr <- residual_diff_lr
    residual_reduc_lr[residual_reduc_lr > 0] <- 0
    residual_reduc_lr <- abs(residual_reduc_lr)

    # calculate the allocation weight for the residual land reduction:
    # divide temporay non-croplad pools at high resolution by sum of the temporary at low resolution
    crop_reduc_weight <- cropland_hr[, t, "crop"] / (tmp_cropland_lr_dagg)
    crop_reduc_weight[is.na(crop_reduc_weight) | is.infinite(crop_reduc_weight)] <- 0

    # allocate the residual non-cropland pool reduction
    cropland_hr[, t, "crop"] <- cropland_hr[, t, "crop"] - residual_reduc_lr * crop_reduc_weight

    # calculate the residual land expansion
    residual_expan_lr <- residual_diff_lr
    residual_expan_lr[residual_expan_lr < 0] <- 0

    # calculate the remaining available land at high resolution
    cropland_remain_hr <- (avl_cropland_hr[, t, ] - cropland_hr[, t, "crop"])

    # sum remaining cropland in the current time step to calculate expansion weight
    tmp_remain_lr <- toolAggregate(cropland_remain_hr * 1e+10, map, from = "cell", to = "cluster")
    tmp_remain_lr_dagg <- toolAggregate(tmp_remain_lr, map, from = "cluster", to = "cell")
    tmp_remain_lr_dagg <- tmp_remain_lr_dagg / 1e+10

    # calculate the expansion weight for the residual expansion
    # divide the available land by the total area of the residual expansion at low resolution
    crop_expan_weight <- cropland_remain_hr / tmp_remain_lr_dagg
    crop_expan_weight[is.na(crop_expan_weight) | is.infinite(crop_expan_weight)] <- 0

    # allocate the residual non-cropland pool expansion
    cropland_hr[, t, "crop"] <- cropland_hr[, t, "crop"] + residual_expan_lr * crop_expan_weight
  }


  # ========================================================================
  # allocate non-cropland pools
  # ========================================================================

  # calculate total land area in each grid cell (high resolution)
  land_tot_hr <- setYears(dimSums(x_ini_hr, dim = 3), NULL)
  # calculate non-cropland vegetation pool after disaggregation
  #sum urban and cropland
  if (is.magpie(urban_land_hr)){
     if (as.integer(gsub("y", "", year_ini)) < getYears(urban_land_hr, as.integer = TRUE)[1]){
     urban_land_hr <- add_columns(urban_land_hr, addnm=year_ini, dim=2.1, fill=0)
     urban_land_hr <- urban_land_hr[,sort(getYears(urban_land_hr, as.integer = TRUE)),]
     urban_land_hr[,year_ini,] <- urban_land_hr[,2,]

   land_tot_nocrop_veg_hr <- setNames(setYears(land_tot_hr - urban_land_hr[,getYears(cropland_hr),] - cropland_hr, getYears(cropland_hr)), NULL)
   }} else {
   land_tot_nocrop_veg_hr <- setNames(setYears(land_tot_hr - x_ini_hr[, , "urban"] - cropland_hr, getYears(cropland_hr)), NULL)
     }

  land_tot_nocrop_veg_hr[land_tot_nocrop_veg_hr < 0] <- 0

  # define sets used during land allocation
  nocrop_pools <- getNames(x_ini_lr)[-which(getNames(x_ini_lr) == "crop")]
  secd_veg <- nocrop_pools[-which(nocrop_pools == "urban" | nocrop_pools == "primforest")]

  # create new magpie object
  land_nocrop_hr <- new.magpie(getCells(cropland_hr), getYears(cropland_hr), nocrop_pools)
  land_nocrop_hr[, year_ini, nocrop_pools] <- x_ini_hr[, , nocrop_pools]


  #------------------------------------------------------------------------
  # allocate urban (exogenous, dataset read in from input data)
  #------------------------------------------------------------------------
  if(is.magpie(urban_land_hr)){
  land_nocrop_hr[, , "urban"] <- urban_land_hr[, getYears(land_nocrop_hr), "urban"]
  } else {  land_nocrop_hr[, , "urban"] <- x_ini_hr[, , "urban"] }

  #------------------------------------------------------------------------
  # allocate primary forest
  #------------------------------------------------------------------------

  for (t in 2:nyears(land_nocrop_hr)) {

    # check where there was cropland and urban expansion
    where_expan <- where(cropland_hr[, t, ] > cropland_hr[, t - 1, ] | land_nocrop_hr[, t, "urban"] >  land_nocrop_hr[, t - 1, "urban"])[["true"]][["regions"]]
    # calculate share of primary forest in previous time step
    shr_prev_primforest_hr <- land_nocrop_hr[, t - 1, "primforest"] / (land_tot_nocrop_veg_hr[, t - 1, ])
    shr_prev_primforest_hr[is.na(shr_prev_primforest_hr) | is.infinite(shr_prev_primforest_hr)] <- 0
    # multiply share of primary forest in previous time step with available land in current time step
    # to allocate primary forest in cells with cropland expansion
    land_nocrop_hr[where_expan, t, "primforest"] <- shr_prev_primforest_hr[where_expan, , ] * land_tot_nocrop_veg_hr[where_expan, t, ]

    # check in which cells there was no cropland nor urban expansion
    where_noexpan <- where(cropland_hr[, t, ] > cropland_hr[, t - 1, ] | land_nocrop_hr[, t, "urban"] >  land_nocrop_hr[, t - 1, "urban"] )[["false"]][["regions"]]
    # allocate the area of primary forest in the previous time step to those cells
    land_nocrop_hr[where_noexpan, t, "primforest"] <- land_nocrop_hr[where_noexpan, t - 1, "primforest"]

    # sum temporary primary forest pool at low resolution to compare it with the pool in x
    tmp_land_primforest_lr <- toolAggregate(land_nocrop_hr[, t, "primforest"] * 1e+10, map, from = "cell", to = "cluster")
    tmp_land_primforest_lr_dagg <- toolAggregate(tmp_land_primforest_lr, map, from = "cluster", to = "cell")
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

    ## calculate the remaining available for primary land at high resolution
     # primary forest cannot be bigger in the current time step than in the previous time step
     primf_avail_hr <- land_nocrop_hr[, t - 1, "primforest"] - land_nocrop_hr[, t, "primforest"]

     #Check if there is still  land available in this cell for the amount of forest required
     land_left_hr <-  setNames(land_tot_hr - setYears(dimSums(land_nocrop_hr[,t , c("urban","primforest")],dim=3),NULL) -
                                                     setNames(setYears(cropland_hr[,t,],NULL), NULL), NULL)
     # weight by either gap between current tmp and previous time step, or by amount of land remaining
     primf_avail_hr_constr <- pmin(primf_avail_hr, land_left_hr)
     primf_avail_hr_constr[primf_avail_hr_constr < 0 ] <- 0

    # calculate the allocatiob weight for the residual primary forest land
    # divide the available land by the total area of the residual allocation at low resolution
    primf_alloc_weight <- primf_avail_hr_constr / (primf_residual_alloc_lr)
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
    tmp_land_nocrop_lr <- toolAggregate(land_nocrop_hr[, t, secd_veg] * 1e+10, map, from = "cell", to = "cluster")
    tmp_land_nocrop_lr_dagg <- toolAggregate(tmp_land_nocrop_lr, map, from = "cluster", to = "cell")
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
    # cell GRL.13164 has a total land area of 0
    out <- toolConditionalReplace(out, "is.na()", replaceby = 0)
  } else if (unit == "Mha") {
    out <- hr
  }

  return(out)
}
