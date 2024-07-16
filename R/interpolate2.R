#' Interpolate 2
#'
#' Disaggregates a cellular MAgPIE output to 0.5 degree based on the given mapping
#' and information about the initial 0.5 degree distribution.
#'
#' There was a now deleted function called interpolate in this package before,
#' hence the name interpolate2.
#'
#' The function is based on the following assumption: \code{x} is an object in
#' low resolution with more than one data dimension and the sum over the data
#' dimensions is constant over time. Example: One column cropland, the other
#' one (cell size - cropland). \code{x_ini} provides the same type of data
#' as \code{x} but in high resolution and for the time step previous the initial
#' time step of x (e.g. if \code{x} goes from t=a to t=a+10, \code{x_ini} must
#' be provided for t=a-1). The function calculates the amount by which the
#' individual data columns of x change in each timestep. The output is based
#' on \code{x_ini} and only the differences in later timesteps to ths starting
#' point are disaggregated by the given mapping. This assures that as little
#' information as possible is lost from the original dataset \code{x_ini}.
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
#' @param x The object to be disaggregated. See details for further important
#' information.
#' @param x_ini The initial distribution of x in high resolution.
#' @param map A relation map between low and high resolution
#' @param x_ini_lr Low resolution version of x_ini. Will be calculated automatically, if not provided.
#' Can speed up computation, if provided.
#' @return The disaggregated MAgPIE object containing x_ini as first
#' timestep
#' @author Jan Philipp Dietrich
#' @importFrom magclass is.magpie nregions nyears getNames getYears<- getCells<- mbind
#' @importFrom magclass dimSums setYears getYears as.magpie new.magpie as.array
#' @importFrom madrat toolAggregate
#' @seealso \code{\link[madrat]{toolAggregate}}
#' @export
interpolate2 <- function(x, x_ini, map, x_ini_lr = NULL) { # nolint: object_name_linter.
  xIni <- x_ini
  xIniLr <- x_ini_lr
  if (!is.magpie(x) || !is.magpie(xIni)) stop("x and x_ini have to be magpie objects")
  if (is.null(getYears(x))) stop("Years must be defined in x.")
  if (nyears(xIni) > 1) stop("Initialization data must only contain one timestep")
  if (!setequal(getNames(x), getNames(xIni))) stop("data names of x and x_ini must be identical")

  x <- x[, sort(getYears(x, as.integer = TRUE)), ]
  if (is.null(getYears(xIni))) {
    getYears(xIni) <- getYears(x, as.integer = TRUE)[1] - 1
  } else if (getYears(xIni, as.integer = TRUE) >=  getYears(x, as.integer = TRUE)[1]) {
    stop("Year of x_ini must be smaller than start year of x")
  }

  if (is.null(xIniLr)) xIniLr <- madrat::toolAggregate(xIni, map, from = "cell", to = "cluster")

  lr <- mbind(xIniLr, x)
  # Test if the total sum is constant
  test <- dimSums(lr, dim = 3)
  if (any(abs(test - test[, 1, ]) > 0.0001)) warning("Total stock is not constant over time. See help for details")

  # calculate reduction and extension shares which then can be disaggregated
  more <- less <- lr[, 2:nyears(lr), ] - setYears(lr[, 1:(nyears(lr) - 1), ], getYears(lr)[2:nyears(lr)])
  less[less > 0] <- 0
  more[more < 0] <- 0
  reduct <- -less / setYears((lr[, 1:(nyears(lr) - 1), ] + 10^-100), getYears(less))
  avail  <- dimSums(-less, dim = 3)
  extent <- more / (avail + 10^-100)

  # disaggregate shares (combined disaggregation is faster than separate disaggregation
  # for both data sets)
  extentHr <- reductHr <- madrat::toolAggregate(extent - reduct, map, from = "cluster", to = "cell")
  extentHr[extentHr < 0] <- 0
  reductHr[reductHr > 0] <- 0
  reductHr <- -reductHr

  # calculate land pools in high res (hr)
  hr <- new.magpie(getCells(xIni), c(getYears(xIni), getYears(reductHr)), getNames(xIni))
  hr[, 1, ] <- xIni

  for (y in 1:nyears(reductHr)) {
    hr[, y + 1, ] <- (1 - reductHr[, y, ]) * setYears(hr[, y, ], NULL) +
      dimSums(reductHr[, y, ] * setYears(hr[, y, ], NULL), dim = 3) * extentHr[, y, ]
  }
  return(hr)
}
