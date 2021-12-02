#' Disaggregate managed pastures areas
#'
#' This function calculates pasture management patters based on linear programing accounting for
#' costs of management, area available and production demands.
#'
#' @usage GrassProdDisagg(grass_prod_lr,land_hr, weight_factor, map_file, lpjml_yields, ite)
#' @param grass_prod_lr grassland production in low resolution to be disaggregated
#' @param land_hr Magpie object with celular cost for each pasture management type
#' @param weight_factor Magpie object with celular yields for each management type
#' @param map_file cluster mappint to cells and regions
#' @param lpjml_yields Lpjml yeilds used for magpie optimization
#' @param ite number of redistribution iterations
#' @return Mapie object with the total pasture areas optimized for fulfulling pasture demand in the available areas. An extra dimention is added to
#' capture cells that are infeasible.
#' @author Marcos Alves
#' @examples
#' \dontrun{
#' GrassProdDisagg(grass_prod_lr,land_hr, weight_factor, map_file, lpjml_yields, ite)
#' }
#'
#' @importFrom lpSolve lp
#' @export
#'

GrassProdDisagg <- function(grass_prod_lr,land_hr, weight_factor, map_file, lpjml_yields, ite = 10) {
  poten_prod <- lpjml_yields * land_hr
  
  name <- getItems(land_hr, dim = 3)
  if(length(name) > 1){
    stop("Cannot handle more than one type of grassland at the same time (land_hr)..")
  }
  
  if(dim(weight_factor)[3] > 1){
    stop("Cannot handle more than one type of weight factors at the same time.")
  }
  
  if(dim(grass_prod_lr)[3] > 1){
    stop("Cannot handle more than one type of grassland at the same time (grass_prod_lr).")
  }
  
  if(dim(lpjml_yields)[3] > 1){
    stop("Cannot handle more than one type of grassland at the same time (lpjml_yields).")
  }

  if(!getItems(grass_prod_lr, dim = 3) == getItems(land_hr, dim = 3)){
    stop("Datasets with different names")
    if(!getItems(land_hr, dim = 3) == getItems(lpjml_yields, dim = 3)){
      stop("Datasets with different names")
    }
  }
  
  if (name == "pastr") {
    weight <- land_hr * weight_factor # pasture suitability
    print("Pasture")
  } else {
    weight <- land_hr / weight_factor # accessibility
    weight[is.nan(weight) | is.infinite(weight)] <- 0
    print("rangelands")
  }

  for (i in getYears(weight)) {
    weight[, i, ] <- (weight[, i, ] - min(weight[, i, ])) / (max(weight[, i, ]) - min(weight[, i, ]))
  }
  weight[is.nan(weight) | is.infinite(weight)] <- 0

  grass_prod_hr <- toolAggregate(grass_prod_lr, rel = map_file, weight = weight, from = "cluster", to = "cell")

  # biomass adjustment
  excess_prod <- grass_prod_hr - poten_prod
  excess_prod[excess_prod < 0] <- 0
  count <- 0

  while (sum(excess_prod) > 5 & count < ite) {
    remain_prod <- grass_prod_hr - excess_prod
    excess_prod_cluster <- toolAggregate(excess_prod, rel = map_file, from = "cell", to = "region")
    red_prod_cell <- toolAggregate(excess_prod_cluster,
      rel = map_file,
      weight = weight,
      from = "region", to = "cell"
    )
    grass_prod_hr <- red_prod_cell + remain_prod
    excess_prod <- grass_prod_hr - poten_prod
    excess_prod[excess_prod < 0] <- 0
    print(sum(excess_prod))
    count <- count + 1
  }
  return(grass_prod_hr)
}
