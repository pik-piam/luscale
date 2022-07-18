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
  area  <-  land_hr > 0
  remain_prod <- grass_prod_hr - excess_prod

  print(paste("Excess poll"," -> ","Regular pool"," -> ", "Erased biomass"))
  print(paste(sum(excess_prod)," -> ",sum(remain_prod), " -> ", 0))
  while (sum(excess_prod) > 1 & count < ite) {
    reserve  <-  poten_prod > remain_prod * (1 + 0.01) # solution cells (partial)
    avl_reseve_areas  <-  area * reserve * (excess_prod == 0) # all areas where there is production but no excess and productivity and reserve capacity
    avl_area_weight <- weight * (avl_reseve_areas > 0)



    excess_prod_cluster <- toolAggregate(excess_prod, rel = map_file, from = "cell", to = "country")
      red_prod_cell <- toolAggregate(excess_prod_cluster,
                                     rel = map_file,
                                     weight = avl_area_weight,
                                     from = "country", to = "cell"
      )

    #reserve_cty <- (avl_reseve_areas > 0) * (poten_prod - remain_prod)
    #reserve_cty <- sapply(sort(getRegions(reserve_cty)), function(x) dimSums(reserve_cty[x,,], dim=c(1)))
    #excess_prod_cty <- sapply(sort(getRegions(excess_prod_cluster)), function(x) dimSums(excess_prod_cluster[x,,], dim=c(1)))
    #print(round(data.frame(reserve_cty-excess_prod_cty)))
    #print(paste0("total reserve: ", sum(reserve_cty)))
    #print(paste0("total excess: ", sum(excess_prod_cty)))

    grass_prod_hr_tmp <- red_prod_cell + remain_prod
    excess_prod <- grass_prod_hr_tmp - poten_prod
    excess_prod[excess_prod < 0] <- 0
    remain_prod <- grass_prod_hr_tmp - excess_prod

    print(paste(sum(excess_prod)," -> ",sum(remain_prod), " -> ", sum(grass_prod_hr - grass_prod_hr_tmp)))
    count <- count + 1
  }
  return(remain_prod)
}
