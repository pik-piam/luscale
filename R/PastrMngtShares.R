#' Disaggregate managed pastures areas
#'
#' This function calculates pasture management patters based on linear programing accounting for
#' costs of management, area available and production demands.
#'
#' @usage PastrMngtShares(cost, yields, avl_area, production)
#' @param cost Magpie object with cellular cost for each pasture management type
#' @param yields Magpie object with cellular yields for each management type
#' @param avl_area Magpie object with cellular data on the total managed pasture areas available
#' @param production Magpie object with cellular data on the total pasture demand
#' @return Mapie object with the total pasture areas optimized to fulfill pasture demand in the available areas. An extra dimension is added to
#' capture cells that are unfeasible.
#' @author Marcos Alves
#' @examples
#' \dontrun{
#' PastrMngtShares(cost, yields, avl_area, production)
#' }
#'
#' @importFrom lpSolve lp
#' @export
#'
PastrMngtShares <- function(cost, yields, avl_area, production) {
  y <-
    array(
      data = 0,
      dim = c(dim(yields)[c(1, 2)], dim(yields)[3] + 1),
      dimnames = c(dimnames(yields)[c(1, 2)], list("out" = append(
        dimnames(yields)[[3]], "infes"
      )))
    )
  if (!length(getItems(cost, dim = 3)) == length(getItems(yields, dim = 3))) {
    stop("costs and yeilds have different sizes")
  }
  if (length(getItems(avl_area, dim = 3)) > 1) {
    stop("Area has more than one data dimention")
  }
  if (length(getItems(production, dim = 3)) > 1) {
    stop("production has more than one data dimention")
  }

  x <- mbind(cost, yields, avl_area, production)
  start <- Sys.time()
  for (year in getItems(x, dim = 2)) {
    print(paste0("Optimizing year: ", year))
    core_opt <- function(x) {
      infes <- 0
      i <- length(x)
      n <- (i - 2) / 2
      if ((x[i - 1] != 0)) {
        f.obj <- x[seq_len(n)]
        f.con <- rbind(
          rep(1, n),
          x[seq(n + 1, n * 2, 1)]
        )
        f.dir <- c(
          "==",
          "=="
        )
        f.rhs <- c(
          x[i - 1],
          x[i]
        )
        res <- lp("min", f.obj, f.con, f.dir, f.rhs)
        if (res$status == 2) {
          j <- match(max(x[seq(n + 1, n * 2, 1)]), x[seq(n + 1, n * 2, 1)])
          res$solution[j] <- x[i - 1]
          infes <- 1
        }
        return(c(res$solution, infes))
      }
      return(rep(0, n + 1))
    }
    z <- unlist(apply(x[, year, ], 1, core_opt))
    z <- t(z)
    y[, year, ] <- z
  }
  y <- as.magpie(y)
  end <- Sys.time()
  print(end - start)
  return(y)
}