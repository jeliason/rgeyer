#' Compute the components of papangelou for the stepper model, from a grid to data
#'
#' Will compute the components from all grid points to data, for each type of grid point.
#'
#' @param from coordinates of the data
#' @param to coordinate matrix for locations for evaluation of papangelou. Assumed to be different from 'from'
#' @param K the number of steps
#' @param R maximum range
#' @param r alternative to K and R, the r vector itself. Increasing, 0 shouldn't be included.
#' @param sat the saturation levels, >0 integers.
#' @param bbox bounding box
#' @param dbg dbg level
#' @param toroidal >0 use toroidal distances
#' @param ... ignored
#'
#' @return Matrix with dimension (nrow(to), K) with
#' d_k(to_i|from) = 1( ne_k(to_i, from)>0) with ne_k(x,from) = sum_j 1(r[k-1] < ||to_i - from_j|| < r[k])
#'
#' If 'to' is NULL, computes the components for 'from' locations, i.e. the papangelou for each data point.
#' @export
#' @useDynLib rgeyer

stepper_multi_components_grid <- function(from, to, ranges1, ranges2, sat1, sat2, bbox, dbg=0, toroidal=0, ...){
    rstepper_components_multi_grid_c(
      from,
      to,
      ranges1,
      ranges2,
      sat1,
      sat2,
      bbox,
      dbg,
      toroidal)
}
