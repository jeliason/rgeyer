#' Compute the components of papangelou for the stepper model
#'
#' @param from coordinates of the data
#' @param to coordinate matrix for locations for evaluation of papangelou. Assumed to be different from 'from'
#' @param K the number of steps
#' @param R maximum range
#' @param r alternative to K and R, the r vector itself. Increasing, 0 shouldn't be included.
#' @param bbox bounding box
#' @param dbg dbg level
#' @param toroidal >0 use toroidal distances
#'
#' @return Matrix with dimension (nrow(to), K) with
#' d_k(to_i|from) = 1( ne_k(to_i, from)>0) with ne_k(x,from) = sum_j 1(r[k-1] < ||to_i - from_j|| < r[k])
#'
#' If 'to' is NULL, computes the components for 'from' locations, i.e. the papangelou for each data point.
#' @export
#' @useDynLib rgeyer

stepper_components <- function(from, to, R, K, r, bbox, dbg=0, toroidal=0){
  if(missing(r)) r <- (1:K) * R /K
  if(is.null(to) | missing(to)){
    rstepper_components_at_data_c(
      from,
      r,
      bbox,
      dbg,
      toroidal)
  }
  else{
  rstepper_components_c(
    from,
    to,
    r,
    bbox,
    dbg,
    toroidal)
  }
}
