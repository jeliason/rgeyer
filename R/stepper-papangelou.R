#' Compute the papangelou conditional intensity for the multivariate stepper model
#'
#' @param from coordinates of the data, 3rd column 0-based integer marks
#' @param to coordinate matrix for locations for evaluation of papangelou, 0-based integer marks in 3rd column. Assumed to be different from 'from'. If NULL or missing, compute at "from".
#' @param theta1 list of intra-type parameters, as in simulation (see rstepper_multi)
#' @param theta2 list of inter-type parameters
#' @param bbox bounding box, defines the window.
#' @param dbg dbg level
#' @param toroidal >0 use toroidal distances
#' @param multi_to If TRUE, treat 'to' as unmarked locations at which to evaluate the papangelou for each mark level.
#' @details The first order effects (trends, baselines) are not included, make sure add them afterwards.
#'
#' @return log-conditional intensity l(to_i|from) = theta[0] + sum( theta[k] * d_k(to_i,from) )
#' where d_k(to_i|from) = 1( ne_k(to_i, from)>0) and ne_k(x,from) = sum_j 1(r[k-1] < ||to_i - from_j|| < r[k])
#'
#' If multi_to = TRUE, a matrix with each column corresponding to each mark level.
#'
#' @export
#' @useDynLib rgeyer


stepper_multi_log_papangelou <- function(from, to, theta1, theta2, bbox, dbg=0, toroidal=0, multi_to = FALSE){
  #
  if(missing(to) || is.null(to)){
    rstepper_multi_log_papangelou_at_data_c(
    theta1,
    theta2,
    from,
    bbox,
    dbg,
    toroidal)
  }
  else if(!multi_to){
    rstepper_multi_log_papangelou_c(
    theta1,
    theta2,
    from,
    to,
    bbox,
    dbg,
    toroidal)
  }
  else{
    types <- sort( unique(from[,3]) )
    out <- rstepper_multi_log_papangelou_grid_c(
      theta1,
      theta2,
      from,
      to,
      types,
      bbox,
      dbg,
      toroidal)
    colnames(out) <- types
    out
  }
}

#' Compute the papangelou conditional intensity for the stepper model
#'
#' @param from coordinates of the data
#' @param to coordinate matrix for locations for evaluation of papangelou. Assumed to be different from 'from'
#' @param theta vector of weights, K+1 long
#' @param r steps, increasing, does not include 0
#' @param sat saturation level
#' @param bbox bounding box
#' @param dbg dbg level
#' @param toroidal >0 use toroidal distances
#'
#' @return log-conditional intensity l(to_i|from) = theta[0] + sum( theta[k] * d_k(to_i,from) )
#' where d_k(to_i|from) = 1( ne_k(to_i, from)>0) and ne_k(x,from) = sum_j 1(r[k-1] < ||to_i - from_j|| < r[k])
#'
#' @export
#' @useDynLib rgeyer


stepper_log_papangelou <- function(from, to, theta, r, sat, bbox, dbg=0, toroidal=0){
  rstepper_log_papangelou_c(
    from,
    to,
    theta,
    r,
    sat,
    bbox,
    dbg,
    toroidal)

}
