#' Compute the papangelou conditional intensity for the stepper model
#'
#' @param from coordinates of the data
#' @param to coordinate matrix for locations for evaluation of papangelou. Assumed to be different from 'from'
#' @param theta vector of weights, K+1 long
#' @param r steps, increasing, does not include 0
#' @param bbox bounding box
#' @param dbg dbg level
#' @param toroidal >0 use toroidal distances
#'
#' @return log-conditional intensity l(to_i|from) = theta[0] + sum( theta[k] * d_k(to_i,from) )
#' where d_k(to_i|from) = 1( ne_k(to_i, from)>0) and ne_k(x,from) = sum_j 1(r[k-1] < ||to_i - from_j|| < r[k])
#'
#' @export
#' @useDynLib rgeyer


stepper_log_papangelou <- function(from, to, theta, r, bbox, dbg=0, toroidal=0){
  rstepper_log_papangelou_c(
    from,
    to,
    theta,
    r,
    bbox,
    dbg,
    toroidal)

}
