#' Simulate the bivariate Step-Geyer model
#'
#' Symmetric interaction weights
#'
#' @param theta weights, length K+2, first element controls intensity
#' @param R maximum intereaction
#' @param K number of steps
#' @param r alternative to R,K. Should be increasing, not containing 0.
#' @param bbox bounding box
#' @param iter iteration of bd
#' @param x0 the starting pattern location matrix
#' @param toroidal toroidal distances?
#' @param dbg dbg
#'
#' @details
#' the steps will be (1:K)*R/K. 0 will be the beginning of the first interval, so theta[2] is mapped to range (0,r[1]).
#' @useDynLib rgeyer
#' @import Rcpp
#' @export

rstepper_biv_obs <- function(theta, R, K, r, bbox = cbind(0:1, 0:1), iter = 10, x0 = NULL, toroidal = 0, dbg = 0) {
  if(missing(r)) r <- (1:K) * R / K
  K <- length(r)
  if(length(theta) != K+2) stop("theta should be of length K+2")

  if(is.null(x0)){
    V <- prod( apply(bbox, 2, diff))
    n <- max(5, rpois(1, (exp(theta[1])+exp(theta[2])) * V))
    x0 <- apply(bbox, 2, function(a) runif(n, a[1], a[2]))
    prob <- exp(theta[2])/ sum(exp(theta[1:2]))
    m0 <- rbinom(n, 1, prob)
    x0 <- cbind(x0, m0)
  }
  mm <- table(m0)
  if(dbg>10) cat(sprintf("starting pattern: (%i, %i)\n", mm[1],mm[2]))
  out <- rstepper_biv_c(theta, r, bbox, iter, x0, dbg, toroidal)

  out
}
