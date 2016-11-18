#' Simulate the univariate Step-Geyer model
#'
#' Simulate Geyer's saturation with K different annulus and saturation fixed to 1.
#'
#' @param theta weights, length K+1, first element controls intensity
#' @param n alternative to theta0: fixed point counts. Then theta should be K-length. see details.
#' @param R maximum intereaction
#' @param K number of steps
#' @param r alternative to R,K. Should be increasing, not containing 0. Length will be taken as K.
#' @param bbox bounding box
#' @param iter iteration of bd
#' @param x0 the starting pattern location matrix
#' @param toroidal toroidal distances?
#' @param trend Optional im-object defining a trend.
#' @param dbg dbg
#'
#' @details
#' If r not given, the steps will be (1:K)*R/K. 0 will be the beginning of the first interval, so theta[2] is mapped to range (0,r[1]).
#' If trend is given, it should be a 'spatstat' im-object covering the bbox. Will not be checked (at the moment).
#' @useDynLib rgeyer
#' @import Rcpp
#' @export

rstepper <- function(theta, R, K, r, n=NULL,
                     bbox = cbind(0:1, 0:1), iter = 10,
                     x0 = NULL, toroidal = 0, dbg = 0, trend = NULL) {
  # check range vector
  if(missing(r)) r <- (1:K) * R / K
  #
  K <- length(r)
  if(K==1) if(any(diff(r))<0) stop("r should increasing.")
  #
  n1 <- length(theta)
  if(!is.null(n)){
    if(length(n) != 1) stop("n should be single integer > 0")
    if(n1 != K) stop("n given: length of theta should be length of r (=K)")
    fixed <- TRUE
  }
  else{
    fixed <- FALSE
    if((K+1) != n1 ) stop("Length of theta should be [length of r (=K) + 1] ")
  }

  # check trend
  if(!is.null(trend)){
    if(!is.im(trend)) stop("trend should be im-object")
  }
  else{
    trend <- list()
  }

  # generate starting pattern
  if(is.null(x0)){
    V <- prod( apply(bbox, 2, diff) )
    if(fixed){
      nstart <- n
      x0 <- apply(bbox, 2, function(a) runif(nstart, a[1], a[2]))
    }
    else{
      nstart <- max(10,  rpois(1,exp(theta[1])/V  )   )
      x0 <- apply(bbox, 2, function(a) runif(nstart, a[1], a[2]))
    }
  }
  if(fixed) out <- rstepper_univ_fixed_c(n, theta, r, bbox, iter, x0, dbg, toroidal, trend)
  else out <- rstepper_univ_c(theta, r, bbox, iter, x0, dbg, toroidal, trend)

  out
}
