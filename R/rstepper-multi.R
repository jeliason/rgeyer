#' Simulate the multivariate Step-neighbour model
#'
#' Symmetric interaction weights
#'
#' @param theta0 first order parameters
#' @param theta1 intra-type interaction parameters. see Details
#' @param theta2 inter-type interaction parameters. see Details
#' @param n alternative to theta0: fixed point counts
#' @param bbox bounding box
#' @param iter iteration of bd
#' @param x0 the starting pattern location matrix (not necessary)
#' @param toroidal toroidal distances? 1=yes, 0=no
#' @param dbg dbg verbosity level.
#' @param trend optional list of im-objects to be used as trend components for each type.
#'
#' @details
#' Let's say we want p types.
#'
#' Then length(theta0) = p, length(theta1) = p, length(theta2) = p(p-1)/2.
#'
#' theta1 should be a list of lists, each of which specify the grid "r" and the heights "theta" of the stepfunctions, and a "c", non-negative integer vector of saturation levels.
#' "r" should not contain 0 as it will be added. The defaults have no interactions.
#'
#' theta2 should be like theta1, but now giving the parameters for each pair of types (i,j), i < j.
#'
#' x0, if given, should be a matrix cbind(x_coords, y_coords, marks) where marks are integers between 0 and p-1.
#'
#' @examples
#' theta0 <- log( c(200, 200, 200))
#' theta1 <- list(list(r= 0.02, theta = -1, c=1),
#'                list(r= 0.02, theta = -1, c=2),
#'                list(r= 0.02, theta = -1, c=3))
#' theta2 <- list(list(r= 0.01, theta=1, c=1),
#'                list(r= 0.01, theta=1, c=2),
#'                list(r= 0.02, theta=0, c=1))
#' bbox <- cbind(c(0,.5), c(-1,1))
#' x <- rstepper_multi(theta0, theta1, theta2, iter = 1e5, dbg=200, toroidal=T, bbox = bbox)
#' plot(x[,-3], col=x[,3]+1, pch=1+2*x[,3], cex=.5, asp=1)
#' # Fixed point counts
#' nvec <- c(200,200,200)
#' y <- rstepper_multi(n=nvec, theta1, theta2, iter = 1e5, dbg=200, toroidal=T, bbox = bbox)
#' plot(y[,-3], col=y[,3]+1, pch=1+2*y[,3], cex=.5, asp=1)
#'
#' @useDynLib rgeyer
#' @import Rcpp
#' @export

rstepper_multi <- function(n = c(100,100,100),
                           theta0 = rep(5, length(n)),
                           theta1 = rep(list(list(r = 0, theta = 0, c=0)), length(theta0)),
                           theta2 = rep(list(list(r = 0, theta = 0, c=0)), length(theta0)*(length(theta0)-1)/2),
                           bbox = cbind(0:1, 0:1),
                           iter = 10, x0 = NULL,
                           toroidal = 0, dbg = 0, trend = NULL) {
  n0 <- length(theta0)
  n1 <- length(theta1)
  n2 <- length(theta2)

  if(!is.null(n)) {
    nf <- length(n)
    if(nf != n1 | n2 != (nf*(nf-1)/2)) stop("check your parameters, length of n not the same as theta1.")
    n0 <- length(n)
    fixed <- TRUE
  }
  else{
    fixed <- FALSE
    if(n0 != n1 | n2 != (n0*(n0-1)/2) ) stop("check your parameters, lengths don't match.")
  }

  # check missing saturaionts for backward compatibility
  for(i in 1:n1){
    if(is.null(theta1[[i]][["c"]])) theta1[[i]][["c"]] <- rep(1, length(theta1[[i]][["r"]]))
  }
  for(i in 1:n2){
    if(is.null(theta2[[i]][["c"]])) theta2[[i]][["c"]] <- rep(1, length(theta2[[i]][["r"]]))
  }
  #browser()

  # check trend
  if(is.null(trend)){
    trendl <- rep(list(list()), n0)
  }
  else{ # just one image for all
    if(is.im(trend)){
      trendl <- rep(list(trend), n0)[1:n0]
    }
    else if(is.list(trend)){
        if(any(!sapply(trend, is.im)))
          stop("trend should be a list of im-objects, one for each type")
        if(length(trend) < n0) stop("trend should be a list of im-objects, one for each type")
        trendl <- trend
    }
    else stop("can't interpret trend")
  }
  #browser()
  # generate starting pattern
  if(is.null(x0)){
    V <- prod( apply(bbox, 2, diff))
    if(fixed){
      nstart <- sum(n)
      x0 <- apply(bbox, 2, function(a) runif(nstart, a[1], a[2]))
      m0 <- rep(1:n0-1, n)
    }
    else{
      nstart <- max(5 * n0,  rpois(1, sum(exp(theta0))/V  )   )
      x0 <- apply(bbox, 2, function(a) runif(nstart, a[1], a[2]))
      m0 <- sample(1:n0-1, nstart, replace=T, prob = exp(theta0))
    }
    x0 <- cbind(x0, m0)
  }
  else{
    m0 <- x0[,3]
  }
  # mm <- table(m0)
  # #browser()
  # if(dbg>10) cat(sprintf("starting pattern, n=%i \n", sum(mm)))
  if(fixed) out <- rstepper_multi_fixed_c(theta1, theta2, bbox, iter, x0, dbg, toroidal, trendl)
  else out <- rstepper_multi_c(theta0, theta1, theta2, bbox, iter, x0, dbg, toroidal, trendl)
  out
}
