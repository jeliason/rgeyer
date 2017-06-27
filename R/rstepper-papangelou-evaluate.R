#' Compute the papangelou conditional intensity for multi-type multi-range saturation model
#'
#'
#' @param x the point pattern data location matrix
#' @param m mark vector, converted to 0-based integers by as.integer(factor(m))-1
#' @param at which ones to compute: "data": at data points, "grid": at grid locations
#' @param grid_dimyx if at="grid", the (nrow, ncol) dimensions of the grid
#' @param bbox if at="grid", the bounding box for the grid
#' @param theta0 intensity parameters, see Details.
#' @param theta1 intra-type parameters
#' @param theta2 inter-type parameters
#' @param covariates an optional list of im-objects
#'
#' @details
#' The parametrisation is much like in \link{rstepper-multi}. Let's say nlevels(factor(m)) = p, and length(covariates)=k:
#'
#' theta0: list(c(intercept_1, cov_11, cov_12, ..., cov_1k), ..., c(intercept_p, ..., cov_pk)) with individual coefficients being real numbers.
#'
#' theta1: list(list(theta1_1, r1_1, c1_1), ..., c(theta1_p, r1_p, c1_p)) where theta1_i are the coefficients for the step functions defined by ranges r1_1, saturation levels per step in c1_i. So length(theta1) = p.
#'
#' theta2: like theta1 for each pair i=(l,k), l < k in 1,...,p. Length(theta2) = p*(p-1)/2.
#'
#' covariates: list of im-objects, length(covariates) = k. If missing, the cov_il coefficients in theta0 should be omitted.
#'
#' If covariates are not given, theta0 can be a simple vector of length p.
#'
#' @export
#' @import spatstat
#' @useDynLib rgeyer

rstepper_papangelou <- function(x, m, at = "data", grid_dimyx, bbox, theta0, theta1, theta2, covariates) {
  # check input:
  covs_given <- !missing(covariates)
  if(covs_given) {
    if(!is.list(covariates)) stop("covariates should be a list")
    if(!all(sapply, covariates, is.im)) stop("covariate elements should be im-objects.")
    ncovs <- length(covariates)
  }
  else ncovs <- 0
  if(!is.matrix(x)) stop("x should be coordinate matrix")
  n <- nrow(x)
  d <- ncol(x)
  m <- as.integer(factor(m)) - 1
  p <- length(unique(m))
  if(length(m) != n) stop("nrow(x) != length(m)")
  #
  l0 <- length(theta0)
  if(l0 != p) stop("length(theta0) != nlevels(m)")
  if(ncovs){
    if(!is.list(theta0)) stop("covariates given, theta0 should be a list")
  }
  else{
    covariates <- list()
    if(!is.list(theta0)) theta0 <- as.list(theta0)
  }
  if(all(sapply(theta0, length) != (ncovs + 1))) stop("theta0 elements wrong size")
  #
  l1 <- length(theta1)
  if(l1 != p) stop("length(theta1) != nlevels(m)")
  if(!is.list(theta1)) stop("theta1 should be a list")
  req <- c("theta", "r", "c")
  if(!all(sapply(theta1, function(z) all(names(z)%in%req) )))
    stop("theta1 elements should be list(theta=c(), r=c(), c=c())")
  if(!all( sapply(theta1, function(z) all(sapply(z, length)==length(z$r)) ) ) )
    stop("theta1 elements should be list(theta=c(), r=c(), c=c()) with equal length vectors")
  #
  l2 <- length(theta2)
  P <- p*(p-1)/2
  if(l2 != P) stop("length(theta2) != nlevels(m)*(nlevels(m)-1)/2")
  if(!is.list(theta2)) stop("theta2 should be a list")
  req <- c("theta", "r", "c")
  if(!all(sapply(theta2, function(z) all(names(z)%in%req) )))
    stop("theta2 elements should be list(theta=c(), r=c(), c=c())")
  if(!all( sapply(theta1, function(z) all(sapply(z, length)==length(z$r)) ) ) )
    stop("theta2 elements should be list(theta=c(), r=c(), c=c()) with equal length vectors")
  #

  #
  if(at == "data") {
    out <- log_papangelou_at_data(x, m, theta0, theta1, theta2, covariates)
  }
  else if(at == "grid"){
    if((missing(grid_dimyx)|missing(bbox)))
      stop("at='grid' requires that also grid_dimyx and bbox are given.")
    grid <- cbind(seq(bbox[1,1], bbox[2,1], l = grid_dimyx[2]),
                  seq(bbox[1,2], bbox[2,2], l = grid_dimyx[1]))
    out <- log_papangelou_at_other(x, m, theta0, theta1, theta2, covariates, grid)
  }
  # done.
  out
}




