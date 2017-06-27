# speeding up the components calculation
# Can't be sped up.
# double check the lambda calculations

  library(devtools)
  load_all(".")


  r <- c(.4, .5)
  sat <- c(1,1)
  K <- length(sat)
  n <- 2000
  set.seed(1)
  x <- matrix(runif(2*n), nc=2)
  u <- matrix(runif(2*n), nc=2)
  bb <- cbind(0:1, 0:1)

  if(!exists("a")){
    t0 <- system.time( a <- stepper_components(x, to=u, r=r, sat=sat, bbox = bb) )
    rr <- c(0,r)
    # by hand:
    # tf <- function(u, x, r1, r2, sat) {
    #   d <- sqrt(rowSums(t(t(x)-u)^2))
    #   min(sum(d < r2 & d >= r1), sat)
    # }
    #
    # b3 <- sapply(1:2, function(k) sapply(1:n, function(i){
    #   y <- rbind(x,u[i,])
    #   sum( sapply(1:nrow(y), function(j) tf(y[j,], y[-j,], rr[k],rr[k+1],sat[k]) ) )
    # }))
    # b4 <- sapply(1:2, function(k) sapply(1:n, function(i) sum( sapply(1:n, function(j) tf(x[j,], x[-j,], rr[k],rr[k+1],sat[k]) ) ) ))
    # b <- b3-b4

  #  print(cbind(a,b))
  }


t1 <- system.time( b <- stepper_components(x, to=u, r=r, sat=sat, bbox = bb) )

print(all.equal(a,b))


print(rbind(t0, t1))
