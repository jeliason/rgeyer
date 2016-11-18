# test multi-simulator saturation
library(devtools)
load_all(".")

# check basic poisson
if(0){
  theta0 <- log(c(50, 50, 300, 10))
  bbox <- cbind(c(-1,.4), c(0,1.5))
  V <- prod(apply(bbox,2, diff))
  x <- rstepper_multi(theta0, dbg = 1100, iter=2e4, bbox = bbox)
}
# fixed n
if(0){
  n <- c(50, 50, 300, 10)
  bbox <- cbind(c(-1,.4), c(0,1.5))
  x <- rstepper_multi(n=n, dbg = 1100, iter=2e4, bbox = bbox)
}

# fixed simplest bivariate
if(1){
  n <- c(150, 150)
  bbox <- cbind(c(-1,.4), c(0,1.5))
  theta1 <- list(list(r=0.2, theta=1, c=3), list(r=0.1, theta=-1))
  x <- rstepper_multi(n=n, dbg = 1100, iter=2e4, bbox = bbox, theta1=theta1)
  plot(x,col=x[,3]+1, asp=1, cex=.5)
}
