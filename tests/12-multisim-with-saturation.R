# test multi-simulator saturation varition
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
  n <- c(150, 10)
  bbox <- cbind(c(-1,.4), c(0,1.5))
  theta1 <- list(list(r=c(.02), theta=-5, c=1), list(r=0.1, theta=-100, c = 1))
  theta2 <- list(list(r=c(.02, 0.1), theta=c(-5, 10), c=c(20,1)))
  x <- rstepper_multi(n=n, dbg = 11, iter=1e5, bbox = bbox, theta1=theta1, theta2 = theta2)
  plot(x,col=x[,3]+1, asp=1, cex=.5, xlim=bbox[,1], ylim=bbox[,2])
  for(i in 1:2) symbols(x[x[,3]==1,], circles=rep(theta2[[1]]$r[i], sum(x[,3])), inches=F, add=T, col="gray90")
}
