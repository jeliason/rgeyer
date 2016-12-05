# Test papangelou with different c

library(devtools)
load_all(".")


set.seed(1)
L <- 10
xy <- matrix(runif(600, 0, L),ncol=2)
xs <- seq(0, L, l = 128)
bbox <- cbind(c(0,L), c(0,L))
grid <- as.matrix(expand.grid(xs, xs))
R <- 1
K <- 3
r <- 1:K * R/K
theta <- c(1, K:1)

toroidal <- 1
cv1 <- rep(1, K)


par(mfrow=c(2,2))
for(i in c(0,2,5,10)){
  out <- stepper_log_papangelou(xy, grid, theta, r, bbox, toroidal = toroidal, sat = cv1+i)
  M <- matrix(out, ncol=length(xs))
  image(xs, xs, M, asp=1, col=gray.colors(1000,0), main=i+1, zlim = c(0,50))
  points(xy, pch=3, cex=0.5)
}

