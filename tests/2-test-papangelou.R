# Test papangelou

library(devtools)
load_all(".")


set.seed(1)
L <- 10
xy <- matrix(runif(200, 0, L),ncol=2)
xs <- seq(0, L, l = 128)
bbox <- cbind(c(0,L), c(0,L))
grid <- as.matrix(expand.grid(xs, xs))
R <- 1
K <- 3
r <- 1:K * R/K
theta <- c(1, K:1)

toroidal <- 1

out <-  rstepper_log_papangelou_c(
  xy, grid, theta,
  r,
  bbox,
  dbg = 0,
  toroidal = toroidal)


M <- matrix(out, ncol=length(xs))


out2 <- stepper_log_papangelou(xy, grid, theta, r, bbox, toroidal = toroidal)
M2 <- matrix(out2, ncol=length(xs))

par(mfrow=c(2,1))
image(xs, xs, M, asp=1, col=gray.colors(1000,0), zlim=c(0, max(M)))
points(xy, pch=3)

image(xs, xs, M2, asp=1, col=gray.colors(1000,0), zlim=c(0, max(M)))
points(xy, pch=3)
