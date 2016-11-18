# Test papangelou components

library(devtools)
load_all(".")


set.seed(1)
L <- 10
xy <- matrix(runif(40, 0, L),ncol=2)
xs <- seq(0, L, l = 200)
bbox <- cbind(c(0,L), c(0,L))
grid <- as.matrix(expand.grid(xs, xs))
R <- 2
K <- 3
r <- R/K * (1:K)
toroidal <- 1

out <-  stepper_components(
  xy, grid,
  R, K,
  bbox,
  dbg = 0,
  toroidal = toroidal)


par(mfrow=c(ceiling(K/2), 2))
asm <- function(k) matrix(out[,k], ncol=length(xs))

for(k in 1:K) {
  image(xs,xs,asm(k), col=gray.colors(120, 0), zlim=c(0, max(out)), asp=1)
  points(xy, pch=3, col="white")
  symbols(xy, inches=F, add=T, circles=rep(r[k], nrow(xy)), fg=rgb(0,0,0,0.5))
}
# ok

##

outx <-  stepper_components(
  xy, to = NULL,
  R, K,
  bbox,
  dbg = 0,
  toroidal = toroidal)

k <-2
image(xs,xs,asm(k), col=gray.colors(120, 0), zlim=c(0, max(out)), asp=1)
points(xy, pch=3, col=-outx[,k]+1)
symbols(xy, inches=F, add=T, circles=rep(r[k], nrow(xy)), fg=rgb(0,0,0,0.5))




