# Test papangelou components

library(devtools)
load_all(".")


set.seed(1)
L <- 10
p <- 3
xy <- matrix(runif(2*p*20, 0, L),ncol=2)
# mark
m <- sample(0:1, replace=T, nrow(xy))
xym <- cbind(xy, m)

xs <- seq(0, L, l = 200)
bbox <- cbind(c(0,L), c(0,L))
grid <- as.matrix(expand.grid(xs, xs))
R <- 2
K <- 3
r <- R/K * (1:K)
toroidal <- 1

gridm <- cbind(grid, 0)

out <-  stepper_biv_components(
  xym, gridm, r = r, sat = rep(1,K),
  bbox=bbox,
  dbg = 0,
  toroidal = toroidal)


par(mfrow=c(ceiling(K/2), 2))
asm <- function(k) matrix(out[,k], ncol=length(xs))

for(k in 1:K) {
  image(xs,xs,asm(k), col=gray.colors(120, 0), zlim=c(0, max(out)), asp=1)
  points(xy, pch=1+m, col=c("cyan","magenta")[1+m])
  symbols(xy, inches=F, add=T, circles=rep(r[k], nrow(xy)), fg=rgb(m,.5,1-m,0.5))
}
# ok

# ##
#
# outx <-  stepper_components(
#   xy, to = NULL,
#   R, K,
#   bbox,
#   dbg = 0,
#   toroidal = toroidal)
#
# k <-2
# image(xs,xs,asm(k), col=gray.colors(120, 0), zlim=c(0, max(out)), asp=1)
# points(xy, pch=3, col=-outx[,k]+1)
# symbols(xy, inches=F, add=T, circles=rep(r[k], nrow(xy)), fg=rgb(0,0,0,0.5))
#
#
#
#
