# Test papangelou components varying c

library(devtools)
load_all(".")


set.seed(1)
L <- 10
xy <- matrix(runif(20, 0, L),ncol=2)
# mark
m <- sample(0:1, replace=T, nrow(xy))
xym <- cbind(xy, m)

xs <- seq(0, L, l = 200)
bbox <- cbind(c(0,L), c(0,L))
grid <- as.matrix(expand.grid(xs, xs))
R <- 3
K <- 3
r <- R/K * (1:K)
toroidal <- 1
cv0 <- rep(0, K)
gridm <- cbind(grid, 0)

# out <-  stepper_biv_components(
#   xym, gridm,
#   r=r,
#   sat=cv0 + 1,
#   bbox = bbox,
#   dbg = 0,
#   toroidal = toroidal)
#
# par(mfrow=c(ceiling(K/2), 2))
# asm <- function(k) matrix(out[,k], ncol=length(xs))
#
# for(k in 1:K) {
#   image(xs,xs,asm(k), col=gray.colors(120, 0), zlim=c(0, max(out)), asp=1)
#   points(xy, pch=1+m, col=c("cyan","magenta")[1+m])
#   symbols(xy, inches=F, add=T, circles=rep(r[k], nrow(xy)), fg=rgb(m,.5,1-m,0.5))
# }
# ok
out1 <-  stepper_biv_components(
  xym, NULL,
  r=r,
  sat=cv0 + 1,
  bbox = bbox,
  dbg = 0,
  toroidal = toroidal)

out2 <-  stepper_biv_components(
  xym, NULL,
  r=r,
  sat=cv0 + 2,
  bbox = bbox,
  dbg = 0,
  toroidal = toroidal)

xyv <- xy

par(mfrow=c(2,2))
for(k in 1:K) {
  plot(xyv, cex=0)#, xlim=c(5,10), ylim=c(3,9))
  # plot annulus
  r1 <- c(0,r)[k]
  d <- r[k]-r1
  for(i in 1:100){
    symbols(xyv, inches=F, add=T, circles=rep(r1+(i-1)/100*d, nrow(xyv)), fg=gray(.9))
  }
  symbols(xyv, inches=F, add=T, circles=rep(r1, nrow(xyv)), fg=gray(.7))
  symbols(xyv, inches=F, add=T, circles=rep(r1+d, nrow(xyv)), fg=gray(.7))

  #
  text(xy, labels = 1:nrow(xy), col=m+1)
}


