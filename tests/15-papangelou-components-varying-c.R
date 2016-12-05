# Test papangelou components various saturations

library(devtools)
load_all(".")


set.seed(1)
L <- 10
xy <- matrix(runif(40, 0, L),ncol=2)
xs <- seq(0, L, l = 140)
bbox <- cbind(c(0,L), c(0,L))
grid <- as.matrix(expand.grid(xs, xs))
R <- 2
K <- 3
r <- R/K * (1:K)
toroidal <- 1
cv0 <- rep(0, K)

asm <- function(k) matrix(out[,k], ncol=length(xs))
par(mfrow=c(2,3))

for(i in c(1,4)){
out <-  stepper_components(xy, grid, r = r, sat = cv0+i, bbox = bbox, dbg = 0, toroidal = toroidal)
  for(k in 1:K) {
    image(xs,xs,asm(k), col=gray.colors(120, 0), zlim=c(0, max(out)), asp=1, main=paste0("c=",i,": k=",k))
    points(xy, pch=3, col="white")
    symbols(xy, inches=F, add=T, circles=rep(r[k], nrow(xy)), fg=rgb(0,0,0,0.5))
  }
}
# ok
