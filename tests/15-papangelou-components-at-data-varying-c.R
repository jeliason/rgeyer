# Test papangelou components various saturations, at data

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

out <-  stepper_components(xy, NULL, r = r, sat = cv0+10, bbox = bbox, dbg = 0, toroidal = toroidal)
# hard to check... check with simulations
