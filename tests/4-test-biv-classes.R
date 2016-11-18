# testing of the basic functions of the c-classes
library(devtools)
load_all(".")

set.seed(1)
xx0 <- matrix(runif(30),nc=2)
m0 <- sample(0:1, nrow(xx0), replace = T)
x0 <- cbind(xx0, m0)
bbox <- cbind(0:1, 0:1)
iter <- 10
R <- 0.2
K <- 3
r <- (1:K)*R/K
theta <- rep(0, K+2)

out <- rtest_classes_biv_c(theta, r, bbox, iter, x0, dbg=11, toroidal=F)

#

print(colSums(out[,3:(K+3)]))

plot(out[,1:2], asp=1, cex=0, xlim=c(0,1), ylim=c(0,1))

for(k in K:1){
  v <- (k/(K+1))
  v1 <- v * (out[,3] == 0)
  v2 <- v * (out[,3] == 1)
  symbols(out[,1], out[,2], circles = rep(r[k],nrow(out)), inches=F, add=T, bg=rgb(v1,0,v2,0.2), fg=NA)
}
text(out[,1:2], col="white", pch=3)
# ok
