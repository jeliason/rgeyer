# testing of the basic functions of the c-classes
library(devtools)
load_all(".")

set.seed(2)
x0 <- matrix(runif(20),nc=2)
bbox <- cbind(0:1, 0:1)
iter <- 10
R <- 0.2
K <- 5
r <- (1:K)*R/K
theta <- c(  log(100), rep(.1, K))



out <- rtest_classes_c(theta, iter=1, toroidal =0, r = r, bbox = bbox, dbg = 1000, x0 = x0)
#

print(colSums(out[,3:(K+2)]))

plot(out[,1:2], asp=1, cex=0)

for(k in K:1){
  v <- (k/(K+1))
  symbols(out[,1], out[,2], circles = rep(r[k],nrow(out)), inches=F, add=T, bg=rgb(v,v,v,0.2), fg=NA)
}
text(out[,1:2], col="white", pch=3)
# ok
