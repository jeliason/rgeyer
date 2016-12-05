# Check univ simulator, varying c

library(devtools)
load_all(".")
r <- c(.1,.5)
bbox <- cbind(0:1,0:1)*10
theta <- c(-1, 1)
cv <- c(1,1)*10

# Simulate free n
if(0){
x <- rstepper(theta = c(log(5), theta), r = r, sat = cv, dbg=100, iter=1e5, bbox=bbox, toroidal=F)
print(table(x[,3]))
}
if(1){ # fixed n
x <- rstepper(theta = theta, n=400, r = r, sat = cv, dbg=100, iter=1e6, bbox=bbox, toroidal=F)
}
#

plot(x[,1:2], asp=1, cex=.2)
