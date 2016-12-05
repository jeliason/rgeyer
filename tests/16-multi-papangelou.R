#' Test multi-stepper papangelou

library(devtools)
load_all(".")


set.seed(1)
L <- 10
# mark
p <- 3
# m <- sample(0:(p-1), replace=T, nrow(xy))
# xym <- cbind(xy, m)

# grid
xs <- seq(0, L, l = 150)
bbox <- cbind(c(0,L), c(0,L))
grid <- as.matrix(expand.grid(xs, xs))
gridm <- cbind(grid, 0)
gridm1 <- cbind(grid, 1)
# parameters
R <- .8
K <- 1
r <- R/K * (1:K)
cv0 <- rep(0, K)
# the parameters
theta1 <- rep(list(list(r=r, theta=c(1,.5,.2)[1:K], c=cv0+1)), p)
theta2 <- rep(list(list(r=r, theta = c(-1,.2,0)[1:K], c=cv0+1)), p*(p-1)/2)

xym <- rstepper_multi(n=rep(20, p), theta1=theta1, theta2=theta2, bbox=bbox, iter=1e4)
xym <- xym[order(xym[,1]),]
xy <- xym[,-3]
m <- xym[,3]

## and at datapoints.
vd <- stepper_multi_log_papangelou(xym, NULL, theta1, theta2, bbox, dbg=111, toroidal = 0)
vd <- round(vd,1)

# on the grid
vl <- list()
for(i in 1:p)
  vl[[i]] <- stepper_multi_log_papangelou(xym, cbind(grid,i-1), theta1, theta2, bbox, dbg=1, toroidal = 0)



# check
par(mfrow=c(2,p/2+1))
for(i in 1:p){

  image(xs, xs, matrix(vl[[i]], length(xs)), asp=1, col=gray.colors(1000), main=i)
  points(xy, pch=m*2+1)
}

plot(xy, pch=m*2+1, cex=abs(vd)/max(abs(vd))+1, col = 3-(vd<0)+(vd==0), asp=1)
for(i in 1:K) symbols(xy, circles=rep(r[i],nrow(xy)), inches=F, add=T, fg="gray95")
text(xy-.15, labels = vd, cex=.8)

