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
xs <- seq(0, L, l = 50)
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

if(!exists("xym")) xym <- rstepper_multi(n=rep(20, p), theta1=theta1, theta2=theta2, bbox=bbox, iter=1e4)
xym <- xym[order(xym[,1]),]
xy <- xym[,-3]
m <- xym[,3]

## and at datapoints.
vd <- stepper_multi_log_papangelou(xym, NULL, theta1, theta2, bbox, dbg=111, toroidal = 0)
vd <- round(vd,1)


dbg <- 0

# on the grid, faster
t1 <- Sys.time()
vl2m <- stepper_multi_log_papangelou(xym, grid, theta1, theta2, bbox, dbg=dbg, toroidal = 0, multi_to = TRUE)
t1 <- Sys.time() - t1
vl2 <- split(c(vl2m), rep(1:p, each=nrow(vl2m)))


# on the grid
vl <- list()
t0 <- Sys.time()
for(i in 1:p)
  vl[[paste0(i-1)]] <- stepper_multi_log_papangelou(xym, cbind(grid,i-1), theta1, theta2, bbox, dbg=dbg, toroidal = 0)
t0 <- Sys.time()-t0
vlm <- do.call(cbind, vl)


# check
if(p < 4){
par(mfrow=c(3,p+1))
for(i in 1:p){
  image(xs, xs, matrix(vl[[i]], length(xs)), asp=1, col=gray.colors(1000), main=i)
  points(xy, pch=m*2+1)
}
plot(xy, pch=m*2+1, cex=abs(vd)/max(abs(vd))+1, col = 3-(vd<0)+(vd==0), asp=1)
for(i in 1:K) symbols(xy, circles=rep(r[i],nrow(xy)), inches=F, add=T, fg="gray95")
text(xy-.15, labels = vd, cex=.8)
for(i in 1:p){
  image(xs, xs, matrix(vl2[[i]], length(xs)), asp=1, col=gray.colors(1000), main=i)
  points(xy, pch=m*2+1)
}
plot(xy, pch=m*2+1, cex=abs(vd)/max(abs(vd))+1, col = 3-(vd<0)+(vd==0), asp=1)

for(i in 1:p){
  image(xs, xs, matrix(vl2[[i]]-vl[[i]], length(xs)), asp=1, col=gray.colors(1000), main=i)
  points(xy, pch=m*2+1)
}
}

print(all.equal(vlm, vl2m) )
print(rbind( t0, t1) )


