# Test basic simulator for the relative stepper, univariate

library(devtools)
load_all(".")

set.seed(40)
L <- 2
bbox <- cbind(c(0,L), c(0,L))
R <- L/20
K <- 3
r <- (1:K)*R/K
theta <- c(-1,1,.5)*1 #(K:1) / K
iter <- 1e5


# Free n
if(0){
theta0 <- log(300/prod(apply(bbox,2,diff)))

xo <- rrelativestepper(theta=c(theta0, theta), r = r, bbox = bbox, dbg=11, iter = iter, toroidal = T)
}

# fixed n
if(1){
n <- 229
xo <- rrelativestepper(theta=theta, n=n, r = r, bbox = bbox, dbg=11, iter = iter, toroidal = T)
}


par(mfrow=c(2,1))
plot(xo, asp=1, cex=.3, xlim=bbox[,1], ylim=bbox[,2])
plot(xo, asp=1, cex=.1, xlim=bbox[,1], ylim=bbox[,2])
for(i in 1:K) symbols(xo, circles=rep(r[i],nrow(xo))/2, inches=F, add=T, fg="gray90")

# crap
