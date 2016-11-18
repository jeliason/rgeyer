# Test multiv simulator with intensity given by image.
library(spatstat)
library(devtools)
load_all(".")

set.seed(2)
L <- 5
bbox <- cbind(c(0,2*L), c(0,L))
W <- as.owin(c(bbox))
# trend
nr <- 100
nc <- 100

trend1 <- as.im(function(x,y) exp(x/(2*L)), W = W, dimyx = c(nr,nc))
trend2 <- 2-trend1
trend <- list(trend1, trend2, trend1*0)


iter <- 1e4
p <- 3
np <- p*(p-1)/2

K <- 5
n <- rep(100, p)

# generate randomm parameters
theta1 <- list()
for(i in 1:p) {
  K <- sample(1:10, 1)
  R <- runif(1, L*0.002, L*0.1)
  r <- 1:K/K * R
  theta1[[i]] <- list(r=r, theta = -abs(rnorm(K)))
}

theta2 <- list()
for(i in 1:np) {
  K <- sample(0:10, 1)
  R <- runif(1, L*0.002, L*0.05)
  r <- if(K==0) 0 else 1:K/K * R
  theta2[[i]] <- list(r=r, theta = -abs(rnorm(K)))
}
t0 <- system.time(
o <- rstepper_multi(n=n, theta1=theta1, theta2=theta2,
                    bbox=bbox, iter = iter, dbg=1000, tor=0, trend = trend)
)

print(t0)


## at least it runs
par(mfrow=c(3,1))
for(i in 1:3){
  plot(trend[[i]], col=gray.colors(120,.4,.7))
  points(o[o[,3]==i-1,-3])
}

