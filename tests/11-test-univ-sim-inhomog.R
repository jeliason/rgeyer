# Test basic simulator with intensity given by image.

library(devtools)
load_all(".")

set.seed(40)
L <- 2
bbox <- cbind(c(0,L*2), c(0,L))
W <- as.owin(c(bbox))
nc <- 128
nr <- nc/2
x0 <- matrix(runif(200,0, L), nc=2)
R <- L/10
K <- 5
r <- (1:K)*R/K
theta <- -(K:1)
iter <- 1e5
n <- 100

# trend
#z <- runifpoint(15, win = as.owin(c(bbox)))
#trend <- density(z, adjust = 1)
#trend <- trend/max(trend) * 10

trend <- as.im(function(x,y) exp(x), W = W, dimyx = c(nr,nc))
trend <- trend/2

# t0 <- system.time(
#   out <-  rstepper(theta = theta, n = n,
#                         r = r,
#                         bbox = bbox,
#                         iter = iter,
#                         dbg = 11,
#                         toroidal = 1)
# )
# with trend
t1 <- system.time(
  out2 <-  rstepper(theta = theta, n = n,
                   r = r,
                 bbox = bbox,
                 iter = iter,
                 dbg = 1,
                 trend = trend,
                 toroidal = 0)
)

#print(rbind(t0,t1))

par(mfrow=c(2,1))
#plot(trend*0, col=gray.colors(120, .4,.8))
#points(out[,1:2])
plot(trend, col=gray.colors(120, .4,.8), ribbon=F)
points(out2[,1:2])
plot(density(ppp(out2[,1], out2[,2], window = W), adjust=2))
