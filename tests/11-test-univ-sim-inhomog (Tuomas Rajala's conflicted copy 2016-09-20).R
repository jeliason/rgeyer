# Test basic simulator with intensity given by image.

library(devtools)
load_all(".")

set.seed(40)
L <- 5
bbox <- cbind(c(0,L*2), c(0,L))
x0 <- matrix(runif(200,0, L),nc=2)
R <- 0.2
K <- 5
r <- (1:K)*R/K
theta <- -(K:1)
iter <- 1e5
n <- 200

# trend
z <- runifpoint(5, win = as.owin(c(bbox)))
trend <- density(z, adjust = 1)
trend <- trend/max(trend) * 10

t0 <- system.time(
  out <-  rstepper(theta = theta, n = n,
                        r = r,
                        bbox = bbox,
                        iter = iter,
                        dbg = 11,
                        toroidal = 1)
)
# with trend
t1 <- system.time(
  out2 <-  rstepper(theta = theta, n = n,
                   r = r,
                 bbox = bbox,
                 iter = iter,
                 dbg = 11,
                 trend = trend,
                 toroidal = 1)
)

print(rbind(t0,t1))

par(mfrow=c(2,1))
plot(trend*0, col=gray.colors(120, .4,.8))
points(out[,1:2])
plot(trend, col=gray.colors(120, .4,.8), ribbon=F)
points(out2[,1:2])
