# Check fixed n simulator

library(devtools)
load_all(".")
r <- c(.1,.15,.25)
bbox <- cbind(0:1,0:1)*10
theta1 <- rep(list(list(r=2*r,   theta = -2 *c(1,1,1) )), 3)
theta2 <- rep(list(list(r=r, theta = 0.5*c(1,1,1) )), 3)


# Simulate
x <- rstepper_multi(n = c(100,100,100), theta1=theta1,
                    theta2 = theta2,
                    dbg=100, iter=1e5, bbox=bbox, toroidal=T)
(table(x[,3]))

#
## Compute pcfs
library(Kcross)
library(spatstat)
p <- ppp(x[,1],x[,2], window = as.owin(c(bbox)), marks=factor(x[,3]))
rv <- seq(0, 2, l = 50)
g <- pcf_cross_all_box(p, r = rv, adjust=1.5)




# check
par(mfrow=c(2,1))
plot(x, col=x[,3]+1, asp=1, pch=x[,3]+1, cex=.5)
plot(rv,g[1,1,], ylim=c(0,max(g)), type="l")
lines(rv, g[2,2,], col=2)
lines(rv, g[3,3,], col=3)
lines(rv, g[1,2,], lty=2)
lines(rv, g[1,3,], col=2, lty=2)
lines(rv, g[2,3,], col=3, lty=2)
abline(h=1)
