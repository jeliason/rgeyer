# Test bivariate simulator

library(devtools)
load_all(".")

set.seed(40)
L <- 5
bbox <- cbind(c(0,L), c(0,L))
R <- 0.2
K <- 3
r <- (1:K)*R/K
theta0 <- log(30/L^2)
theta <- c(theta0, theta0, 1*(K:1)/K)
iter <- 1e5

xx0 <- matrix(runif(400,0, L), nc=2)
prob <- exp(theta[1])/(exp(theta[1])+exp(theta[2]))
m0 <- rbinom(nrow(xx0), 1, prob)
x0 <- cbind(xx0, m0)

# out <-  rstepper_biv_c(theta,
#                         r,
#                         bbox,
#                         iter,
#                         x0,
#                         dbg = 1,
#                         toroidal = 1)

out <- rstepper_biv(theta, R, K, bbox, iter = iter, dbg = 100)

plot(out[,1:2], asp=1, cex=0)
symbols(out[,1], out[,2], inches=F, add=T, circles = rep(R, nrow(out)), fg=NA, bg = rgb(1-out[,3], 0, out[,3], 0.2))
points(out[,1:2], pch=19, col=2+2*out[,3])
print(table(out[,3]))
print(colSums(out[,-c(1:3)]))


# ok!
