# test the rstepper papangelou evaluation

library(devtools)
load_all(".")




p <- 5
P <- p*(p-1)/2
n <- 200

x <- matrix(runif(2 * n), nc=2)
m <- sample(1:p, n, replace=T)

rv <- c(0.1, 0.2)
the <- c(0.1, 0.05)
sat <- c(1,1)

theta0 <- rep(0.1, p)
theta1 <- rep(list(list(r=rv, theta=the, c=sat)), p)
theta2 <- rep(list(list(r=rv, theta=the, c=sat)), P)


rstepper_papangelou(x, m, theta0=theta0, theta1 = theta1, theta2 = theta2)
