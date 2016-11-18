# testing of the basic functions of the c-classes
library(devtools)
load_all(".")

set.seed(2)

x00 <- matrix(runif(500),nc=2)

bbox <- cbind(0:1, 0:1)
iter <- 10
p <- 4
np <- p*(p-1)/2
m <- sample(0:(p-1), rep=T, nrow(x00))

x0 <- cbind(x00, m)

K <- 5
R <- 0.1
r <- 1:K/K * R

theta0 <- rep(2.1, p)

# generate randomm parameters
theta1 <- list()
for(i in 1:p) {
  K <- sample(1:10, 1)
  R <- runif(1, 0.06, 0.5)
  r <- 1:K/K * R
  theta1[[i]] <- list(r=r, theta = rnorm(K))
}

theta2 <- list()
for(i in 1:np) {
  K <- sample(1:10, 1)
  R <- runif(1, 0.06, 0.5)
  r <- 1:K/K * R
  theta2[[i]] <- list(r=r, theta = rnorm(K))
}

o <- rstepper_multi_c(theta0, theta1, theta2, bbox, iter = 1000, x0, 1000, 1)




## at least it runs
plot(o[,1:2], col=1+o[,3], asp=1, pch=1+o[,3])

