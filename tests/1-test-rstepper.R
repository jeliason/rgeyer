# Test basic simulator

library(devtools)
load_all(".")

set.seed(40)
L <- 5
bbox <- cbind(c(0,L), c(0,L))
x0 <- matrix(runif(200,0, L),nc=2)
iter <- 10
R <- 0.2
K <- 5
r <- (1:K)*R/K
theta <- c(  log(100/L^2), 0.7*(K:1)/K)
iter <- 1e5


out <-  rstepper_univ_c(theta,
                        r,
                        bbox,
                        iter,
                        x0,
                        dbg = 0,
                        toroidal = 1)



plot(out[,1:2], asp=1)
