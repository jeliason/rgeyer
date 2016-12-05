# Test basic simulator for the relative stepper, univariate, model 3: n(u,x)/c(x)

library(devtools)
load_all(".")

set.seed(40)
L <- 2
bbox <- cbind(c(0,L), c(0,L))
R <- L/10
K <- 1
r <- (1:K)*R/K
theta <- c(1,-2)[1:K]
iter <- 1e6


# Free n


theta0 <- log(100/prod(apply(bbox,2,diff)))

xo <- rrelativestepper3(theta=c(theta0, theta), r = r, bbox = bbox, dbg=11, iter = iter)


plot(xo, asp=1, cex=.5)


# crap. effectively theta --> theta/n(x) so larger n less interaction.

