# test multi-simulator
library(devtools)
load_all(".")
library(Kcross)
library(spatstat)

# check basic poisson
if(0){
  theta0 <- log(c(50, 50, 300, 10))
  bbox <- cbind(c(-1,.4), c(0,1.5))
  V <- prod(apply(bbox,2, diff))
  x <- rstepper_multi(theta0, dbg = 11, iter=2e4, bbox = bbox)
  print(rbind(exp=exp(theta0)*V, table(x[,3])))
  plot(x[,-3], col=1+x[,3], pch=1+2*x[,3], asp=1)
}

# check bivariate
if(0){
  set.seed(1)
  # Fiddle with these
  theta0 <- log( c(100, 100))
  theta1 <- list(list(r= 0.04, theta = -5),
                 list(r= 0.04, theta = -5))
  theta2 <- list(list(r= 0.02, theta = 2.5))

  x <- rstepper_multi(theta0, theta1, theta2, iter = 1e5, dbg=1, toroidal=T)
  xx <- x[,-3]
  m <- x[,3]
  print(table(m))
  pp <- ppp(xx[,1],xx[,2], marks=factor(m))
  ppl <- split(pp)
  g12 <- pcf_cross_all_box(pp, r = r<-seq(0, .1, l =50), adjust=1.5)
  #
  par(mfrow=c(2,1))
  plot(xx, col=m+1, pch=1+2*m, asp=1)
  plot(r, g12[1,1,], col=1, ylim=c(0,2), type="l")
  lines(r,g12[2,2,], col=2)
  lines(r,g12[1,2,], col=3)
  abline(h=1)
}
# seems to work.


# check trivariate
if(0){
  set.seed(1)
  # Fiddle with these
  theta0 <- log( c(200, 200, 200))
  theta1 <- list(list(r= 0.02, theta = -1),
                 list(r= 0.02, theta = -1),
                 list(r= 0.02, theta = -1))
  theta2 <- list(list(r= 0.01, theta=1),
                 list(r= 0.01, theta=1),
                 list(r= 0.02, theta=0))
  # simulate
  x <- rstepper_multi(theta0, theta1, theta2, iter = 1e5, dbg=200, toroidal=T)
  xx <- x[,-3]
  m <- x[,3]
  print(table(m))
  pp <- ppp(xx[,1],xx[,2], marks=factor(m))
  ppl <- split(pp)
  g12 <- pcf_cross_all_box(pp, r = r<-seq(0, .1, l =50), adjust=1.5)
  #
  par(mfrow=c(2,1))
  plot(xx, col=m+1, pch=1+2*m, asp=1, cex=.6)
  plot(r, g12[1,1,], col=1, ylim=c(0,2), type="l")
  lines(r,g12[2,2,], col=2)
  lines(r,g12[3,3,], col=3)
  lines(r,g12[1,2,], col=2, lty=2)
  lines(r,g12[1,3,], col=3, lty=2)
  lines(r,g12[2,3,], col=4, lty=2)
  abline(h=1)

}

# Works!
