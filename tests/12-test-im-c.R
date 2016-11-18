# Test the im-handler in c

library(devtools)
load_all(".")

library(spatstat)
bbox <- cbind(c(0,3), c(0,1.5))
#x <- rpoispp(10, win = W <- as.owin(c(bbox)))
#nc <- 50
#y <- density(x, dimyx=c(nc,nc))
#y <- as.im(function(x,y) y^2+x^2, W = W, dimyx = c(nc,nc))

y <- exp(readRDS("~/temp.rds"))
bbox <- cbind(y$xrange, y$yrange)
W <- as.owin(c(bbox))
V <- area(W)
xs <- y$xcol #seq(0,bbox[2,1], l = nc)
ys <- y$yrow #seq(0,bbox[2,2], l = nc)
nc <- length(xs)
g <- as.matrix( expand.grid(xs,ys) )

v <- rtest_im_c(y, g)
#v0 <- rtest_im_c(list(), g)

set.seed(15)
t0 <- rpoispp(10/V, win = W)

#y0 <- im(matrix(v0, ncol=nc), xs, ys)
y2 <- im(matrix(v, ncol=nc, byrow=T), xs, ys)
y3 <- im(matrix(y[ppp(g[,1],g[,2], window = W)], byrow=T, ncol=nc), xs, ys)

par(mfrow=c(2,2))

m0 <- y0[t0]
m1 <- y[t0]
m2 <- y2[t0]
m3 <- y3[t0]
print(rbind(m1,m2,m3))
z <- 1
plot(y); points(t0, cex=m1/z)
plot(y2); points(t0, cex=m2/z)
plot(y3); points(t0, cex=m3/z)
plot(y0); points(t0, cex=m0/z)
