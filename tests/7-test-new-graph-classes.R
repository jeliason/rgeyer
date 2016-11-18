# testing of the basic functions of the c-classes
library(devtools)
load_all(".")

set.seed(2)
x0 <- matrix(runif(500),nc=2)
bbox <- cbind(0:1, 0:1)
iter <- 10
p <- 9
m <- sample(0:p, rep=T, nrow(x0))

x <- cbind(x0, m)

K <- 5
R <- 0.1
r <- 1:K/K * R

o <- rtest_graph_c(x, bbox, r)


## Plot the pattern
plot(x0, asp=1, col = 1+m, cex=0)
text(x0[,1],x0[,2], col=1+m, 1:nrow(x0))
for(s in r) symbols(x0, circles = rep(s, nrow(x0)), add=T, inches=F, fg=1+m)


## Overall, makes sense only for K=1
pairs <- unlist(sapply(0:(p-1), function(i) paste(i, (i+1):p, sep="-")))
# print("Pairs")
# pv <- do.call(cbind,o[[2]])
# colnames(pv) <- pairs
# print(pv)



# for one point, return neighbour counts, works for p < 10
fori <- function(i, ...) {
  mi <- m[i]
  # i-i
  cat("Point", i , "of type", mi, " neighbour counts:\n")
  n <- NULL
  for(k in 1:K){
    # intra
    nk <- rep(0, p+1)
    nk[mi+1] <- o[[1]][[mi+1]][i,k]
    j <- grep(mi, pairs)
    mp <- pairs[j]
    mp <- as.integer(gsub(paste0("-",mi, "|",mi,"-"), "", mp))
    nk[mp+1] <- sapply(o[[2]][j], function(z) z[i,k])
    n <- rbind(n, nk)
  }
  colnames(n) <- paste0("type", 0:p)
  rownames(n) <- paste0("k",1:K)
  print(n)
  plot(x0, asp=1, col = 1+m, cex=0)
  text(x0[,1],x0[,2], col=1+m, 1:nrow(x0), ...)
  for(s in r) symbols(x0[i,1],x0[i,2], circles = s, add=T, inches=F, fg=1+m[i])
}

