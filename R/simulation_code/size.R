library(putils)
library(covercorr)

params <- sim.params(n=c(10, 100, 1000), dx=c(1,2), repeats=1:5000, shuffle=TRUE)
for (i in 1:nrow(params)){
  bunch(n, dx, repeats) %=% params[i, ]
  x <- matrix(runif(n*dx), nrow=n)
  y <- matrix(runif(n*dx), nrow=n)
  d <- dx * 2
  
  eps <- n^(-1/d) / 2
  zmin <- cbind(x - eps, y - eps)
  zmax <- cbind(x + eps, y + eps)
  
  ret <- split_rectangles(zmin, zmax)
  zmin_s <- ret$zmin
  zmax_s <- ret$zmax
  
  total_volume <- covered_volume_partitioned(zmin_s, zmax_s)
  excess_vacancy <- 1 - exp(-1) - total_volume
  kappa <- excess_vacancy / (1 - exp(-1))
  sd <- sqrt(variance_formula(n, d)) 
  Z <- excess_vacancy * sqrt(n) / sd 
  pval <- 1 - pnorm(Z)

  println(show.params(n, dx, pval))
}

