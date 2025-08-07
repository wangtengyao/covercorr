library(covercorr)
library(XICOR)
library(putils)


generate_plot <- function(counter){
  # raw plot
  pdf(paste0('fig/fig1', letters[counter*3-2], '.pdf'), width=4, height=4)
  par(tcl=-0.3, mar=c(3.5,3.5,1,1), mgp=c(2.1,0.6,0), cex.lab=1.5, cex.axis=1.3, family='serif')
  plot(x, y, xlab='X', ylab='Y')
  dev.off()
  
  # chatterjee's line plot
  pdf(paste0('fig/fig1', letters[counter*3-1], '.pdf'), width=4, height=4)
  par(tcl=-0.3, mar=c(3.5,3.5,1,1), mgp=c(2.1,0.6,0), cex.lab=1.5, cex.axis=1.3, family='serif')
  ord <- order(x)
  r <- rank(y[ord])
  plot((1:n)/n, r/n, type='l', xlab='X ranks', ylab='Y ranks')
  dev.off()
  
  # coverage plot
  pdf(paste0('fig/fig1', letters[counter*3], '.pdf'), width=4, height=4)
  par(tcl=-0.3, mar=c(3.5,3.5,1,1), mgp=c(2.1,0.6,0), cex.lab=1.5, cex.axis=1.3, family='serif')
  u <- runif(n)
  v <- runif(n)
  x_rank <- sort(u)[rank(x)]
  y_rank <- sort(v)[rank(y)]
  
  eps <- n^(-1/2) / 2
  zmin <- cbind(x_rank - eps, y_rank - eps)
  zmax <- cbind(x_rank + eps, y_rank + eps)
  
  bunch(zmin_splitted, zmax_splitted) %=% split_rectangles(zmin, zmax)
  xmin <- zmin_splitted[, 1]
  ymin <- zmin_splitted[, 2]
  xmax <- zmax_splitted[ ,1]
  ymax <- zmax_splitted[, 2]

  plot(x_rank, y_rank, pch=20, cex=0.3, asp=1, frame=FALSE, xlab='X ranks', ylab='Y ranks')
  plot_rectangles(xmin, xmax, ymin, ymax, add=TRUE)
  dev.off()
  
  bunch(stat_xi, p_xi) %=% xicor(x, y, pvalue=TRUE)[c('xi','pval')] # Chatterjee's
  bunch(stat_kappa, p_kappa) %=% coverage_correlation(x, y)[c('stat','pval')] # cover corr
  println(show.params(stat_xi, p_xi))
  println(show.params(stat_kappa, p_kappa))
}

set.seed(42)
# example 1
n <- 1000
x <- rnorm(n); y <- rnorm(n)
generate_plot(1)

set.seed(42)
# example 2
x <- rnorm(n)
y <- sin(x * 10) + rnorm(n) * 0.5
generate_plot(2)

set.seed(42)
# example 3
x <- rnorm(n)
y <- rnorm(n)
idx <- sample(c(TRUE,FALSE), n, replace=TRUE, prob=c(0.5, 0.5))
y[idx] <- x[idx]
generate_plot(3)

set.seed(42)
# example 4
t <- seq(0, 1, length=n)
x <- t * sin(t * 30) + rnorm(n) * 0.01
y <- t * cos(t * 30) + rnorm(n) * 0.01
generate_plot(4)

