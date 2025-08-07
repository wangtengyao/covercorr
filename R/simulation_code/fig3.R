# ##### analyse #####
library(putils)
library(covercorr)
library(MASS)
data("CD8T")
dat <- CD8T

# compute pairwise correlations
bunch(n, p) %=% dim(dat)
stat <- kappa <- xi <- rho <- r <- matrix(0, p, p)
genes <- colnames(dat)
counter <- 0
total <- p * (p-1) / 2
for (i in 1:(p-1)){
  for (j in (i+1):p){
    x <- dat[, i]
    y <- dat[, j]
    r[i, j] <- cor.test(x, y, method='pearson')$p.value
    rho[i, j] <- cor.test(x, y, method='spearman', exact=FALSE)$p.value
    xi[i, j] <- XICOR::xicor(x, y, pvalue=TRUE)$pval
    ret <- covercorr::coverage_correlation(x, y)
    kappa[i, j] <- ret$pval
    stat[i, j] <- ret$stat
    counter <- counter + 1
    if (counter %% (total %/% 1000) == 0){
      println(counter %/% (total %/% 1000), '/1000 completed.')
    }
  }
}

# Bonferonni
r <- pmin(r * (p*(p-1)/2), 1)
rho <- pmin(rho * (p*(p-1)/2), 1)
xi <- pmin(xi * (p*(p-1)/2), 1)
kappa <- pmin(kappa * (p*(p-1)/2), 1)

# generate plot 
selected <- which(r > 0.05 & rho > 0.05 & xi > 0.05 & kappa <= 0.05, 
                  arr.ind=TRUE)
df <- data.frame(i= selected[,1], j=selected[,2], pval=kappa[selected])
df <- sort_data_frame(df, 'pval')


wrap_points <- function(df) {
  shifts <- expand.grid(dx = -1:1, dy = -1:1)
  do.call(rbind, lapply(1:nrow(shifts), function(i) {
    shift <- shifts[i, ]
    data.frame(x = df$x + shift$dx, y = df$y + shift$dy)
  }))
}

visualise_density <- function(x_rank, y_rank){
  points <- data.frame(x=x_rank, y=y_rank)
  wrapped_points <- wrap_points(points)
  dens_wrap <- kde2d(wrapped_points$x, wrapped_points$y, n = 100, lims = c(0, 1, 0, 1))
  normaliser <- mean(dens_wrap$z[dens_wrap$x>=0 & dens_wrap$y >=0 & dens_wrap$x <= 1 & dens_wrap$y <= 1])
  log_ratio_wrap <- log(dens_wrap$z / normaliser)  # since uniform density is 1
  col_palette <- hcl.colors(100, "RdBu", rev = TRUE)
  dens_range <- range(log_ratio_wrap)
  zlim = c(min(dens_range,-1/2), max(dens_range,1/2))
  breaks <- seq(zlim[1], zlim[2], length.out = 101)  # 100 intervals
  fields::image.plot(dens_wrap$x, dens_wrap$y, log_ratio_wrap,
        col = col_palette, breaks = breaks, zlim = zlim,
        xlab = "X rank", ylab = "Y rank", asp=1, 
        mar=c(3.5,3.5,2.5,1), mgp=c(2.1,0.6,0),
        tcl=-0.3, oma=c(0,0,0,0))
}



for (k in 1:2){
  i <- df$i[k]; j <-df$j[k]
  x <- dat[, i]
  y <- dat[, j]
  x_rank <- MK_rank(as.matrix(x), runif(n))
  y_rank <- MK_rank(as.matrix(y), runif(n))
  pdf(paste0('fig/fig3', letters[k*2-1], '.pdf'), width=4.5, height=4.5)
  par(tcl=-0.3, mar=c(3.5,3.5,2.5,1), mgp=c(2.1,0.6,0), cex.lab=1.5, cex.axis=1.3, family='serif')
  plot(x,y, pch='.', xlab=genes[i], ylab=genes[j])
  dev.off()
  pdf(paste0('fig/fig3', letters[k*2], '.pdf'), width=4.75, height=4.5)
  visualise_density(x_rank, y_rank)
  dev.off()
}

