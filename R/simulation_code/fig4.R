library(putils)
library(minerva)
data("Spellman")
X <- as.matrix(Spellman[, -1])
bunch(n, p) %=% dim(X)

# compute pairwise correlations
stat <- kappa <- r <- rho <- xi <- matrix(0, p, p)
total <- p * (p-1) / 2
counter <- 0
for (i in 1:(p-1)){
  for (j in (i+1):p){
    x <- X[, i]
    y <- X[, j]
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
kappa <- pmin(kappa * (4381*4380/2), 1)
r <- pmin(r * (4381*4380/2), 1)
rho <- pmin(rho * (4381*4380/2), 1)
xi <- pmin(xi * (4381*4380/2), 1)

# generate plot
selected <- which(kappa <= 0.05 & r > 0.05 & rho > 0.05 & xi > 0.05, 
                  arr.ind = TRUE)
df <- data.frame(i=selected[,1], j=selected[,2])
df$stat <- stat[selected]
df$pval <- kappa[selected]
df <- sort_data_frame(df, 'pval')

for (k in 1:4){
  i <- df$i[k]; j <- df$j[k]
  pdf(paste0('fig/fig4', letters[k], '.pdf'), width=4, height=4)
  par(tcl=-0.3, mar=c(3.5,3.5,2.5,1), mgp=c(2.1,0.6,0), cex.lab=1.5, cex.axis=1.3, family='serif')
  plot(X[,c(i,j)], xlab=colnames(X)[i], ylab=colnames(X)[j])
  dev.off()
}

