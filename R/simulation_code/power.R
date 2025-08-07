library(putils) # remotes::install_github('wangtengyao/putils')
library(USP)
library(XICOR)
library(Rfast)
library(dHSIC)
library(covercorr) # remotes::install_github('wangtengyao/covercorr')
library(KPC)
source('perm.R')
source('data_generation.R')

params <- sim.params(type=c('sinusoidal', 'zigzag', 'circle', 'spiral', 'Lissajous', 'local'), 
                     n=c(1000,2000), noise=seq(0, 2, length=11), repeats=1:10, shuffle=TRUE)

for (i in 1:nrow(params)){
  bunch(type, n, noise, repeats) %=% params[i, ]
  dx <- ifelse(n == 1000, 1, 2)
  bunch(x, y) %=% data_generation(n, type, noise, dx)
  bunch(stat_kappa, p_kappa) %=% coverage_correlation(x, y)[c('stat','pval')] # cover corr
  bunch(p_USP, stat_USP) %=% USPFourier(x, y, M=3, B=100) # USP
  bunch(stat_dcor, p_dcor) %=% dcor.test(x, y, B=100)[c('stat', 'pval')] # dCor
  if (dx == 1){
    bunch(stat_xi, p_xi) %=% xicor(x, y, pvalue=TRUE)[c('xi','pval')] # Chatterjee's
  } else {
    bunch(stat_xi, p_xi) %=% c(NA, NA)
  }
  bunch(stat_hsic, p_hsic) %=% dhsic.test(x,y,B=100)[c('statistic', 'p.value')] # HSIC
  bunch(stat_kmac, p_kmac) %=% kmac.test(y, x, B=100)
  println(show.params(type, n, dx, noise, stat_kappa, p_kappa, stat_xi, p_xi, stat_dcor, p_dcor, 
                      stat_hsic, p_hsic, stat_kmac, p_kmac, stat_USP, p_USP))
}
