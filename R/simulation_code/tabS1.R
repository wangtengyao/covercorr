library(putils)
library(USP)
library(XICOR)
library(Rfast)
library(dHSIC)
library(covercorr)
library(KPC)
source('perm.R')
source('data_generation.R')

n <- 1000
dx <- 1
x <- matrix(rnorm(n*2),n)
y <- matrix(rnorm(n*2),n)
system.time(bunch(stat_kappa, p_kappa) %=% coverage_correlation(x, y)[c('stat','pval')]) # cover corr
system.time(bunch(stat_xi, p_xi) %=% xicor(x, y, pvalue=TRUE)[c('xi','pval')]) # Chatterjee's
system.time(bunch(stat_dcor, p_dcor) %=% dcor.test(x, y, B=100)[c('stat', 'pval')]) # dCor
system.time(bunch(stat_hsic, p_hsic) %=% dhsic.test(x,y,B=100)[c('statistic', 'p.value')]) # HSIC
system.time(bunch(stat_kmac, p_kmac) %=% kmac.test(y, x, B=100))
system.time(bunch(p_USP, stat_USP) %=% USPFourier(x, y, M=3, B=100)) # USP
