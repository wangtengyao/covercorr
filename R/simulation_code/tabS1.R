library(putils)
library(USP)
library(XICOR)
library(Rfast)
library(dHSIC)
library(covercorr)
library(KPC)
source('perm.R')
source('data_generation.R')


# Define a list of test methods as named functions, with skip conditions
method_list <- list(
  covercorr = list(
    run = function(x, y) coverage_correlation(x, y)[c("stat", "pval")],
    skip = function(n, dx) FALSE
  ),
  xicor = list(
    run = function(x, y) xicor(x, y, pvalue = TRUE)[c("xi", "pval")],
    skip = function(n, dx) dx == 2
  ),
  dcor = list(
    run = function(x, y) dcor.test(x, y, B = 100)[c("stat", "pval")],
    skip = function(n, dx) FALSE
  ),
  hsic = list(
    run = function(x, y) dhsic.test(x, y, B = 100)[c("statistic", "p.value")],
    skip = function(n, dx) n == 8000
  ),
  kmac = list(
    run = function(x, y) kmac.test(y, x, B = 100),
    skip = function(n, dx) n >= 4000
  ),
  USP = list(
    run = function(x, y) USPFourier(x, y, M = 3, B = 100),
    skip = function(n, dx) n == 8000
  )
)

results_df <- sim.params(dx=c(1,2),n=c(125, 250, 500, 1000, 2000, 4000, 8000),repeats=1:10)
for (method_name in names(method_list)){
  results_df[[method_name]] <- NA
}

for (i in 1:nrow(settings)){
  bunch(dx, n, r) %=% results_df[i, 1:3]
  x <- matrix(rnorm(n * dx), n)
  y <- matrix(rnorm(n * dx), n)
  for (method_name in names(method_list)){
    method <- method_list[[method_name]]
    if (!method$skip(n, dx)) {
      results_df[i, method_name] <- system.time(method$run(x, y))["user.self"]
    }
  }
  printPercentage(i, nrow(settings))
}

head(results_df)
aggregate(cbind(covercorr, dcor, hsic, kmac, USP) ~ n + dx, FUN=mean, data=results_df)
