## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(covercorr)

## -----------------------------------------------------------------------------
n <- 100
p <- 1
X <- rnorm(n)
Y <- rnorm(n)
result <- coverage_correlation(Y, X, visualise = TRUE)
str(result)

## -----------------------------------------------------------------------------
n <- 100
p <- 1
X <- rnorm(n)
Z <- rnorm(n)
rho <- 0.9
Y <- rho * X + sqrt(1 - rho^2) * Z
result <- coverage_correlation(Y, X, visualise = TRUE)
str(result)

## -----------------------------------------------------------------------------
n <- 100
p <- 2
X <- matrix(rnorm(p * n), ncol = p)
Y <- matrix(0, nrow = n, ncol = p)
Y[, 1] <- X[, 1]^2
Y[, 2] <- X[, 1] * X[, 2]
result <- coverage_correlation(Y, X)
str(result)

## -----------------------------------------------------------------------------
n <- 50
p <- 2
X <- matrix(rnorm(p * n), ncol = p)
Y <- matrix(rnorm(p * n), ncol = p)
result <- coverage_correlation(Y, X, method = 'approx')
str(result)

