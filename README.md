# covercorr

This repository contains implementations of the **coverage correlation coefficient**, a statistical measure for comparing two multivariate samples, based on geometric probability and optimal transport ideas.

The repository supports both **Python** and **R**, with structurally separated packages in the `python/` and `R/` subdirectories.

---

## 🔍 What is Coverage Correlation?

The **coverage correlation coefficient** is a statistical measure that quantifies dependence between two random vectors by computing the union volume of data-centered hyperrectangles in a uniform space. It is especially useful in multivariate, nonlinear, or non-Euclidean contexts.

For a sample of paired observations $(X_i, Y_i)$, the statistic:

- Transforms each sample to multivariate ranks using optimal transport.
- Places a small cube around each ranked point and compute the union volume
- Standardizes this to obtain a correlation-like measure and a p-value.

---

## Python Package

All Python code is under the `python/` subdirectory. To install directly from GitHub:

```bash
pip install git+https://github.com/wangtengyao/covercorr.git#subdirectory=python
```

### Usage:

```python
from covercorr import coverage_correlation
import numpy as np

x = np.random.rand(100, 2)
y = np.random.rand(100, 2)

kappa, pval = coverage_correlation(x, y)
print(f"Coverage correlation: {kappa}, p-value: {pval}")
```

---

## R package

All R code is under the `R/` subdirectory. To install from the GitHub:

```r
remotes::install_github('wangtengyao/covercorr', subdir='R/covercorr')
```

### Usage

```r
library(covercorr)

x <- matrix(runif(200), ncol = 2)
y <- matrix(runif(200), ncol = 2)

result <- coverage_correlation(x, y)
print(result)
```

## Repository Structure

```bash
covercorr/
│
├── R_package/     # R package
│
├── python/        # Python package
│
└── README.md
```
