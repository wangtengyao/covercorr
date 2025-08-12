# covercorr

This repository contains implementations of the **coverage correlation coefficient**, a statistical measure for comparing two multivariate samples, based on geometric probability and optimal transport ideas.

The repository supports both **Python** and **R**, with structurally separated packages in the `python/` and `R/` subdirectories.

---

## What is Coverage Correlation?

The **coverage correlation coefficient** is a statistic that measures dependence between two random variables/vectors by quantifying the extent to which they have a joint distribution concentrated on a singular subset relative to the product of the marginals. It converges to 0 if and only if the variables are independent and 1 if and only if the copula is singular. It is especially useful in picking up dependencies where both $X$ and $Y$ can be described approximately as functions of a latent variable $U$. It is distribution-free, admits an analytically tractable asymptotic null distribution, and can be computed efficiently, making it well-suited for detecting complex, potentially nonlinear associations in large-scale pairwise testing.

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

---

## Repository Structure

```bash
covercorr/
│
├── R/                           # R package and simulation
│   ├── covercorr/               # R package source code
│   └── simulation_code/         # Simulation scripts and related files
│
├── python/                      # Python package
│   ├── covercorr/               # Python package source code
│   ├── setup.py                 # Packaging configuration
│   └── pyproject.toml           # Build system configuration
│
└── README.md
```

---

## Reference

Yang, X., Azadkia, M. and Wang, T. (2025) Coverage correlation: detecting singular dependencies between random variables. _Preprint_. arxiv:2508.06402.
