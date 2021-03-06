
# solid: A Divide-and-Conquer Method for Sparse Risk Prediction and Evaluation

[![CRAN](https://www.r-pkg.org/badges/version/solid)](https://CRAN.R-project.org/package=solid)

## Overview

Divide-and-conquer (DAC) is a commonly used strategy to overcome the
challenges of extraordinarily large data, by first breaking the dataset
into series of data blocks, then combining results from individual data
blocks to obtain a final estimation. We propose a screening and one-step
linearization infused DAC (SOLID) algorithm to fit sparse logistic
regression to massive datasets, by integrating the DAC strategy with a
screening step and sequences of linearization.

The package `solid` consists 1) a screening and one-step linearization
infused DAC (SOLID) algorithm to fit sparse logistic regression to
massive datasets, and 2) a modified cross-validation (MCV) that utilizes
the side products of the SOLID hence substantially reduce the
computational burden.

## Installation

Install development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("celehs/solid")
```

# Citation

Hong C, Wang Y, Cai T. A divide-and-conquer method for sparse risk
prediction and evaluation. Biostatistics (Oxford, England). 2020
Sep. <https://doi.org/10.1093/biostatistics/kxaa031>
