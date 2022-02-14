
<!-- README.md is generated from README.Rmd. Please edit that file -->

# melt - Multiple Empirical Likelihood Tests

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/melt)](https://CRAN.R-project.org/package=melt)
[![R-CMD-check](https://github.com/markean/melt/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/markean/melt/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

The R package **melt** provides a unified framework for data analysis
with empirical likelihood methods. A collection of functions are
available for basic regression analysis and hypothesis testing. Much of
its functionality and syntax are designed to mimic the corresponding
base R functions. The core routines are written in C++ and utilize
OpenMP for parallelization.

## Installation

You can install the latest stable release from
[CRAN](https://cran.r-project.org/package=melt).

``` r
install.packages("melt", dependencies = TRUE)
```

You can install the latest development version from
[Github](https://github.com/markean/melt).

``` r
# install.packages("devtools")
devtools::install_github("markean/melt")
```

## Usage

``` r
library(melt)
# one sample test (mean)
el_mean(par = 0, x = rnorm(n = 100))

# linear regression
fit <- el_lm(formula = mpg ~ wt, data = mtcars)
summary(fit)

# analysis of variance
el_aov(formula = Sepal.Length ~ Species, data = iris)
```
