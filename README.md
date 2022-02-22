
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
#> 
#>  One sample EL test
#> 
#> data: rnorm(n = 100)
#> Chisq = 1.594, df = 1, p-value = 0.2068
#> alternative hypothesis: true  is not equal to 0
#> maximum EL estimates:
#> [1] 0.1186762

# linear regression
fit <- el_lm(formula = mpg ~ wt, data = mtcars)
summary(fit)
#> 
#> Call:
#> el_lm(formula = mpg ~ wt, data = mtcars)
#> 
#>     Min      1Q  Median      3Q     Max 
#> -4.5432 -2.3647 -0.1252  1.4096  6.8727 
#> 
#> Coefficients:
#>             estimate chisq-value p-value    
#> (Intercept)   37.285       88.41  <2e-16 ***
#> wt            -5.344       84.82  <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Multiple R-squared:  0.7528, Adjusted R-squared:  0.7446 
#> Chisq-statistic: 84.82 on 1 DF, p-value: < 2.2e-16

# analysis of variance
el_aov(formula = Sepal.Length ~ Species, data = iris)
#> Warning in el_aov(formula = Sepal.Length ~ Species, data = iris): convergence failed
#> Call:
#> el_aov(formula = Sepal.Length ~ Species, data = iris)
#> 
#> minimizer:
#> 5.8433 5.8433 5.8433
#> 
#> statistic:
#> 133.2364
```
