
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
install.packages("melt")
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
# one sample test for mean
fit1 <- el_mean(par = 0, x = rnorm(n = 100))
confint(fit1)
#>        lower      upper
#> 1 -0.3192432 0.09835835


# linear regression
fit2 <- el_lm(formula = mpg ~ wt, data = mtcars)
summary(fit2)
#> 
#> Call:
#> el_lm(formula = mpg ~ wt, data = mtcars)
#> 
#> Coefficients:
#>             estimate chisq-value  p-value    
#> (Intercept)   37.285       44.11 3.10e-11 ***
#> wt            -5.344       41.67 1.08e-10 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Chisq: 41.67, df: 1, p-value: 1.08e-10
confint(fit2)
#>                lower     upper
#> (Intercept) 34.49139 41.180745
#> wt          -6.21663 -4.413012


# analysis of variance 
data("clothianidin")
fit3 <- el_lm(clo ~ -1 + trt, clothianidin)
summary(fit3)
#> 
#> Call:
#> el_lm(formula = clo ~ -1 + trt, data = clothianidin)
#> 
#> Coefficients:
#>              estimate chisq-value  p-value    
#> trtNaked       -4.479      60.781 6.38e-15 ***
#> trtFungicide   -3.427      56.844 4.72e-14 ***
#> trtLow         -2.800      41.587 1.13e-10 ***
#> trtHigh        -1.307       4.653    0.031 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Chisq: 163.9, df: 4, p-value: < 2.2e-16
confint(fit3)
#>                  lower      upper
#> trtNaked     -5.002118 -3.9198229
#> trtFungicide -4.109816 -2.6069870
#> trtLow       -3.681837 -1.9031795
#> trtHigh      -2.499165 -0.1157222


# test of no treatment effect
lhs <- matrix(c(1, -1, 0, 0,
                0, 1, -1, 0,
                0, 0, 1, -1), byrow = TRUE, nrow = 3)
lht(fit3, lhs)
#> 
#> Empirical Likelihood Linear Hypothesis Test: lm 
#> 
#> Chisq: 26.6, df: 3, p-value: 7.148e-06
#> maximum EL estimates:
#> [1] -3.368 -3.368 -3.368 -3.368
```
