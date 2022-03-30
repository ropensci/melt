
<!-- README.md is generated from README.Rmd. Please edit that file -->

# melt - Multiple Empirical Likelihood Tests

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/melt)](https://CRAN.R-project.org/package=melt)
[![R-CMD-check](https://github.com/markean/melt/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/markean/melt/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/markean/melt/branch/master/graph/badge.svg)](https://app.codecov.io/gh/markean/melt?branch=master)
<!-- badges: end -->

## Overview

The R package **melt** provides a unified framework for data analysis
with empirical likelihood methods. A collection of functions are
available for basic regression analysis and hypothesis testing. Much of
its functionality and syntax mimics the corresponding base R functions.
The core computational routines are implemented with the ‘Eigen’ C++
library and ‘RcppEigen’ interface, with OpenMP for parallel computation.
Additional functions are available for multiple testing for the analysis
of experimental designs. Details of the testing procedures are given in
[Kim et al. (2021)](https://arxiv.org/abs/2112.09206).

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
#> 1 -0.3495232 0.06062649


# linear regression
fit2 <- el_lm(formula = mpg ~ wt, data = mtcars)
summary(fit2)
#> 
#> Call:
#> el_lm(formula = mpg ~ wt, data = mtcars)
#> 
#> Coefficients:
#>             estimate chisq-value p-value    
#> (Intercept)   37.285       443.3  <2e-16 ***
#> wt            -5.344       439.1  <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Chisq: 439.1, df: 1, p-value: < 2.2e-16
confint(fit2)
#>                 lower     upper
#> (Intercept) 33.175865 41.865426
#> wt          -6.746896 -4.149511


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
#> trtNaked       -4.479     411.072  < 2e-16 ***
#> trtFungicide   -3.427      59.486 1.23e-14 ***
#> trtLow         -2.800      62.955 2.11e-15 ***
#> trtHigh        -1.307       4.653    0.031 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Chisq: 1628, df: 4, p-value: < 2.2e-16
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
