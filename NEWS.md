# melt 1.6.0
### BREAKING CHANGES
* `el_mean()` takes arguments in a different order to comply with the 'tidyverse' style. It takes the data argument `x` first, followed by the parameter specification `par` as `el_mean(x, par)`.

* `lht()` is renamed to `elt()`. 

* `model` argument in `el_mean()`, `el_lm()`, and `el_glm()` are removed. Use `keep_data` in `el_control()`.

### NEW FEATURES
* New package dependencies are added (BH, dqrng, and graphics).

* New `elt()` replaces `lht()`. It accepts additional arguments `alpha` and `calibrate`.
 
* New `el_sd()` performs empirical likelihood test for the standard deviation.

* New `elmt()` tests multiple hypotheses with empirical likelihood.

* New `weights()` extracts the re-scaled weights from a model.

* New `formula()` extracts the model formula used from a model.

* New `nobs()` extracts the number of observations from a model.

* New `conv()` extracts the convergence status from a model.

* New `logLR()` extracts the log empirical likelihood ratio from a model.

* `el_control()` gains additional arguments `verbose`, `keep_data`, `seed`, `b`, and `m`.

### MINOR IMPROVEMENTS
* `cv` argument in `confint()` and `confreg()` defaults to `NULL`. If non-`NULL`, `level` is ignored.

* `probit` link produces more accurate result in `el_glm()`

* `print()` method for an `EL` object shows whether the data are weighted or not.

* All row or column names (if any) of input data are preserved in a fitted `EL` object.

### BUG FIXES
* `confint()` and `confreg()` check if the `cv` argument is compatible with the `th` value set by `control_el()`.

### DEPRECATED AND DEFUNCT
* `el_pairwise()` and `lht()` are removed along with the dependency on the RcppProgress package.


# melt 1.5.2
### NEW FEATURES
* `lht()` accepts both numeric vector and matrix for `lhs` and `rhs` arguments.

* OpenMP parallelization is available for `confint()` by specifying `nthreads` through `control` argument.

### DEPRECATED AND DEFUNCT
* `el_test()` is removed.

* `el_pairwise()` is deprecated and will be removed in a future release. 


# melt 1.5.1
### BUG FIXES
* Unit test errors are fixed.


# melt 1.5.0
### NEW FEATURES
* S4 classes, generics, and methods are adopted throughout the package.

* New `confreg()` constructs confidence regions.

* New `eld()` computes empirical likelihood displacement values.

* New `el_control()` the specifies `control` argument. 

* New `el_glm()` performs empirical likelihood tests to generalized linear models. More families and link functions will be supported in a future release.

* `confint()` gains `cv` argument for a user-supplied critical value.

### DEPRECATED AND DEFUNCT
* `el_aov()` is removed. 

* `el_test()` is deprecated and will be removed in a future release. 


# melt 1.4.0
### NEW FEATURES
* New `lht()` performs linear hypothesis testing.

* New `confint()` constructs confidence intervals.

* New `logLik()` extracts empirical log-likelihood.

### DEPRECATED AND DEFUNCT
* `el_aov()` is deprecated in favor of `el_lm()`. It will be removed in a future release. 


# melt 1.3.0
### NEW FEATURES
* `el_eval()` is added for direct computation with custom estimating functions.

* `el_mean()` and `el_lm()` accepts an optional `weights` argument for weighted EL. Arguments on optimization are now handled by a new `control` argument. It will be used in other functions in future releases.


# melt 1.2.0
### NEW FEATURES
* New `el_lm()` performs empirical likelihood tests for linear models.


# melt 1.1.0
### NEW FEATURES
* New `el_aov()` performs one-way analysis of variance. 


# melt 1.0.1
### BUG FIXES
* Header file issues related to OpenMP and C++ array class are fixed.


# melt 1.0.0
* Released on CRAN.
