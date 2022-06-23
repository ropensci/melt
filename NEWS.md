# melt 1.5.2.9000 (development version)
### BREAKING CHANGES
* The order of arguments of `el_mean()` changed to comply with the 'tidyverse' style. It takes the data argument `x` first, followed by the parameter specification `par` as `el_mean(x, par)`.

* `lht()` is renamed to `elt()`.

### NEW FEATURES
* New package dependencies added (`BH`, `dqrng`, and `graphics`).

* New `elmt()` performs multiple testing with empirical likelihood.

* New `weights()` extracts rescaled weights.

* New `elt()` replaces `lht()`. It accepts additional arguments `alpha` and `calibrate`.

### MINOR IMPROVEMENTS
* `cv` argument in `confint()` and `confreg()` defaults to `NULL`. If non-`NULL`, `level` is ignored.

### BUG FIXES
* `confint()` and `confreg()` check if the `cv` argument is compatible with the `th` value set by `control_el()`.

* unit test errors fixed.

### DEPRECATED AND DEFUNCT
* `el_pairwise()` is removed.


# melt 1.5.2
### NEW FEATURES
* `lht()` accepts both numeric vector and matrix for `lhs` and `rhs` arguments.

* OpenMP parallelization is available for `confint()` by specifying `nthreads` through `control` argument.

### DEPRECATED AND DEFUNCT
* `el_test()` is removed.

* `el_pairwise()` is deprecated and will be removed in a future release. 


# melt 1.5.1
### BUG FIXES
* unit test error fixed.


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

* `el_test` objects additionally return `npar`, `log.prob`, and `loglik`.

### DEPRECATED AND DEFUNCT
* `el_aov()` is deprecated in favor of `el_lm()`. It will be removed in a future release. 


# melt 1.3.0
### NEW FEATURES
* `el_eval()` is added for direct computation with custom estimating functions.

* `melt` class replaced by `el_test` class.

* `el_mean()` and `el_lm()` accepts an optional `weights` argument for weighted EL. Arguments on optimization are now handled by a new `control` argument. It will be used in other functions in future releases.


# melt 1.2.0
### NEW FEATURES
* New `el_lm()` performs empirical likelihood tests for linear models.


# melt 1.1.0
### NEW FEATURES
* New `el_aov()` performs one-way analysis of variance. 


# melt 1.0.1
### BUG FIXES
* Fixed header file issues related to OpenMP and C++ array class.


# melt 1.0.0
### NEW FEATURES
* Currently, 4 functions are available: `el_mean()`, `el_test()`, `el_pairwise()`, and `el_mht()`.

* `clothianidin` data set is added.
