melt 1.6.0
=========================
### BREAKING CHANGES
* The order of arguments of `el_mean()` changed to comply with the 'tidyverse' style. It takes the data argument `x` first, followed by the parameter specification `par`.

### NEW FEATURES
* New package dependencies added (`graphics`, `lifecycle`, `BH`, and `dqrng`).

* Generic function `weights()` is added.

* `lht()` is renamed to `elt()`, which accepts additional arguments `alpha` and `calibrate`.

### MINOR IMPROVEMENTS
* `cv` argument in `confint()` and `confreg()` defaults to NULL.

### BUG FIXES
* unit test error fixed.

### DEPRECATED AND DEFUNCT
* `el_pairwise()` is removed.


melt 1.5.2
=========================
### NEW FEATURES
* `lht()` accepts both numeric vector and matrix for `lhs` and `rhs` arguments.

* OpenMP parallelization is available for `confint()` by specifying `nthreads` through `control` argument.

### DEPRECATED AND DEFUNCT
* `el_test()` is removed.

* `el_pairwise()` is deprecated and will be removed in a future release. 


melt 1.5.1
=========================
### BUG FIXES
* unit test error fixed.


melt 1.5.0
=========================
### NEW FEATURES
* The package depends on more recent version of R (>= 4.0.0).

* S4 classes, generics, and methods are adopted throughout the package.

* Generic functions `confreg()` and `eld()` are added.

* `el_control()` function is added for specifying `control` argument. Plain list it no longer accepted.

* `el_glm()` function is added for generalized linear models. More families and link functions will be supported in a future release.

* `confint()` has additional `cv` argument for a user-supplied critical value.

### DEPRECATED AND DEFUNCT
* `el_aov()` is removed. 

* `el_test()` is deprecated and will be removed in a future release. 


melt 1.4.0
=========================
### NEW FEATURES
* C++14 standards are used for the package.

* `lht()` is added for linear hypothesis testing.

* Generic functions `confint()` and `logLik()` are added for `el_test` class.

* `el_test` objects additionally return `npar`, `log.prob`, and `loglik`.

### DEPRECATED AND DEFUNCT
* `el_aov()` deprecated in favor of `el_lm()`. It will be removed in a future release. 


melt 1.3.0
=========================
### NEW FEATURES
* `el_eval()` function added for direct computation with custom estimating functions.

* `melt` class replaced by `el_test` class.

* `el_mean()` and `el_lm()` accepts an optional `weights` argument for weighted EL. Arguments on optimization are now handled by a new `control` argument. It will be used in other functions in future releases.


melt 1.2.0
=========================
### NEW FEATURES
* Dependence on R version updated to '3.6.0'.

* `el_lm()` function added for linear regression analysis.


melt 1.1.0
=========================
### NEW FEATURES
* `el_aov()` function added for one-way analysis of variance. It only supports one variable at the moment.


melt 1.0.1
=========================
### BUG FIXES
* Fixed header file issues related to OpenMP and C++ array class.


melt 1.0.0
=========================
### NEW FEATURES
* Currently, 4 functions are available: `el_mean()`, `el_test()`, `el_pairwise()`, and `el_mht()`.

* `clothianidin` data set added.
