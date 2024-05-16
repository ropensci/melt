# melt 1.11.4
## MINOR IMPROVEMENTS
* The internal pseudo-random number generator has been switched from Xoshiro256+ to Xoshiro256++ with the release of the new dqrng package.


# melt 1.11.3
## NEW FEATURES
* New `"ael"` option has been added in the `calibrate` argument of `elt()` for adjusted empirical likelihood calibration.

## MINOR IMPROVEMENTS
* The package vignette has been updated.


# melt 1.11.2
## MINOR IMPROVEMENTS
* The package vignette has been updated.


# melt 1.11.1
## MINOR IMPROVEMENTS
* Some function arguments now utilize the checkmate package for validation.

* The package vignette has been updated.


# melt 1.11.0
## MINOR IMPROVEMENTS
* Updated package vignette with the publication in the Journal of Statistical Software.

## DEPRECATED AND DEFUNCT
* Removed `el_pairwise()` and associated methods.

* Removed `sigTests()` for objects inheriting from `SummaryLM`.


# melt 1.10.0
## NEW FEATURES
* `el_glm()` accepts `quasipoisson` family with `"sqrt"` link function for the argument `family`.

## DEPRECATED AND DEFUNCT
* `sigTests()` is deprecated in favor of `coef()` for an object that inherits from `SummaryLM` and will be removed in a future release. 

* `logLik()` is removed.


# melt 1.9.0
## NEW FEATURES
* `confint()` is applicable to an `EMLT` object to produce simultaneous confidence intervals.

* All model objects gain `control` slot of `ControlEL` class. All methods that apply to these objects inherit `control` unless it is overwritten by the user explicitly.

## MINOR IMPROVEMENTS
* `summary()` is applicable to an object that inherits from `EL`, `ELT`, and `EMLT`.

* A more informative message is printed regarding the convergence status.

* `optim` slot in all model or summary objects gains a single logical element `cstr` that shows whether a constrained EL computation is involved or not. 

## DEPRECATED AND DEFUNCT
* `logLik()` is deprecated and will be removed in a future release. 

## BUG FIXES
* `confreg()` checks whether `parm` matches the parameters in `object` when a `character` vector is specified for `parm`.


# melt 1.8.0
## NEW FEATURES
* New accessor method `logProb()` extracts a model's log probabilities of empirical likelihood.

* `el_lm()` and `el_glm()` gain an argument `offset`.

* `el_glm()` accepts `quasipoisson` family with `"identity"` link function for the argument `family`.

* `elt()` accepts a character vector for the argument `lhs`, allowing a symbolic description of a hypothesis.

* `eltmt()` accepts a character vector as an element of the argument `lhs`, allowing a symbolic description of hypotheses.

* `plot()` applies to an object that inherits from `EL` to plot empirical likelihood displacement values versus observation index.

* New dataset `thiamethoxam` added.

## MINOR IMPROVEMENTS
* `coef()` and `getDF()` is applicable to an object of class `EMLT`.

* `print()` shows the tested hypothesis when applied to an object of class `ELT`.

* `print()` shows the tested hypotheses, the estimates, and marginal degrees of freedom when applied to an object of class `ELMT`. The description of the hypotheses and the estimates are printed only when the marginal degrees of freedom are all one.

* `"boot"` option in the `calibrate` argument of `elt()` yields a more reliable result when applied to an object that inherits from `LM`.

* Internal routines for projection operation do not compute an explicit inverse (thanks to @awstringer1).

## BUG FIXES
* `elmt()` returns a correct critical value when applied to an object of class `QGLM`.

* `"boot"` option in the `calibrate` argument of `elt()` works with an object of class `SD`.


# melt 1.7.0
## NEW FEATURES
* `el_glm()` accepts `quasipoisson` family with `"log"` link function for the argument `family`.

* New accessor methods added (`chisq()`, `critVal()`, `getDF()`, `getOptim()`, `sigTests()`, `logL()`, and `pVal()`).

* `conv()` is applicable to an object returned by `summary()`.

## MINOR IMPROVEMENTS
* `print()` shows class-specific information.

* `p.value` returned by `el_eval()` is renamed to `pval` for consistency with other functions.

## BUG FIXES
* `confint()` and `confreg()` are not applicable to an object whose `data` is `NULL`.


# melt 1.6.0
## BREAKING CHANGES
* `el_mean()` takes arguments in a different order to comply with the 'tidyverse' style. It takes the data argument `x` first, followed by the parameter specification `par` as `el_mean(x, par)`.

* `lht()` is renamed to `elt()`. 

* `model` argument in `el_mean()`, `el_lm()`, and `el_glm()` are removed. Use `keep_data` in `el_control()`.

## NEW FEATURES
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

## MINOR IMPROVEMENTS
* `cv` argument in `confint()` and `confreg()` defaults to `NULL`. If non-`NULL`, `level` is ignored.

* `probit` link produces a more accurate result in `el_glm()`.

* `print()` for an `EL` object shows whether the data are weighted or not.

* All row or column names (if any) of input data are preserved in a fitted `EL` object.

## BUG FIXES
* `confint()` and `confreg()` check if the `cv` argument is compatible with the `th` value set by `control_el()`.


# melt 1.5.2
## NEW FEATURES
* `lht()` accepts both numeric vector and matrix for `lhs` and `rhs` arguments.

* OpenMP parallelization is available for `confint()` by specifying `nthreads` through `control` argument.

## DEPRECATED AND DEFUNCT
* `el_test()` is removed.

* `el_pairwise()` is deprecated and will be removed in a future release. 


# melt 1.5.1
## BUG FIXES
* Unit test errors are fixed.


# melt 1.5.0
## NEW FEATURES
* S4 classes, generics, and methods are adopted throughout the package.

* New `confreg()` constructs confidence regions.

* New `eld()` computes empirical likelihood displacement values.

* New `el_control()` the specifies `control` argument. 

* New `el_glm()` performs empirical likelihood tests to generalized linear models. More families and link functions will be supported in a future release.

* `confint()` gains `cv` argument for a user-supplied critical value.

## DEPRECATED AND DEFUNCT
* `el_aov()` is removed. 

* `el_test()` is deprecated and will be removed in a future release. 


# melt 1.4.0
## NEW FEATURES
* New `lht()` performs linear hypothesis testing.

* New `confint()` constructs confidence intervals.

* New `logLik()` extracts empirical log-likelihood.

## DEPRECATED AND DEFUNCT
* `el_aov()` is deprecated in favor of `el_lm()`. It will be removed in a future release. 


# melt 1.3.0
## NEW FEATURES
* `el_eval()` is added for direct computation with custom estimating functions.

* `el_mean()` and `el_lm()` accepts an optional `weights` argument for weighted EL. Arguments on optimization are now handled by a new `control` argument. It will be used in other functions in future releases.


# melt 1.2.0
## NEW FEATURES
* New `el_lm()` performs empirical likelihood tests for linear models.


# melt 1.1.0
## NEW FEATURES
* New `el_aov()` performs a one-way analysis of variance. 


# melt 1.0.1
## BUG FIXES
* Header file issues related to OpenMP and C++ array class are fixed.


# melt 1.0.0
* Released on CRAN.
