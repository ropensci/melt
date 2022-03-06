# melt 1.3.0
* `el_eval` function added for direct computation with custom estimating functions.

* `melt` class replaced by `el_test` class.

* `el_mean` and `el_lm` accepts an optional `weights` argument for weighted EL. Arguments on optimization are now handled by a new `control` argument. It will be used in other functions in future releases.


# melt 1.2.0
* `el_lm` function added for linear regression analysis.

* Dependence on R version updated to '3.6.0'.

* Reference style streamlined across functions and data.


# melt 1.1.0

* DESCRIPTION file modified. 

* `el_aov` function added for one-way analysis of variance. It only supports one variable at the moment.


# melt 1.0.1

* Fixed header file issues related to OpenMP and C++ array class.

* `el_pairwise` documentation specifies the use of OpenMP for NB procedure.


# melt 1.0.0
* Currently, 4 functions are available: `el_mean`, `el_test`, `el_pairwise`, and `el_mht`.

* `clothianidin` data set added.
