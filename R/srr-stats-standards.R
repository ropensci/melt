#' srr_stats
#'
#' All of the following standards initially have `@srrstatsTODO` tags.
#' These may be moved at any time to any other locations in your code.
#' Once addressed, please modify the tag from `@srrstatsTODO` to `@srrstats`,
#' or `@srrstatsNA`, ensuring that references to every one of the following
#' standards remain somewhere within your code.
#' (These comments may be deleted at any time.)
#'
#' @srrstatsVerbose FALSE
#'
#' @srrstats {G1.2} Life cycle statement is included in
#'   `.github/CONTRIBUTING.md`.
#' @srrstats {G2.15} All user-facing functions in the package inspect missing
#'   values in inputs. Some functions produce an error and some functions allow
#'   users to apply appropriate treatments to the missing values. These inputs
#'   are not passed to base routines like `mean()`, `sd()`, `cor()`, etc. See
#'   `R/validate.R` as well.
#' @srrstats {G5.2, G5.2a, G5.2b} Every message produced within the package is
#'   unique. Unit tests in `test/testthat` cover all the conditions that trigger
#'   any message. See `R/validate.R` as well.
#' @srrstats {G5.12} Some unit tests compare parallel computation using OpenMP
#'   with single threaded computation. Those tests are run regardless of the
#'   availability of OpenMP without affect the test results. See
#'   `.github/CONTRIBUTING.md`.
#' @srrstats {RE5.0} Documented in the package [website](https://markean.github.io/melt/articles/performance.html).
#' @noRd
NULL


#' NA_standards
#'
#' Any non-applicable standards can have their tags changed from `@srrstatsTODO`
#' to `@srrstatsNA`, and placed together in this block, along with explanations
#' for why each of these standards have been deemed not applicable.
#' (These comments may also be deleted at any time.)
#' @srrstatsNA {G1.5} The package does not have any associated publication.
#' @srrstatsNA {G1.6} The package does not make any performance claim,
#'   especially compared to other R packages.
#' @srrstatsNA {G2.4d, G2.4e} No function in the package expects a factor
#'   variable as an input and hence conversion to or from a factor is not
#'   required.
#' @srrstatsNA {G2.5} The package does not use any `factor` type input.
#' @srrstatsNA {G2.10} No function in the package extracts or filters a column
#'   from tabular inputs.
#' @srrstatsNA {G2.14c} Any form of imputation not only affects the fitting
#'   process but also the interpretation of outcomes. Imputation is beyond the
#'   scope of the package. We believe that users should implement imputation
#'   themselves in a way that empirical likelihood methods and the package are
#'   applicable.
#' @srrstatsNA {G3.1, G3.1a, RE4.6} Empirical likelihood methods (and thus
#'   the package) in general do not require covariance calculations. The
#'   "implicit studentization" is one of the merits of empirical likelihood.
#' @srrstatsNA {G4.0} The package does not allow outputs to be written to local
#'   files.
#' @srrstatsNA {G5.4b, G5.4c} The package concerns G5.4 and G5.4a.
#' @srrstatsNA {G5.10} Unit tests take less then 30 seconds in all platforms we
#'   tested. We avoid long running tests in compliance with the CRAN guidelines.
#'   We believe that the separate flags are not necessary for the package.
#' @srrstatsNA {G5.11} Unit tests do not require large data sets.
#' @srrstatsNA {G5.11a} Unit tests do not download additional data sets.
#' @srrstatsNA {RE2.3} It is beyond the scope of estimating functions that the
#'   package considers. There seems to be no simple way to pass this information
#'   to the internal routines. Users can use `el_eval()` to manually construct
#'   estimating functions with the offset included, but any further method not
#'   applicable to the output of `el_eval()`.
#' @srrstatsNA {RE4.1} The focus of the package is not on fitting but on
#'   inference through hypothesis testing. Internal routines rely on the
#'   coefficient estimates from functions such as `lm.fit()` and `glm.fit()` as
#'   a starting point to produce values need to construct various S4 class
#'   objects. It is certainly possible to generate a model object without
#'   fitting values, but those objects cannot be applied to any of the methods
#'   that the package provides. Also, `lm.fit()` and `glm.fit()` are fast enough
#'   for most use cases. We believe that the feature to generate a model without
#'   coefficients are not useful for the package.
#' @srrstatsNA {RE4.10} The package does not compute the fitted values. Hence,
#'   the residuals are unavailable.
#' @srrstatsNA {RE4.12} No function in the package transforms input data, except
#'   for the use of `model.response()` or `model.matrix()` in `el_lm()` and
#'   `el_glm()`. The same operations are performed by `lm()` and `glm()`, so
#'   including such methods in the package seems redundant.
#' @srrstatsNA {RE4.14, RE4.15, RE4.16, RE7.4} Empirical likelihood methods (and
#'   thus the package) are inherently not suitable for prediction,
#'   extrapolation, or forecasting. Especially, extrapolation can often violate
#'   the convex hull constraint, leading to an invalid empirical likelihood
#'   value. For this reason, the package does provide any `predict()` method.
#' @srrstatsNA {RE6.1} All `plot()` methods in the package are S4 generics.
#' @srrstatsNA {RE6.2} No function in the package produces `fitted` values.
#' @srrstatsNA {RE4.8, RE4.9, RE4.13} The formula interface used by `el_lm()`
#'   and `el_glm()` are just a `lm()` like wrapper to transform the input `data`
#'   to an appropriate estimating equations and pass them to the internal
#'   routines. The distinction between the response and the predictors is not
#'   important in terms of estimating functions and empirical likelihood. In
#'   addition, the package has no `fitted()` or `resid()` methods where the
#'   extraction can be useful. This is intentional since the empirical
#'   likelihood methods for linear models yield the exact same fitted values and
#'   residuals as those produced by `lm()`. Users may resort to other functions
#'   such as `model.frame()` with `lm()` and `glm()` to extract the variables
#'   and metadata.
#' @srrstatsNA {RE6.3} Empirical likelihood methods are almost never used in a
#'   prediction or a forecast context, in part due to the convex hull
#'   constraint. The package does not provide any forecast method (e.g.,
#'   `predict()`), as well as associated `plot()` methods.
#' @noRd
NULL
