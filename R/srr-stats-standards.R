#' srr_stats
#'
#' All of the following standards initially have `@srrstatsTODO` tags.
#' These may be moved at any time to any other locations in your code.
#' Once addressed, please modify the tag from `@srrstatsTODO` to `@srrstats`,
#' or `@srrstatsNA`, ensuring that references to every one of the following
#' standards remain somewhere within your code.
#' (These comments may be deleted at any time.)
#'
#' @srrstatsVerbose TRUE
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





#' @srrstatsTODO {G2.11} *Software should ensure that `data.frame`-like tabular objects which have columns which do not themselves have standard class attributes (typically, `vector`) are appropriately processed, and do not error without reason. This behaviour should be tested. Again, columns created by the [`units` package](https://github.com/r-quantities/units/) provide a good test case.*

#' @srrstatsTODO {G2.12} *Software should ensure that `data.frame`-like tabular objects which have list columns should ensure that those columns are appropriately pre-processed either through being removed, converted to equivalent vector columns where appropriate, or some other appropriate treatment such as an informative error. This behaviour should be tested.*



#' @srrstatsTODO {G5.4} **Correctness tests** *to test that statistical algorithms produce expected results to some fixed test data sets (potentially through comparisons using binding frameworks such as [RStata](https://github.com/lbraglia/RStata)).*
#' @srrstatsTODO {G5.4a} *For new methods, it can be difficult to separate out correctness of the method from the correctness of the implementation, as there may not be reference for comparison. In this case, testing may be implemented against simple, trivial cases or against multiple implementations such as an initial R implementation compared with results from a C/C++ implementation.*
#' @srrstatsTODO {G5.4b} *For new implementations of existing methods, correctness tests should include tests against previous implementations. Such testing may explicitly call those implementations in testing, preferably from fixed-versions of other software, or use stored outputs from those where that is not possible.*
#' @srrstatsTODO {G5.4c} *Where applicable, stored values may be drawn from published paper outputs when applicable and where code from original implementations is not available*



#' @srrstatsTODO {RE1.3} *Regression Software which passes or otherwise transforms aspects of input data onto output structures should ensure that those output structures retain all relevant aspects of input data, notably including row and column names, and potentially information from other `attributes()`.*
#' @srrstatsTODO {RE1.3a} *Where otherwise relevant information is not transferred, this should be explicitly documented.*


#' @srrstatsTODO {RE2.0} *Regression Software should document any transformations applied to input data, for example conversion of label-values to `factor`, and should provide ways to explicitly avoid any default transformations (with error or warning conditions where appropriate).*


#' offset!!
#' @srrstatsTODO {RE2.3} *Where applicable, Regression Software should enable data to be centred (for example, through converting to zero-mean equivalent values; or to z-scores) or offset (for example, to zero-intercept equivalent values) via additional parameters, with the effects of any such parameters clearly documented and tested.*


#' @srrstatsTODO {RE4.8} *Response variables, and associated "metadata" where applicable.*
#' @srrstatsTODO {RE4.9} *Modelled values of response variables.*


#' @srrstatsTODO {RE4.13} *Predictor variables, and associated "metadata" where applicable.*







#' @srrstatsTODO {RE7.1} *Tests with noiseless, exact relationships between predictor (independent) and response (dependent) data.*
#' @srrstatsTODO {RE7.1a} *In particular, these tests should confirm that model fitting is at least as fast or (preferably) faster than testing with equivalent noisy data (see RE2.4b).*

#' @srrstatsTODO {RE7.2} Demonstrate that output objects retain aspects of input data such as row or case names (see **RE1.3**).
#' @srrstatsTODO {RE7.3} Demonstrate and test expected behaviour when objects returned from regression software are submitted to the accessor methods of **RE4.2**--**RE4.7**.
#' @srrstatsTODO {RE7.4} Extending directly from **RE4.15**, where appropriate, tests should demonstrate and confirm that forecast errors, confidence intervals, or equivalent values increase with forecast horizons.
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
#' @srrstatsNA {G5.10} Unit tests take less then 30 seconds in all platforms we
#'   tested. We avoid long running tests in compliance with the CRAN guidelines.
#'   We believe that the separate flags are not necessary for the package.
#' @srrstatsNA {G5.11} Unit tests do not require large data sets.
#' @srrstatsNA {G5.11a} Unit tests do not download additional data sets.
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
#' @srrstatsNA {RE4.14, RE4.15, RE4.16} Empirical likelihood methods (and thus
#'   the package) are inherently not suitable for prediction, extrapolation, or
#'   forecasting. Especially, extrapolation can often violate the convex hull
#'   constraint, leading to an invalid empirical likelihood value. For this
#'   reason, the package does provide `predict()` methods.
#' @srrstatsNA {RE6.1} All `plot()` methods in the package are S4 generics.
#' @srrstatsNA {RE6.2} No function in the package produces `fitted` values.
#' @srrstatsNA {RE6.3} Empirical likelihood methods are almost never used in a
#'   prediction or a forecast context, in part due to the convex hull
#'   constraint. The package does not provide any forecast method (e.g.,
#'   `predict()`), as well as associated `plot()` methods.
#' @noRd
NULL
