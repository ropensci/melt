#' Chi-square statistic
#'
#' Extracts chi-square statistic from a model.
#'
#' @param object An object that inherit from \linkS4class{EL},
#'   \linkS4class{ELT}, \linkS4class{ELMT}, or \linkS4class{SummaryLM}.
#' @param ... Further arguments passed to methods.
#' @return The form of the value returned by [chisq()] depends on the class of
#'   its argument.
#' @seealso [pVal()]
#' @usage NULL
#' @examples
#' data("precip")
#' fit <- el_mean(precip, par = 40)
#' chisq(fit)
#' @exportMethod chisq
setGeneric("chisq", function(object, ...) standardGeneric("chisq"))


#' Model coefficients
#'
#' Extracts maximum empirical likelihood estimates from a model.
#'
#' @param object An object that inherit from \linkS4class{EL}, including
#'   \linkS4class{CEL}, \linkS4class{LM}, and \linkS4class{GLM}.
#' @param ... Further arguments passed to methods.
#' @return A numeric vector of the maximum empirical likelihood estimates.
#' @usage NULL
#' @examples
#' data("mtcars")
#' fit <- el_lm(mpg ~ wt, data = mtcars)
#' coef(fit)
#' @exportMethod coef
setGeneric("coef", function(object, ...) standardGeneric("coef"))


#' Confidence interval for model parameters
#'
#' Computes confidence intervals for one or more parameters in a model.
#'
#' @param object An object that inherit from \linkS4class{EL}, including
#'   \linkS4class{CEL}, \linkS4class{LM}, and \linkS4class{GLM}.
#' @param parm A specification of which parameters are to be given confidence
#'   intervals, either a vector of numbers or a vector of names. If missing, all
#'   parameters are considered.
#' @param level A single numeric for the confidence level required. Defaults to
#'   `0.95`.
#' @param cv A single numeric for the critical value for calibration of
#'   empirical likelihood ratio statistic. Defaults to `NULL` and set to
#'   `qchisq(level, 1L)`. If non-`NULL`, `level` is ignored.
#' @param control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @return A matrix with columns giving lower and upper confidence limits for
#'  each parameter. In contrast to other methods that rely on studentization,
#'  the lower and upper limits obtained from empirical likelihood do not
#'  correspond to the `(1 - level) / 2` and `1 - (1 - level) / 2` in %,
#'  respectively.
#' @references Owen A (1990).
#'   “Empirical Likelihood Ratio Confidence Regions.”
#'   The Annals of Statistics, 18(1), 90–120.
#'   \doi{10.1214/aos/1176347494}.
#' @seealso [confreg()], [el_control()], [elt()]
#' @usage NULL
#' @examples
#' data("mtcars")
#' fit <- el_lm(mpg ~ ., data = mtcars)
#' confint(fit, parm = c(2, 3))
#' @exportMethod confint
setGeneric("confint", function(object, parm, level = 0.95, ...)
  standardGeneric("confint")
)


#' Confidence region for model parameters
#'
#' Computes boundary points of a two-dimensional confidence region for model
#'   parameters.
#'
#' @param object An object that inherit from \linkS4class{EL}, including
#'   \linkS4class{CEL}, \linkS4class{LM}, and \linkS4class{GLM}.
#' @param parm A specification of which parameters are to be given a confidence
#'   region, either a vector of numbers or a vector of names. It should be a
#'   vector of length two of the form `c(x, y)`. If missing, the first two
#'   parameter in `object` are considered.
#' @param level A single numeric for the confidence level required. Defaults to
#'   `0.95`. It is ignored if `cv` is non-`NULL`.
#' @param cv A single numeric for the critical value for calibration of
#'   empirical likelihood ratio statistic. Defaults to NULL and set to
#'   `qchisq(level, 2L)`. It must be compatible with the `th` value in
#'   `control`.
#' @param npoints A single integer for the number of boundary points to compute.
#'   Defaults to `50`.
#' @param control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @return An object of class \linkS4class{ConfregEL}.
#' @references Owen A (1990).
#'   “Empirical Likelihood Ratio Confidence Regions.”
#'   The Annals of Statistics, 18(1), 90–120.
#'   \doi{10.1214/aos/1176347494}.
#' @seealso [confint()], [el_control()], [elt()], [plot()]
#' @usage NULL
#' @examples
#' data("mtcars")
#' fit <- el_lm(mpg ~ wt + qsec, data = mtcars)
#' cr <- confreg(fit, parm = c(2, 3), cv = qchisq(0.90, 2))
#' plot(cr)
#' @exportMethod confreg
setGeneric("confreg", function(object,
                               parm,
                               level = 0.95,
                               cv = NULL,
                               npoints = 50L,
                               control = el_control()) {
  standardGeneric("confreg")
})


#' Convergence check
#'
#' Extracts convergence status from a model.
#'
#' @param object An object that inherit from \linkS4class{EL}, including
#'   \linkS4class{CEL}, \linkS4class{LM}, and \linkS4class{GLM}.
#' @param ... Further arguments passed to methods.
#' @return A single logical.
#' @usage NULL
#' @examples
#' ## Convergence check for the overall model test
#' data("mtcars")
#' fit <- el_lm(mpg ~ ., data = mtcars)
#' conv(fit)
#' @exportMethod conv
setGeneric("conv", function(object, ...) standardGeneric("conv"))


#' Critical value
#'
#' Extracts critical value from a model.
#'
#' @param object An object that inherit from \linkS4class{ELT} or
#'   \linkS4class{ELMT}.
#' @param ... Further arguments passed to methods.
#' @return A single numeric.
#' @usage NULL
#' @examples
#' ## F-calibrated critical value
#' set.seed(533414)
#' x <- rnorm(100)
#' fit <- el_mean(x, 0)
#' elt <- elt(fit, rhs = 0.3, calibrate = "f")
#' critVal(elt)
#' @exportMethod critVal
setGeneric("critVal", function(object, ...) standardGeneric("critVal"))


#' Empirical likelihood displacement
#'
#' Computes empirical likelihood displacement for model diagnostics and outlier
#'   detection.
#'
#' @param object An object that inherit from \linkS4class{EL}, including
#'   \linkS4class{CEL}, \linkS4class{LM}, and \linkS4class{GLM}.
#' @param control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @details Let \eqn{L(\theta)} be the empirical log-likelihood function based
#'   on the full sample with \eqn{n} observations. The maximum empirical
#'   likelihood estimate is denoted by \eqn{\hat{\theta}}. Consider a reduced
#'   sample with the \eqn{i}th observation deleted and the corresponding
#'   estimate \eqn{\hat{\theta}_{(i)}}. The empirical likelihood displacement is
#'   defined by
#'   \deqn{\textrm{ELD}_i = 2\{L(\hat{\theta}) - L(\hat{\theta}_{(i)})\}.}
#'   If \eqn{\textrm{ELD}_i } is large, then the \eqn{i}th observation is an
#'   influential point and can be inspected as a possible outlier. `eld`
#'   computes \eqn{\textrm{ELD}_i } for \eqn{i = 1, \dots, n }.
#' @return An object of class \linkS4class{ELD}.
#' @references Lazar NA (2005).
#'   “Assessing the Effect of Individual Data Points on Inference From Empirical
#'   Likelihood.”
#'   Journal of Computational and Graphical Statistics, 14(3), 626–642.
#'   \doi{10.1198/106186005X59568}.
#' @references Zhu H, Ibrahim JG, Tang N, Zhang H (2008).
#'   “Diagnostic Measures for Empirical Likelihood of General Estimating
#'   Equations.” Biometrika, 95(2), 489–507.
#'   \doi{10.1093/biomet/asm094}.
#' @seealso [el_control()], [el_eval()], [plot()]
#' @usage NULL
#' @examples
#' data("precip")
#' fit <- el_mean(precip, par = 30)
#' eld <- eld(fit)
#' plot(eld)
#' @exportMethod eld
setGeneric("eld", function(object, control = el_control()) {
  standardGeneric("eld")
})


#' Empirical likelihood multiple tests
#'
#' Tests multiple linear hypotheses simultaneously.
#'
#' @param object An object that inherit from \linkS4class{EL}, including
#'   \linkS4class{CEL}, \linkS4class{LM}, and \linkS4class{GLM}.
#' @param rhs A numeric vector (column matrix) or a list of numeric vectors for
#'   the right-hand sides of hypotheses. Defaults to `NULL`. See ‘Details’.
#' @param lhs A numeric matrix or a list of numeric matrices for the left-hand
#'   sides of hypothesis. Each row of the matrices gives a linear combination of
#'   the parameters in `object`. The number of columns should be equal to the
#'   number of parameters. Defaults to `NULL`. See ‘Details’.
#' @param alpha A single numeric for the overall significance level. Defaults to
#'   `0.05`.
#' @param control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @details [elmt()] tests multiple hypotheses simultaneously. Each hypothesis
#'   corresponds to the constrained empirical likelihood ratio described in
#'   \linkS4class{CEL}. `rhs` and `lhs` cannot be both `NULL`. The right-hand
#'   side and left-hand side of each hypothesis must be specified as described
#'   in [elt()].
#'
#'   For specifying linear contrasts more conveniently, `rhs` and `lhs` also
#'   take a numeric vector and a numeric matrix, respectively. Each element of
#'   `rhs` and each row of `lhs` correspond to a contrast (hypothesis).
#'
#'   The vector of empirical likelihood ratio statistics asymptotically follows
#'   a multivariate chi-square distribution under the complete null hypothesis.
#'   The multiple testing procedure asymptotically controls the family-wise
#'   error rate at the level `alpha`. Based on the distribution of the maximum
#'   of the test statistics, the adjusted p-values are estimated by Monte Carlo
#'   simulation.
#' @return An object of class of \linkS4class{ELMT}.
#' @references Kim E, MacEachern S, Peruggia M (2021).
#'   “Empirical Likelihood for the Analysis of Experimental Designs.”
#'   arxiv:2112.09206. URL <https://arxiv.org/abs/2112.09206>.
#' @seealso [el_control()], [elt()]
#' @usage NULL
#' @examples
#' ## Example 1: bivariate mean (list `rhs` & no `lhs`)
#' data("women")
#' fit <- el_mean(women, par = c(65, 135))
#' rhs <- list(c(64, 133), c(66, 140))
#' set.seed(143)
#' elmt(fit, rhs = rhs)
#'
#' ## Example 2: pairwise comparison (no `rhs` & matrix `lhs`)
#' data("clothianidin")
#' fit2 <- el_lm(clo ~ -1 + trt, clothianidin)
#' lhs <- matrix(c(
#'   1, -1, 0, 0,
#'   0, 1, -1, 0,
#'   0, 0, 1, -1
#' ), byrow = TRUE, nrow = 3)
#' set.seed(629)
#' elmt(fit2, lhs = lhs)
#'
#' ## Example 3: arbitrary hypotheses (list `rhs` & list `lhs`)
#' data("mtcars")
#' fit <- el_lm(mpg ~ wt + qsec, data = mtcars)
#' lhs <- list(rbind(c(1, 4, 0)), rbind(c(0, 1, 0), c(0, 0, 1)))
#' rhs <- list(0, c(-6, 1))
#' elmt(fit, rhs = rhs, lhs = lhs)
#' @exportMethod elmt
setGeneric("elmt", function(object,
                            rhs = NULL,
                            lhs = NULL,
                            alpha = 0.05,
                            control = el_control()) {
  standardGeneric("elmt")
})


#' Empirical likelihood test
#'
#' Tests a linear hypothesis.
#'
#' @param object An object that inherit from \linkS4class{EL}, including
#'   \linkS4class{CEL}, \linkS4class{LM}, and \linkS4class{GLM}.
#' @param rhs A numeric vector or a column matrix for the right-hand side of
#'   hypothesis, with as many entries as the rows in `lhs`. Defaults to `NULL`.
#'   See ‘Details’.
#' @param lhs A numeric matrix or a vector (treated as a row matrix) for the
#'   left-hand side of hypothesis. Each row gives a linear combination of the
#'   parameters in `object`. The number of columns should be equal to the number
#'   of parameters. Defaults to `NULL`. See ‘Details’.
#' @param alpha A single numeric for the significance level. Defaults to `0.05`.
#' @param calibrate A single character for the calibration method. It is
#'   case-insensitive and must be one of `"chisq"`, `"boot"`, or `"f"`.
#'   Defaults to `"chisq"`. See ‘Details’.
#' @param control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @details [elt()] performs the constrained minimization of \eqn{l(\theta)}
#'   described in \linkS4class{CEL}. `rhs` and `lhs` cannot be both `NULL`. For
#'   non-`NULL` `lhs`, it is required that `lhs` have full row rank
#'   \eqn{q \leq p} and \eqn{p} be equal to the number of parameters in the
#'   `object`.
#'
#'   Depending on the specification of `rhs` and `lhs`, we have the following
#'   three cases:
#'   1. If both `rhs` and `lhs` are non-`NULL`, the constrained minimization
#'   is performed with the right-hand side \eqn{r} and the left-hand side
#'   \eqn{L} as
#'   \deqn{\inf_{\theta: L\theta = r} l(\theta).}
#'   1. If `rhs` is `NULL`, \eqn{r} is set to the zero vector as
#'   \eqn{\inf_{\theta: L\theta = 0} l(\theta)}.
#'   1. If `lhs` is `NULL`, \eqn{L} is set to the identity matrix and the
#'   problem reduces to evaluating at \eqn{r} as \eqn{l(r)}.
#'
#'   `calibrate` specifies the calibration method used. Three methods are
#'   available: `"chisq"` (chi-square calibration), `"boot"` (bootstrap
#'   calibration), and `"f"` (\eqn{F} calibration). `"boot"` is applicable only
#'   when `lhs` is `NULL`. The `nthreads`, `seed`, and `B` slots in `control`
#'   apply to the bootstrap procedure. `"f"` is applicable only to the mean
#'   parameter when `lhs` is `NULL`.
#' @return An object of class of \linkS4class{ELT}. If `lhs` is non-`NULL`, the
#'   `optim` slot corresponds to that of \linkS4class{CEL}. Otherwise, it
#'   corresponds to that of \linkS4class{EL}.
#' @references Adimari G, Guolo A (2010).
#'   “A Note on the Asymptotic Behaviour of Empirical Likelihood Statistics.”
#'   Statistical Methods & Applications, 19(4), 463–476.
#'   \doi{10.1007/s10260-010-0137-9}.
#' @references Qin J, Lawless J (1995).
#'   “Estimating Equations, Empirical Likelihood and Constraints on Parameters.”
#'   Canadian Journal of Statistics, 23(2), 145–159. \doi{10.2307/3315441}.
#' @seealso [el_control()], [elmt()]
#' @usage NULL
#' @examples
#' ## F calibration for the mean
#' set.seed(533414)
#' x <- rnorm(100)
#' fit <- el_mean(x, 0)
#' elt(fit, rhs = 0.3, calibrate = "f")
#'
#' ## Test of no treatment effect
#' data("clothianidin")
#' lhs <- matrix(c(
#'   1, -1, 0, 0,
#'   0, 1, -1, 0,
#'   0, 0, 1, -1
#' ), byrow = TRUE, nrow = 3)
#' fit2 <- el_lm(clo ~ -1 + trt, clothianidin)
#' elt(fit2, lhs = lhs)
#' @exportMethod elt
#' @srrstats {G5.1} `clothianidin` data set is exported.
setGeneric("elt", function(object,
                           rhs = NULL,
                           lhs = NULL,
                           alpha = 0.05,
                           calibrate = "chisq",
                           control = el_control()) {
  standardGeneric("elt")
})


#' Degrees of freedom
#'
#' Extracts degrees of freedom from a model.
#'
#' @param object An object that inherit from \linkS4class{EL},
#'   \linkS4class{ELT}, \linkS4class{logLikEL}, or \linkS4class{SummaryLM}.
#' @return A single integer.
#' @usage NULL
#' @examples
#' data("faithful")
#' fit <- el_mean(faithful, par = c(3.5, 70))
#' getDF(fit)
#' @exportMethod getDF
setGeneric("getDF", function(object) standardGeneric("getDF"))


#' Optimization results
#'
#' Extracts optimization results from a model.
#'
#' @param object An object that inherit from \linkS4class{EL} or
#'   \linkS4class{ELT}.
#' @param ... Further arguments passed to methods.
#' @return A list with the following optimization results:
#'   * `par` A numeric vector of the parameter value. See the documentation of
#'   \linkS4class{EL} and \linkS4class{CEL}.
#'   * `lambda` A numeric vector of the Lagrange multipliers.
#'   * `iterations` A single integer for the number of iterations performed.
#'   * `convergence` A single logical for the convergence status.
#' @seealso [sigTests()]
#' @usage NULL
#' @examples
#' data("precip")
#' fit <- el_mean(precip, par = 40)
#' getOptim(fit)
#' @exportMethod getOptim
setGeneric("getOptim", function(object, ...) standardGeneric("getOptim"))


#' Empirical log-likelihood
#'
#' Extracts empirical log-likelihood from a model.
#'
#' @param object An object that inherit from \linkS4class{EL} or
#'   \linkS4class{ELT}.
#' @param ... Further arguments passed to methods.
#' @return A single numeric.
#' @usage NULL
#' @examples
#' data("precip")
#' fit <- el_mean(precip, par = 40)
#' logL(fit)
#' @exportMethod logL
setGeneric("logL", function(object, ...) standardGeneric("logL"))


#' Maximum empirical log-likelihood
#'
#' Extracts empirical log-likelihood from a model evaluated at the estimated
#'   coefficients.
#'
#' @param object An object that inherit from \linkS4class{EL}.
#' @param ... Further arguments passed to methods.
#' @return An object of class \linkS4class{logLikEL}.
#' @usage NULL
#' @examples
#' data("precip")
#' fit <- el_mean(precip, par = 40)
#' logLik(fit)
#' @exportMethod logLik
setGeneric("logLik", function(object, ...) standardGeneric("logLik"))


#' Empirical log-likelihood ratio
#'
#' Extracts empirical log-likelihood ratio from a model.
#'
#' @param object An object that inherit from \linkS4class{EL} or
#'   \linkS4class{ELT}.
#' @param ... Further arguments passed to methods.
#' @return A single numeric.
#' @usage NULL
#' @examples
#' data("precip")
#' fit <- el_mean(precip, par = 40)
#' logLR(fit)
#' @exportMethod logLR
setGeneric("logLR", function(object, ...) standardGeneric("logLR"))


#' Number of observations in a model
#'
#' Extracts number of observations from a model.
#'
#' @param object An object that inherit from \linkS4class{EL}.
#' @param ... Further arguments passed to methods.
#' @return A single integer.
#' @usage NULL
#' @examples
#' data("precip")
#' fit <- el_mean(precip, par = 40)
#' nobs(fit)
#' @exportMethod nobs
setGeneric("nobs", function(object, ...) standardGeneric("nobs"))


#' Plot methods
#'
#' Provides plot methods for objects.
#'
#' @param x An object to be plotted.
#' @param y Not used.
#' @param ... Further graphical parameters (see [`par`]).
#' @seealso [confreg()], [eld()]
#' @usage NULL
#' @examples
#' ## Model
#' data("mtcars")
#' fit <- el_lm(hp ~ wt, data = mtcars)
#'
#' ## Confidence region
#' out1 <- confreg(fit, npoints = 500)
#' plot(out1)
#'
#' ## Empirical likelihood displacement
#' out2 <- eld(fit)
#' plot(out2)
#' @exportMethod plot
setGeneric("plot", function(x, y, ...) standardGeneric("plot"))


#' Print methods
#'
#' Provides print methods for objects.
#'
#' @param x An object to be printed.
#' @param digits A single integer for the number of significant digits to be
#'   passed to [format()].
#' @param signif.stars A single logical. If `TRUE`, ‘significance stars’
#'   are printed for each parameter.
#' @param ... Further arguments passed to methods.
#' @usage NULL
#' @examples
#' data("precip")
#' fit <- el_mean(precip, par = 40)
#' print(fit)
#' @exportMethod print
setGeneric("print", function(x, ...) standardGeneric("print"))


#' \eqn{p}-value
#'
#' Extracts \eqn{p}-value from a model.
#'
#' @param object An object that inherit from \linkS4class{EL},
#'   \linkS4class{ELT}, or \linkS4class{ELMT}.
#' @param ... Further arguments passed to methods.
#' @return The form of the value returned by [pVal()] depends on the class of
#'   its argument.
#' @seealso [chisq()]
#' @usage NULL
#' @examples
#' data("precip")
#' fit <- el_mean(precip, par = 40)
#' pVal(fit)
#' @exportMethod pVal
setGeneric("pVal", function(object, ...) standardGeneric("pVal"))


#' Significance tests
#'
#' Extracts the results of significance tests from a model.
#'
#' @param object An object that inherit from \linkS4class{LM} or
#'   \linkS4class{SummaryLM}.
#' @param ... Further arguments passed to methods.
#' @return The form of the value returned by [sigTests()] depends on the
#'   class of its argument.
#' @seealso [getOptim()]
#' @usage NULL
#' @examples
#' data("mtcars")
#' fit <- el_lm(mpg ~ ., data = mtcars)
#' sigTests(fit)
#' sigTests(summary(fit))
#' @exportMethod sigTests
setGeneric("sigTests", function(object, ...) standardGeneric("sigTests"))


#' Summary methods
#'
#' Provides summary methods for objects.
#'
#' @param object An object to be summarized.
#' @param ... Further arguments passed to methods.
#' @usage NULL
#' @examples
#' data("mtcars")
#' fit <- el_lm(mpg ~ wt, data = mtcars)
#' summary(fit)
#' @exportMethod summary
setGeneric("summary", function(object, ...) standardGeneric("summary"))


#' Model weights
#'
#' Extracts weights from model objects. The weights are re-scaled to up to the
#'   total number of observations in the fitting procedure.
#'
#' @param object An object that inherit from \linkS4class{EL}, including
#'   \linkS4class{CEL}, \linkS4class{LM}, and \linkS4class{GLM}.
#' @param ... Further arguments passed to methods.
#' @return A numeric vector of the re-scaled weights.
#' @references Glenn N, Zhao Y (2007).
#'   “Weighted Empirical Likelihood Estimates and Their Robustness Properties.”
#'   Computational Statistics & Data Analysis, 51(10), 5130–5141.
#'   \doi{10.1016/j.csda.2006.07.032}.
#' @usage NULL
#' @examples
#' data("airquality")
#' x <- airquality$Wind
#' w <- airquality$Day
#' fit <- el_mean(x, par = 10, weights = w)
#' weights(fit)
#' @exportMethod weights
setGeneric("weights", function(object, ...) standardGeneric("weights"))


setGeneric("getData", function(x) standardGeneric("getData"))
setMethod("getData", "EL", function(x) {
  x@data
})


setGeneric("getMethodEL", function(x) standardGeneric("getMethodEL"))
setMethod("getMethodEL", "EL", function(x) {
  x@method
})


setGeneric("getNumPar", function(x) standardGeneric("getNumPar"))
setMethod("getNumPar", "EL", function(x) {
  x@npar
})


setGeneric("getWeights", function(x) standardGeneric("getWeights"))
setMethod("getWeights", "EL", function(x) {
  x@weights
})
