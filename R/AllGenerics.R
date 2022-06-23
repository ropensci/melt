#' Empirical likelihood displacement
#'
#' Computes empirical likelihood displacement for model diagnostics and outlier
#'   detection.
#'
#' @param object A fitted \linkS4class{EL} object.
#' @param control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @details Let \eqn{L(\theta)} be the empirical log-likelihood function based
#'   on the full sample with \eqn{n} observations. The maximum empirical
#'   likelihood estimate is denoted by \eqn{\hat{\theta}}. Consider a reduced
#'   sample with the \eqn{i}th observation deleted and the corresponding
#'   estimate \eqn{\hat{\theta}_{(i)}}. The empirical likelihood displacement is
#'   defined by
#'   \deqn{\textnormal{ELD}_i = 2\{L(\hat{\theta}) - L(\hat{\theta}_{(i)})\}.}
#'   If \eqn{\textnormal{ELD}_i } is large, then the \eqn{i}th observation is an
#'   influential point and can be inspected as a possible outlier. `eld`
#'   computes \eqn{\textnormal{ELD}_i } for \eqn{i = 1, \dots, n }.
#' @return An object of class \linkS4class{ELD}.
#' @references Lazar, Nicole A. 2005. “Assessing the Effect of Individual Data
#'   Points on Inference From Empirical Likelihood.” Journal of Computational
#'   and Graphical Statistics 14 (3): 626–42.
#'   \doi{10.1198/106186005X59568}.
#' @references Zhu, H., J. G. Ibrahim, N. Tang, and H. Zhang. 2008. “Diagnostic
#'   Measures for Empirical Likelihood of General Estimating Equations.”
#'   Biometrika 95 (2): 489–507.
#'   \doi{10.1093/biomet/asm094}.
#' @seealso \link{el_control}, \link{el_eval}, \link{plot}
#' @usage NULL
#' @examples
#' x <- rnorm(10L)
#' y <- 10
#' fit <- el_mean(c(x, y), 0)
#' eld(fit)
#' @exportMethod eld
setGeneric("eld", function(object, ...) {
  standardGeneric("eld")
})


#' Empirical likelihood test
#'
#' Tests a linear hypothesis for objects that inherit from class
#'   \linkS4class{EL}.
#'
#' @param object An object that inherit from class \linkS4class{EL}, including
#'   \linkS4class{CEL}, \linkS4class{LM}, and \linkS4class{GLM}.
#' @param rhs A numeric vector or a column matrix for the right-hand-side of
#'   hypothesis, with as many entries as the rows in `lhs`. Defaults to `NULL`.
#'   See ‘Details’.
#' @param lhs A numeric matrix or a vector (treated as a row matrix) for the
#'   left-hand-side of hypothesis. Each row gives a linear combination of the
#'   parameters in `object`. The number of columns should be equal to the number
#'   of parameters. Defaults to `NULL`. See ‘Details’.
#' @param alpha A single numeric for the significance level. Defaults to `0.05`.
#' @param calibrate A single character for the calibration method. Defaults to
#'   `"chisq"`. See ‘Details’.
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
#'   \enumerate{
#'   \item If both `rhs` and `lhs` are non-`NULL`, the constrained minimization
#'   is performed with the right-hand-side \eqn{r} and the left-hand-side
#'   \eqn{L} as
#'   \deqn{\inf_{\theta: L\theta = r} l(\theta).}
#'   \item If `rhs` is `NULL`, \eqn{r} is set to the zero vector as
#'   \eqn{\inf_{\theta: L\theta = 0} l(\theta)}.
#'   \item If `lhs` is `NULL`, \eqn{L} is set to the identity matrix and the
#'   problem reduces to evaluating at \eqn{r} as \eqn{l(r)}.
#'   }
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
#' @references Adimari, Gianfranco, and Annamaria Guolo. 2010.
#'   “A Note on the Asymptotic Behaviour of Empirical Likelihood Statistics.”
#'   Statistical Methods & Applications 19 (4): 463–76.
#'   \doi{10.1007/s10260-010-0137-9}.
#' @references Qin, Jing, and Jerry Lawless. 1995.
#'   “Estimating Equations, Empirical Likelihood and Constraints on Parameters.”
#'   Canadian Journal of Statistics 23 (2): 145–59. \doi{10.2307/3315441}.
#' @seealso [el_control()]
#' @usage NULL
#' @examples
#' n <- 100
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' y <- 1 + x1 + x2 + rnorm(n)
#' df <- data.frame(y, x1, x2)
#' fit <- el_lm(y ~ x1 + x2, df)
#' elt(fit, lhs = c(0, 1, -1))
#' elt(fit, lhs = c(0, 1, 1), rhs = 2)
#'
#' # test of no treatment effect
#' data("clothianidin")
#' lhs <- matrix(c(
#'   1, -1, 0, 0,
#'   0, 1, -1, 0,
#'   0, 0, 1, -1
#' ), byrow = TRUE, nrow = 3)
#' fit2 <- el_lm(clo ~ -1 + trt, clothianidin)
#' elt(fit2, lhs = lhs)
#' @exportMethod elt
setGeneric("elt", function(object, ...) {
  standardGeneric("elt")
})


#' Model coefficients
#'
#' Extracts maximum empirical likelihood estimates from a model.
#'
#' @param object An object that inherit from class \linkS4class{EL}, including
#'   \linkS4class{CEL}, \linkS4class{LM}, and \linkS4class{GLM}.
#' @param ... Not used.
#' @return A numeric vector of the maximum empirical likelihood estimates.
#' @examples
#' fit <- el_lm(mpg ~ wt, data = mtcars)
#' coef(fit)
#' @usage NULL
#' @exportMethod coef
setGeneric("coef", function(object, ...) standardGeneric("coef"))


#' Confidence interval for model parameters
#'
#' Computes confidence intervals for one or more parameters in a fitted model.
#' The package \pkg{melt} adds a method for objects inheriting from class
#' \linkS4class{EL}.
#'
#' @param object An object that inherit from class \linkS4class{EL}, including
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
#' @param ... Not used.
#' @importFrom stats qchisq
#' @return A matrix with columns giving lower and upper confidence limits for
#'  each parameter. In contrast to other methods that rely on studentization,
#'  the lower and upper limits obtained from empirical likelihood do not
#'  correspond to the `(1 - level) / 2` and `1 - (1 - level) / 2` in %,
#'  respectively.
#' @references Owen, Art. 1990. “Empirical Likelihood Ratio Confidence Regions.”
#'   The Annals of Statistics 18 (1): 90–120. \doi{10.1214/aos/1176347494}.
#' @seealso \link{confreg}, \link{el_control}, \link{elt}
#' @usage NULL
#' @examples
#' fit <- el_lm(mpg ~ ., data = mtcars)
#' confint(fit, parm = c(2, 3))
#' @exportMethod confint
setGeneric(
  "confint",
  function(object, parm, level = 0.95, ...) standardGeneric("confint")
)


#' Confidence region for model parameters
#'
#' Computes boundary points of a two-dimensional confidence region for model
#'   parameters.
#'
#' @param object An object that inherit from class \linkS4class{EL}, including
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
#' @importFrom stats qchisq
#' @return An object of class \linkS4class{ConfregEL}.
#' @references Owen, Art. 1990. “Empirical Likelihood Ratio Confidence Regions.”
#'   The Annals of Statistics 18 (1): 90–120. \doi{10.1214/aos/1176347494}.
#' @seealso [confint()], [el_control()], [elt()], [plot()]
#' @usage NULL
#' @examples
#' fit <- el_lm(mpg ~ ., data = mtcars)
#' confreg(fit, parm = c(2, 3), level = 0.95, cv = qchisq(0.95, 2))
#' @exportMethod confreg
setGeneric("confreg", function(object, ...) {
  standardGeneric("confreg")
})


#' Empirical log-likelihood
#'
#' Extracts empirical log-likelihood from a model evaluated at the estimated
#'   coefficients.
#'
#' @param object An object that inherit from class \linkS4class{EL}, including
#'   \linkS4class{CEL}, \linkS4class{LM}, and \linkS4class{GLM}.
#' @param ... Not used.
#' @return An object of class \linkS4class{logLikEL}.
#' @examples
#' fit <- el_lm(mpg ~ wt, data = mtcars)
#' logLik(fit)
#' @usage NULL
#' @exportMethod logLik
setGeneric("logLik", function(object, ...) standardGeneric("logLik"))


#' Plot methods
#'
#' Provides plot methods for objects that inherit from class \linkS4class{EL}.
#'
#' @param x An object to be plotted.
#' @param y Not used.
#' @param ... Further graphical parameters (see [`par`]).
#' @seealso [confreg()], [eld()]
#' @usage NULL
#' @examples
#' # model
#' fit <- el_lm(hp ~ wt, data = mtcars)
#'
#' # confidence region
#' out1 <- confreg(fit, npoints = 500)
#' plot(out1)
#'
#' # empirical likelihood displacement
#' out2 <- eld(fit)
#' plot(out2)
#' @exportMethod plot
setGeneric("plot", function(x, y, ...) standardGeneric("plot"))


#' Print methods
#'
#' Provides print methods for the objects from the package \pkg{melt}.
#'
#' @param x An object to be printed.
#' @param ... Further arguments passed to other methods.
#' @param digits The number of significant digits to be passed to [format()].
#' @param signif.stars A single logical. If `TRUE`, ‘significance stars’
#'   are printed for each coefficient.
#' @usage NULL
#' @examples
#' # model
#' fit <- el_lm(mpg ~ wt, data = mtcars)
#'
#' # output
#' fit
#'
#' # empirical log-likelihood
#' logLik(fit)
#'
#' # output summary
#' summary(fit)
#'
#' # equivalent results by testing each parameter separately
#' elt(fit, lhs = c(1, 0))
#' elt(fit, lhs = c(0, 1))
#' @exportMethod print
setGeneric("print", function(x, ...) standardGeneric("print"))


#' Summary methods
#'
#' Provides summary methods for objects that inherit from class
#'   \linkS4class{EL}.
#'
#' @param object An object for which a summary is desired.
#' @param ... Additional arguments affecting the summary produced.
#' @usage NULL
#' @examples
#' fit <- el_lm(mpg ~ wt, data = mtcars)
#' summary(fit)
#' @exportMethod summary
setGeneric("summary", function(object, ...) standardGeneric("summary"))


#' Model weights
#'
#' Extracts model weights from objects that inherit from class
#'   \linkS4class{EL}. The weights are rescaled to up to the total number of
#'   observations in the fitting procedure.
#'
#' @param object An object that inherit from class \linkS4class{EL}, including
#'   \linkS4class{CEL}, \linkS4class{LM}, and \linkS4class{GLM}.
#' @param ... Not used.
#' @usage NULL
#' @return A numeric vector the rescaled weights.
#' @examples
#' data("airquality")
#' x <- airquality$Wind
#' w <- airquality$Day
#' fit <- el_mean(x, par = 10, weights = w)
#' weights(fit)
#' @exportMethod weights
setGeneric("weights", function(object, ...) standardGeneric("weights"))


setGeneric("getMethodEL", function(x) standardGeneric("getMethodEL"))
setMethod("getMethodEL", "EL", function(x) {
  x@method
})

setGeneric("getDataMatrix", function(x) standardGeneric("getDataMatrix"))
setMethod("getDataMatrix", "EL", function(x) {
  x@data
})

setGeneric("getWeights", function(x) standardGeneric("getWeights"))
setMethod("getWeights", "EL", function(x) {
  x@weights
})
