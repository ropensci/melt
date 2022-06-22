#' Empirical likelihood displacement
#'
#' Computes empirical likelihood displacement for model diagnostics and outlier
#'   detection.
#'
#' @param object A fitted \linkS4class{EL} object.
#' @param control An object of class \linkS4class{ControlEL} constructed by
#'   \code{\link{el_control}}.
#' @details Let \eqn{L(\theta)} be the empirical log-likelihood function based
#'   on the full sample with \eqn{n} observations. The maximum empirical
#'   likelihood estimate is denoted by \eqn{\hat{\theta}}. Consider a reduced
#'   sample with the \eqn{i}th observation deleted and the corresponding
#'   estimate \eqn{\hat{\theta}_{(i)}}. The empirical likelihood displacement is
#'   defined by
#'   \deqn{\textnormal{ELD}_i = 2\{L(\hat{\theta}) - L(\hat{\theta}_{(i)})\}.}
#'   If \eqn{\textnormal{ELD}_i } is large, then the \eqn{i}th observation is an
#'   influential point and can be inspected as a possible outlier. \code{eld}
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
setGeneric("eld", function(object, control = el_control()) {
  standardGeneric("eld")
})


#' Model coefficients
#'
#' Extracts maximum empirical likelihood estimates from a model.
#'
#' @param object A fitted \linkS4class{EL} object.
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
#' @param object An object of class \linkS4class{EL}.
#' @param parm A specification of which parameters are to be given confidence
#'   intervals, either a vector of numbers or a vector of names. If missing, all
#'   parameters are considered.
#' @param level A single numeric for the confidence level required. Defaults to
#'   \code{0.95}.
#' @param cv A single numeric for the critical value for calibration of
#'   empirical likelihood ratio statistic. Defaults to \code{NULL} and set to
#'   \code{qchisq(level, 1L)}. If non-\code{NULL}, \code{level} is ignored.
#' @param control An object of class \linkS4class{ControlEL} constructed by
#'   \code{\link{el_control}}.
#' @param ... Not used.
#' @importFrom stats qchisq
#' @return A matrix with columns giving lower and upper confidence limits for
#'  each parameter. In contrast to other methods that rely on studentization,
#'  the lower and upper limits obtained from empirical likelihood do not
#'  correspond to the \code{(1 - level) / 2} and \code{1 - (1 - level) / 2} in
#'  \%, respectively.
#' @references Kim, E., MacEachern, S., and Peruggia, M., (2021),
#'   "Empirical Likelihood for the Analysis of Experimental Designs,"
#'   \href{https://arxiv.org/abs/2112.09206}{arxiv:2112.09206}.
#' @references Owen, Art. 1990. “Empirical Likelihood Ratio Confidence Regions.”
#'   The Annals of Statistics 18 (1): 90–120. \doi{10.1214/aos/1176347494}.
#' @seealso \link{confreg}, \link{el_control}, \link{elt}
#' @usage NULL
#' @examples
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
#' @param object An object of class \linkS4class{EL}.
#' @param parm A specification of which parameters are to be given a confidence
#'   region, either a vector of numbers or a vector of names. It should be a
#'   vector of length two of the form \code{c(x, y)}. If missing, the first two
#'   parameter in \code{object} are considered.
#' @param level A single numeric for the confidence level required. Defaults to
#'   \code{0.95}.
#' @param cv A single numeric for the critical value for calibration of
#'   empirical likelihood ratio statistic. Defaults to NULL and set to
#'   \code{qchisq(level, 2L)}. If non-\code{NULL}, \code{level} is ignored.
#' @param npoints A single integer for the number of boundary points to compute.
#'   Defaults to \code{50}.
#' @param control An object of class \linkS4class{ControlEL} constructed by
#'   \code{\link{el_control}}.
#' @importFrom stats qchisq
#' @return An object of class \linkS4class{ConfregEL}.
#' @references Owen, Art. 1990. “Empirical Likelihood Ratio Confidence Regions.”
#'   The Annals of Statistics 18 (1): 90–120. \doi{10.1214/aos/1176347494}.
#' @seealso \link{confint}, \link{el_control}, \link{elt}, \link{plot}
#' @usage NULL
#' @examples
#' fit <- el_lm(mpg ~ ., data = mtcars)
#' confreg(fit, parm = c(2, 3), level = 0.95, cv = qchisq(0.95, 2L))
#' @exportMethod confreg
setGeneric("confreg", function(object, parm, level = 0.95, cv = NULL,
                               npoints = 50L, control = el_control()) {
  standardGeneric("confreg")
})


#' Empirical log-likelihood
#'
#' Extracts empirical log-likelihood from a model evaluated at the estimated
#'   coefficients.
#'
#' @param object An object of class \linkS4class{EL}.
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
#' @param ... Further graphical parameters (see \code{\link[graphics]{par}}).
#' @seealso \link{confreg}, \link{eld}
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
#' @param digits The number of significant digits to be passed to
#'   \code{\link[base]{format}}.
#' @param signif.stars A single logical. If \code{TRUE}, ‘significance stars’
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
#' @param object An object that inherit from class \linkS4class{EL}.
#' @param ... Not used.
#' @usage NULL
#' @examples
#' data("airquality")
#' x <- airquality$Wind
#' w <- airquality$Day
#' fit <- el_mean(x, par = 10, weights = w)
#' weights(fit)
#' @exportMethod weights
setGeneric("weights", function(object, ...) standardGeneric("weights"))


setGeneric("getMethodEL", function(x) standardGeneric("getMethodEL"))
setMethod("getMethodEL", "EL", function(x) {x@method})

setGeneric("getDataMatrix", function(x) standardGeneric("getDataMatrix"))
setMethod("getDataMatrix", "EL", function(x) {x@data})

setGeneric("getWeights", function(x) standardGeneric("getWeights"))
setMethod("getWeights", "EL", function(x) {x@weights})
