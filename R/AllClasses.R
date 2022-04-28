#' S4 class \linkS4class{EL}
#'
#' S4 class for empirical likelihood.
#'
#' @slot optim A list with the following optimization results:
#'   \itemize{
#'   \item{\code{method } }{Character for method dispatch in internal
#'   functions.}
#'   \item{\code{par } }{Value at which empirical likelihood is evaluated.}
#'   \item{\code{lambda } }{Lagrange multiplier of the dual problem.}
#'   \item{\code{iterations } }{Number of iterations performed.}
#'   \item{\code{convergence } }{Convergence status.}
#'   }
#' @slot logp Log probabilities obtained from empirical likelihood.
#' @slot logl Empirical log-likelihood.
#' @slot loglr Empirical log-likelihood ratio.
#' @slot statistic Minus twice the empirical log-likelihood ratio.
#' @slot df Degrees of freedom of the statistic.
#' @slot pval \eqn{p}-value of the statistic.
#' @slot npar Number of parameters.
#' @slot weights Rescaled weights used for model fitting.
#' @slot dataMatrix Data matrix used for model fitting.
#' @slot coefficients Maximum empirical likelihood estimates of the parameters.
#' @examples
#' showClass("EL")
setClass("EL",
  slots = c(
    optim = "list", logp = "numeric", logl = "numeric", loglr = "numeric",
    statistic = "numeric", df = "integer", pval = "numeric", npar = "integer",
    weights = "numeric", dataMatrix = "matrix", coefficients = "numeric"
  ),
  prototype = list(
    optim = list(), logp = numeric(), logl = numeric(), loglr = numeric(),
    statistic = numeric(), df = 0L, pval = numeric(), npar = 0L,
    weights = numeric(), dataMatrix = matrix(NA_real_, nrow = 0L, ncol = 0L),
    coefficients = numeric()
  )
  # slots = c(
  #   optim = "list", logp = "numeric", logl = "numeric", loglr = "numeric",
  #   statistic = "numeric", df = "integer", pval = "numeric", npar = "integer",
  #   weights = "numeric", dataMatrix = "ANY", coefficients = "numeric"
  # ),
  # prototype = list(
  #   optim = list(), logp = numeric(), logl = numeric(), loglr = numeric(),
  #   statistic = numeric(), df = 0L, pval = numeric(), npar = 0L,
  #   weights = numeric(),
  #   coefficients = numeric()
  # )
)

#' S4 class \linkS4class{MinEL}
#'
#' S4 class for constrained empirical likelihood.
#'
#' @slot optim A list with the following optimization results:
#'   \itemize{
#'   \item{\code{method } }{Character for method dispatch in internal
#'   functions.}
#'   \item{\code{par } }{Value at which empirical likelihood is minimized.}
#'   \item{\code{lambda } }{Lagrange multiplier of the dual problem.}
#'   \item{\code{iterations } }{Number of iterations performed.}
#'   \item{\code{convergence } }{Convergence status.}
#'   }
#' @slot logp Log probabilities obtained from empirical likelihood.
#' @slot logl Empirical log-likelihood.
#' @slot loglr Empirical log-likelihood ratio.
#' @slot statistic Minus twice the empirical log-likelihood ratio.
#' @slot df Degrees of freedom of the statistic.
#' @slot pval \eqn{p}-value of the statistic.
#' @slot npar Number of parameters.
#' @slot weights Rescaled weights used for model fitting.
#' @slot dataMatrix Data matrix used for model fitting.
#' @slot coefficients Maximum empirical likelihood estimates of the parameters.
#' @examples
#' showClass("MinEL")
setClass("MinEL", contains = "EL")

#' S4 class \linkS4class{LM}
#'
#' S4 class for linear models with empirical likelihood.
#'
#' @slot parTests A list with the test results for each parameter:
#'   \itemize{
#'   \item{\code{statistic } }{Numeric vector of chi-squared statistics.}
#'   \item{\code{convergence } }{Logical vector. \code{TRUE} indicates
#'   convergence of the algorithm.}
#'   }
#' @slot misc A list with the following optimization results:
#'   \itemize{
#'   \item{\code{method } }{Character for method dispatch in internal
#'   functions.}
#'   \item{\code{par } }{Value at which empirical likelihood is minimized.}
#'   \item{\code{lambda } }{Lagrange multiplier of the dual problem.}
#'   \item{\code{iterations } }{Number of iterations performed.}
#'   \item{\code{convergence } }{Convergence status.}
#'   }
#' @examples
#' showClass("LM")
setClass("LM", contains = "MinEL", slots = c(parTests = "list", misc = "list"))

#' S4 class \linkS4class{GLM}
#'
#' S4 class for generalized linear models with empirical likelihood.
#'
#' @slot parTests A list with the test results for each parameter:
#'   \itemize{
#'   \item{\code{statistic } }{Numeric vector of chi-squared statistics.}
#'   \item{\code{convergence } }{Logical vector. \code{TRUE} indicates
#'   convergence of the algorithm.}
#'   }
#' @slot optim A list with the following optimization results:
#'   \itemize{
#'   \item{\code{method } }{Character for method dispatch in internal
#'   functions.}
#'   \item{\code{par } }{Value at which empirical likelihood is minimized.}
#'   \item{\code{lambda } }{Lagrange multiplier of the dual problem.}
#'   \item{\code{iterations } }{Number of iterations performed.}
#'   \item{\code{convergence } }{Convergence status.}
#'   }
#' @slot logp Log probabilities obtained from empirical likelihood.
#' @slot logl Empirical log-likelihood.
#' @slot loglr Empirical log-likelihood ratio.
#' @slot statistic Minus twice the empirical log-likelihood ratio.
#' @slot df Degrees of freedom of the statistic.
#' @slot pval \eqn{p}-value of the statistic.
#' @slot npar Number of parameters.
#' @slot weights Rescaled weights used for model fitting.
#' @slot dataMatrix Data matrix used for model fitting.
#' @slot coefficients Maximum empirical likelihood estimates of the parameters.
#' @examples
#' showClass("GLM")
setClass("GLM", contains = "LM")

#' S4 class \linkS4class{ConfregEL}
#'
#' S4 class for confidence region.
#'
#' @slot points A numeric matrix with two columns for boundary points of a
#'   confidence region.
#' @slot estimates A numeric vector of length two for parameter estimates.
#' @slot level A confidence level required.
#' @slot cv A critical value for calibration of empirical likelihood ratio
#'   statistic.
#' @slot pnames A character vector of length two for the name of parameters.
#' @examples
#' showClass("ConfregEL")
setClass("ConfregEL",
  slots = c(
    points = "matrix", estimates = "numeric", level = "numeric", cv = "numeric",
    pnames = "character"
  ),
  prototype = list(
    points = NULL, estimates = NA_real_, level = NA_real_, cv = NA_real_,
    pnames = NA_character_
  )
)

#' S4 class \linkS4class{ELD}
#'
#' S4 class for empirical likelihood displacement.
#'
#' @slot eld A numeric vector of empirical likelihood displacement values.
#' @examples
#' showClass("ELD")
setClass("ELD", slots = c(eld = "numeric"))

#' S4 class \linkS4class{SummaryLM}
#'
#' S4 class for a summary of \linkS4class{LM} objects.
#'
#' @slot statistic Minus twice the constrained empirical log-likelihood ratio
#'   for the overall test of the model.
#' @slot df Degrees of freedom of the statistic.
#' @slot convergence Convergence status of the minimization.
#' @slot parMatrix Numeric matrix of the test results of the parameters.
#' @slot weighted Logical for whether the given model is weighted or not.
#' @slot na.action Information returned by \code{\link[stats]{model.frame}} on
#'   the special handling of NAs.
#' @slot call Matched call.
#' @slot terms \code{\link[stats]{terms}} object used.
#' @slot aliased Named logical vector showing if the original coefficients are
#'   aliased.
#' @examples
#' showClass("SummaryLM")
setClass("SummaryLM", slots = c(
  statistic = "numeric", df = "integer", convergence = "logical",
  parMatrix = "matrix", weighted = "logical", na.action = "ANY", call = "ANY",
  terms = "ANY", aliased = "logical"
))

#' S4 class \linkS4class{ControlEL}
#'
#' S4 class for details of computation of empirical likelihood.
#'
#' @slot maxit Maximum number of iterations for the optimization with
#'   respect to \eqn{\theta}.
#' @slot maxit_l Maximum number of iterations for the optimization with
#'   respect to \eqn{\lambda}.
#' @slot tol Convergence tolerance denoted by \eqn{\epsilon}. The iteration
#'   stops when
#'   \deqn{\|P \nabla l(\theta^{(k)})\| < \epsilon.}
#' @slot tol_l Relative convergence tolerance denoted by \eqn{\delta}. The
#'   iteration stops when
#'   \deqn{\|\lambda^{(k)} - \lambda^{(k - 1)}\| <
#'   \delta\|\lambda^{(k - 1)}\| + \delta^2.}
#' @slot step Step size \eqn{\gamma} for the projected gradient descent
#'   method.
#' @slot th Threshold for the negative empirical log-likelihood ratio value.
#'   The iteration stops if the value exceeds the threshold. Defaults to
#'   \code{NULL} and sets the threshold to \code{200 * d}, where \code{d}
#'   corresponds to the degrees of freedom of the limiting chi-squared
#'   distribution of the statistic.
#' @slot nthreads Number of threads for parallel computation via OpenMP (if
#'   available). Defaults to the half of the available threads. For better
#'   performance, it is recommended to limit the number of threads to the
#'   number of physical cores. Note that it only applies to the following
#'   functions that involve multiple evaluations or minimizations:
#'   \itemize{
#'   \item{\code{\link{confreg}}}
#'   \item{\code{\link{el_lm}}}
#'   \item{\code{\link{el_glm}}}
#'   \item{\code{\link{eld}}}}
#' @examples
#' showClass("ControlEL")
setClass("ControlEL",
  slots = c(
    maxit = "integer", maxit_l = "integer", tol = "numeric", tol_l = "numeric",
    step = "ANY", th = "ANY", nthreads = "integer"
  ),
  prototype = list(
    maxit = 200L, maxit_l = 50L, tol = 1e-06, tol_l = 1e-06,
    step = NULL, th = NULL, nthreads = NULL
  )
)

#' S4 class \linkS4class{logLikEL}
#'
#' S4 class for empirical log-likelihood.
#'
#' @slot logLik Empirical log-likelihood.
#' @slot df Degrees of freedom or the number of (estimated) parameters in the
#'   model.
#' @examples
#' showClass("logLikEL")
setClass("logLikEL", slots = c(logLik = "numeric", df = "integer"))
