#' S4 class \linkS4class{EL}
#'
#' S4 class for empirical likelihood.
#'
#' @slot optim List with the following optimization results:
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
)

#' S4 class \linkS4class{MinEL}
#'
#' S4 class for constrained empirical likelihood.
#'
#' @slot optim List with the following optimization results:
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
#' S4 class for empirical likelihood with linear models.
#'
#' @slot parTests List with the test results for each parameter:
#'   \itemize{
#'   \item{\code{statistic } }{Numeric vector of chi-squared statistics.}
#'   \item{\code{convergence } }{Logical vector. \code{TRUE} indicates
#'   convergence of the algorithm.}
#'   }
#' @slot optim List with the following optimization results:
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
#' showClass("LM")
setClass("LM",
  contains = "MinEL",
  slots = c(
    parTests = "list", na.action = "ANY", xlevels = "ANY", call = "ANY",
    terms = "ANY"
  )
)

#' S4 class \linkS4class{ConfregEL}
#'
#' S4 class for confidence region.
#'
#' @slot points Numeric matrix with two columns for boundary points of a
#'   confidence region.
#' @slot estimates Numeric vector of length two for parameter estimates.
#' @slot level Confidence level required.
#' @slot cv Critical value for calibration of empirical likelihood ratio
#'   statistic.
#' @slot pnames Character vector of length two for the name of parameters.
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
#' @slot eld Numeric vector of empirical likelihood displacement values.
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
