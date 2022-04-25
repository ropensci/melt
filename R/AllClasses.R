#' S4 class \code{\link{EL}}
#'
#' S4 class for empirical likelihood.
#'
#' @aliases EL
#' @slot optim List with the following optimization results:
#'   \itemize{
#'   \item{\code{method } }{Character for method dispatch in internal
#'   functions.}
#'   \item{\code{par } }{Value at which empirical likelihood is evaluated.}
#'   \item{\code{lambda } }{Lagrange multiplier of the dual problem.}
#'   \item{\code{logLR } }{Empirical log-likelihood ratio.}
#'   \item{\code{iterations } }{Number of iterations performed.}
#'   \item{\code{convergence } }{Convergence status.}
#'   }
#' @slot logp Log probabilities obtained from empirical likelihood.
#' @slot logl Empirical log-likelihood.
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
    optim = "list", logp = "numeric", logl = "numeric", statistic = "numeric",
    df = "integer", pval = "numeric", npar = "integer", weights = "numeric",
    dataMatrix = "matrix", coefficients = "numeric"
  ),
  prototype = list(
    optim = list(), logp = numeric(), logl = numeric(), statistic = numeric(),
    df = 0L, pval = numeric(), npar = 0L, weights = numeric(),
    dataMatrix = matrix(NA_real_, nrow = 0L, ncol = 0L),
    coefficients = numeric()
  )
)

#' S4 class \code{ConfregEL}
#'
#' S4 class for confidence region.
#'
#' @aliases ConfregEL
#' @slot points Numeric matrix with two columns for boundary points of a
#'   confidence region.
#' @slot estimates Numeric vector of length two for parameter estimates.
#' @slot level Confidence level required.
#' @slot cv Critical value for calibration of empirical likelihood ratio
#'   statistic.
#' @slot pnames Character vector of length two for the name of parameters.
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

#' S4 class \code{ELD}
#'
#' S4 class for empirical likelihood displacement.
#'
#' @aliases ELD
#' @slot eld Numeric vector of empirical likelihood displacement values.
setClass("ELD",
  slots = c(eld = "numeric"),
  prototype = list(eld = NA_real_)
)


