#' Empirical likelihood for general estimating functions
#'
#' Computes empirical likelihood with general estimating functions.
#'
#' @param g A numeric matrix, or an object that can be coerced to a numeric
#'   matrix. Each row corresponds to an observation of an estimating function.
#' @param weights An optional numeric vector of weights to be used in the
#'   fitting process. Defaults to \code{NULL}, corresponding to identical
#'   weights. If non-\code{NULL}, weighted empirical likelihood is computed.
#' @param control An object of class \linkS4class{ControlEL} constructed by
#'   \code{\link{el_control}}.
#' @details \code{el_eval} evaluates empirical likelihood with a \eqn{n} by
#'   \eqn{p} numeric matrix argument \code{g}, whose \eqn{i}th row is
#'   \eqn{g(X_i, \theta)}. Since the estimating function can be arbitrary,
#'   \code{el_eval} does not return an object of class \linkS4class{EL}, and the
#'   associated generics and methods are not applicable.
#' @return A list with the following components:
#'   \item{optim}{A list with the following optimization results:
#'     \itemize{
#'     \item{\code{lambda } }{Lagrange multiplier of the dual problem.}
#'     \item{\code{iterations } }{Number of iterations performed.}
#'     \item{\code{convergence } }{Convergence status.}
#'     }
#'   }
#'   \item{logp}{Log probabilities obtained from empirical likelihood.}
#'   \item{logl}{Empirical log-likelihood.}
#'   \item{loglr}{Empirical log-likelihood ratio.}
#'   \item{statistic}{Minus twice the empirical log-likelihood ratio statistic
#'   that has an asymptotic chi-square distribution.}
#'   \item{df}{Degrees of freedom of the statistic.}
#'   \item{pval}{\eqn{p}-value of the statistic.}
#'   \item{npar}{Number of parameters.}
#'   \item{weights}{Rescaled weights used for model fitting.}
#' @references Qin, Jin, and Jerry Lawless. 1994.
#'   “Empirical Likelihood and General Estimating Equations.”
#'   The Annals of Statistics 22 (1): 300–325. \doi{10.1214/aos/1176325370}.
#' @seealso \link{el_control}
#' @examples
#' x <- rnorm(50)
#' par <- 0
#' g <- x - par
#' el_eval(g, weights = rep(c(1, 2), each = 25))
#' @importFrom stats pchisq
#' @export
el_eval <- function(g, weights = NULL, control = el_control()) {
  mm <- as.matrix(g)
  n <- nrow(mm)
  p <- ncol(mm)
  stopifnot(
    "not enough observations in 'g' " = (n >= 2L),
    "'g' must have full column rank" = (n > p),
    "'g' must be a numeric matrix" = (is.numeric(mm)),
    "'g' must be a finite numeric matrix" = (all(is.finite(mm))),
    "'g' must have full column rank" = (get_rank_(mm) == p),
    "invalid 'control' specified" = (is(control, "ControlEL"))
  )
  w <- check_weights(weights, n)
  out <- eval_g_(mm, control@maxit_l, control@tol_l, control@th, w)
  out$df <- p
  out$p.value <- pchisq(q = out$statistic, df = out$df, lower.tail = FALSE)
  out$npar <- p
  out$weights <- w
  out
}
