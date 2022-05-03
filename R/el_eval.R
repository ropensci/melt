#' Empirical likelihood for general estimating functions
#'
#' Computes empirical likelihood with general estimating functions.
#'
#' @param g A numeric matrix, or an object that can be coerced to a numeric
#'   matrix. Each row corresponds to an observation of an estimating
#'   function.
#' @param weights An optional numeric vector of weights to be used in the
#'   fitting process. Defaults to \code{NULL}, corresponding to identical
#'   weights. If non-\code{NULL}, weighted empirical likelihood is computed.
#' @param control A list of control parameters set by \code{\link{el_control}}.
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
#'   \item{statistic}{Minus twice the empirical log-likelihood ratio statistic that
#'   has an asymptotic chi-square distribution.}
#'   \item{df}{Degrees of freedom of the statistic.}
#'   \item{pval}{\eqn{p}-value of the statistic.}
#'   \item{npar}{Number of parameters.}
#'   \item{weights}{Rescaled weights used for model fitting.}
#' @references Qin, Jin, and Jerry Lawless. 1994.
#'   “Empirical Likelihood and General Estimating Equations.”
#'   The Annals of Statistics 22 (1): 300–325. \doi{10.1214/aos/1176325370}.
#' @seealso \link{el_control}
#' @examples
#' # test for variance with known mean
#' x <- rnorm(100L)
#' sigma <- 1
#' g <- x^2 - sigma^2
#' el_eval(g)
#' @importFrom stats pchisq
#' @export
el_eval <- function(g, weights = NULL, control = el_control()) {
  mm <- as.matrix(g)
  n <- nrow(mm)
  p <- ncol(mm)
  if (!is.numeric(mm) || !all(is.finite(mm))) {
    stop("'g' must be a finite numeric matrix")
  }
  if (n < 2L) {
    stop("not enough 'g' observations")
  }
  if (get_rank_(mm) != p || n <= p) {
    stop("'g' must have full column rank")
  }
  if (!is(control, "ControlEL")) {
    stop("invalid 'control' specified")
  }
  w <- check_weights(weights, n)
  out <- eval_g_(mm, control@maxit_l, control@tol_l, control@th, w)
  out$df <- p
  out$p.value <- pchisq(out$statistic, df = out$df, lower.tail = FALSE)
  out$npar <- p
  # if (!is.null(weights)) {
  #   out$weights <- w
  # }
  out$weights <- w
  out
}
