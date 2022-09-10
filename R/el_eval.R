#' Empirical likelihood for general estimating functions
#'
#' Computes empirical likelihood with general estimating functions.
#'
#' @param g A numeric matrix, or an object that can be coerced to a numeric
#'   matrix. Each row corresponds to an observation of an estimating function.
#'   The number of rows must be greater than the number of columns.
#' @param weights An optional numeric vector of weights to be used in the
#'   fitting process. The length of the vector must be the same as the number of
#'   rows in `g`. Defaults to `NULL`, corresponding to identical weights. If
#'   non-`NULL`, weighted empirical likelihood is computed.
#' @param control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @details Let \eqn{X_i} be independent and identically distributed
#'   \eqn{p}-dimensional random variable from an unknown distribution \eqn{P}
#'   for \eqn{i = 1, \dots, n}. We assume that \eqn{P} has a positive definite
#'   covariance matrix. For a parameter of interest
#'   \eqn{\theta(F) \in {\rm{I\!R}}^p}, consider a \eqn{p}-dimensional smooth
#'   estimating function \eqn{g(X_i, \theta)} with a moment condition
#'   \deqn{\textrm{E}[g(X_i, \theta)] = 0.}
#'   We assume that there exists an unique \eqn{\theta_0} that solves the above
#'   equation. Given a value of \eqn{\theta}, the (profile) empirical likelihood
#'   ratio is defined by
#'   \deqn{R(\theta) =
#'   \max_{p_i}\left\{\prod_{i = 1}^n np_i :
#'   \sum_{i = 1}^n p_i g(X_i, \theta) = 0, p_i \geq 0, \sum_{i = 1}^n p_i = 1
#'   \right\}.}
#'   [el_mean()] computes the empirical log-likelihood ratio statistic
#'   \eqn{-2\log R(\theta)} with the \eqn{n} by \eqn{p} numeric matrix `g`,
#'   whose \eqn{i}th row is \eqn{g(X_i, \theta)}. Since the estimating function
#'   can be arbitrary, [el_eval()] does not return an object of class
#'   \linkS4class{EL}, and the associated generics and methods are not
#'   applicable.
#' @return A list of the following optimization results:
#'   * `optim` A list with the following optimization results:
#'     * `lambda` A numeric vector of the Lagrange multipliers of the dual
#'     problem.
#'     * `iterations` A single integer for the number of iterations performed.
#'     * `convergence` A single logical for the convergence status.
#'   * `logp` A numeric vector of the log probabilities of the empirical
#'   likelihood.
#'   * `logl` A single numeric of the empirical log-likelihood.
#'   * `loglr` A single numeric of the empirical log-likelihood ratio.
#'   * `statistic` A single numeric of minus twice the empirical log-likelihood
#'   ratio with an asymptotic chi-square distribution.
#'   * `df` A single integer for the degrees of freedom of the statistic.
#'   * `pval` A single numeric for the \eqn{p}-value of the statistic.
#'   * `nobs` A single integer for the number of observations.
#'   * `npar` A single integer for the number of parameters.
#'   * `weights` A numeric vector of the re-scaled weights used for the model
#'   fitting.
#' @references Qin J, Lawless J (1994).
#'   “Empirical Likelihood and General Estimating Equations.”
#'   \emph{The Annals of Statistics}, 22(1), 300--325.
#'   \doi{10.1214/aos/1176325370}.
#' @seealso \linkS4class{EL}, [el_control()]
#' @examples
#' set.seed(123526)
#' mu <- 0
#' sigma <- 1
#' x <- rnorm(100)
#' g <- matrix(c(x - mu, (x - mu)^2 - sigma^2), ncol = 2)
#' el_eval(g, weights = rep(c(1, 2), each = 50))
#' @export
el_eval <- function(g, weights = NULL, control = el_control()) {
  mm <- as.matrix(g, rownames.force = TRUE)
  nm <- rownames(mm)
  n <- nrow(mm)
  p <- ncol(mm)
  stopifnot(
    "`g` must have at least two observations." = (n >= 2L),
    "`g` must be a finite numeric matrix." =
      (isTRUE(is.numeric(mm) && all(is.finite(mm)))),
    "`g` must have full column rank." = (isTRUE(n > p && get_rank(mm) == p)),
    "Invalid `control` specified." = (is(control, "ControlEL"))
  )
  w <- validate_weights(weights, n)
  names(w) <- if (length(w) != 0L) nm else NULL
  out <- compute_generic_EL(mm, control@maxit_l, control@tol_l, control@th, w)
  names(out$logp) <- nm
  out$df <- p
  out$pval <- pchisq(q = out$statistic, df = out$df, lower.tail = FALSE)
  out$nobs <- n
  out$npar <- p
  out$weights <- w
  out
}
