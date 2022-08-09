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
#' @details `el_eval` evaluates empirical likelihood with a \eqn{n} by
#'   \eqn{p} numeric matrix argument `g`, whose \eqn{i}th row is
#'   \eqn{g(X_i, \theta)}. Since the estimating function can be arbitrary,
#'   `el_eval` does not return an object of class \linkS4class{EL}, and the
#'   associated generics and methods are not applicable.
#' @return A list with the following components:
#'   * `optim` A list with the following optimization results:
#'     * `lambda` Lagrange multiplier of the dual problem.
#'     * `iterations` Number of iterations performed.
#'     * `convergence` Convergence status.
#'   * `logp` Log probabilities obtained from empirical likelihood.
#'   * `logl` Empirical log-likelihood.
#'   * `loglr` Empirical log-likelihood ratio.
#'   * `statistic` Minus twice the empirical log-likelihood ratio statistic
#'   that has an asymptotic chi-square distribution.
#'   * `df` Degrees of freedom of the statistic.
#'   * `pval` \eqn{p}-value of the statistic.
#'   * `npar` Number of parameters.
#'   * `weights` Re-scaled weights used for model fitting.
#' @references Qin J, Lawless J (1994).
#'   “Empirical Likelihood and General Estimating Equations.”
#'   The Annals of Statistics, 22(1), 300–325. \doi{10.1214/aos/1176325370}.
#' @seealso [el_control()]
#' @examples
#' set.seed(3271)
#' x <- rnorm(50)
#' par <- 0
#' g <- x - par
#' el_eval(g, weights = rep(c(1, 2), each = 25))
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
