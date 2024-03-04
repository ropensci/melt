#' Empirical likelihood for the mean
#'
#' Computes empirical likelihood for the mean.
#'
#' @param x A numeric matrix, or an object that can be coerced to a numeric
#'   matrix. Each row corresponds to an observation. The number of rows must be
#'   greater than the number of columns.
#' @param par A numeric vector of parameter values to be tested. The length of
#'   the vector must be the same as the number of columns in `x`.
#' @param weights An optional numeric vector of weights to be used in the
#'   fitting process. The length of the vector must be the same as the number of
#'   rows in `x`. Defaults to `NULL`, corresponding to identical weights.
#'   If non-`NULL`, weighted empirical likelihood is computed.
#' @param control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @details Let \eqn{X_i} be independent and identically distributed
#'   \eqn{p}-dimensional random variable from an unknown distribution \eqn{P}
#'   for \eqn{i = 1, \dots, n}. We assume that \eqn{{\textrm{E}[X_i]} =
#'   {\theta_0} \in {\rm{I\!R}}^p} and that \eqn{P} has a positive definite
#'   covariance matrix. Given a value of \eqn{\theta}, the (profile) empirical
#'   likelihood ratio is defined by
#'   \deqn{R(\theta) =
#'   \max_{p_i}\left\{\prod_{i = 1}^n np_i :
#'   \sum_{i = 1}^n p_i X_i = \theta,\
#'   p_i \geq 0,\
#'   \sum_{i = 1}^n p_i = 1
#'   \right\}.}
#'   [el_mean()] computes the empirical log-likelihood ratio statistic
#'   \eqn{-2\log R(\theta)}, along with other values in \linkS4class{EL}.
#' @return An object of class \linkS4class{EL}.
#' @references
#'   Owen A (1990).
#'   ``Empirical Likelihood Ratio Confidence Regions.''
#'   \emph{The Annals of Statistics}, **18**(1), 90--120.
#'   \doi{10.1214/aos/1176347494}.
#' @seealso \linkS4class{EL}, [elt()], [el_eval()], [el_control()]
#' @examples
#' ## Scalar mean
#' data("precip")
#' fit <- el_mean(precip, 30)
#' fit
#' summary(fit)
#'
#' ## Vector mean
#' data("faithful")
#' fit2 <- el_mean(faithful, par = c(3.5, 70))
#' summary(fit2)
#'
#' ## Weighted data
#' w <- rep(c(1, 2), each = nrow(faithful) / 2)
#' fit3 <- el_mean(faithful, par = c(3.5, 70), weights = w)
#' summary(fit3)
#' @export
el_mean <- function(x,
                    par,
                    weights = NULL,
                    control = el_control()) {
  assert_class(control, "ControlEL")
  x <- as.matrix(x, rownames.force = TRUE)
  assert_matrix(x,
    mode = "numeric", any.missing = FALSE, all.missing = FALSE, min.rows = 2L,
    min.cols = 1L
  )
  assert_numeric(x, finite = TRUE)
  stopifnot(
    "`x` must must have larger number of rows than columns." =
      nrow(x) > ncol(x),
    "`x` must have full column rank." = get_rank(x) == ncol(x)
  )
  stopifnot(
    "`par` must be a finite numeric vector." =
      (isTRUE(is.numeric(par) && all(is.finite(par))))
  )
  n <- nrow(x)
  p <- ncol(x)
  if (length(par) != p) {
    stop(gettextf("Length of `par` must be %d.", p, domain = NA))
  }
  w <- validate_weights(weights, n)
  if (length(w) == 0L) {
    est <- colMeans(x)
  } else {
    est <- colSums(x * w) / n
    names(w) <- rownames(x)
  }
  out <- compute_EL(
    "mean", par, x, control@maxit_l, control@tol_l, control@th, w
  )
  optim <- validate_optim(out$optim)
  names(optim$par) <- names(est)
  optim$cstr <- FALSE
  if (control@verbose) {
    message(
      "Convergence ",
      if (out$optim$convergence) "achieved." else "failed."
    )
  }
  new("EL",
    optim = optim, logp = setNames(out$logp, rownames(x)), logl = out$logl,
    loglr = out$loglr, statistic = out$statistic, df = p,
    pval = pchisq(out$statistic, df = p, lower.tail = FALSE), nobs = n,
    npar = p, weights = w, coefficients = est, method = "mean",
    data = if (control@keep_data) x else NULL, control = control
  )
}
