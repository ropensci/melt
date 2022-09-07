#' Empirical likelihood for the standard deviation
#'
#' Computes empirical likelihood for the standard deviation.
#'
#' @param x A numeric vector, or an object that can be coerced to a numeric
#'   vector.
#' @param mean A single numeric for the (known) mean value.
#' @param sd A positive single numeric for the parameter value to be tested.
#' @param weights An optional numeric vector of weights to be used in the
#'   fitting process. The length of the vector must be the same as the length of
#'   `x`. Defaults to `NULL`, corresponding to identical weights. If non-`NULL`,
#'   weighted empirical likelihood is computed.
#' @param control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @details Let \eqn{X_i} be independent and identically random variable from an
#'   unknown distribution \eqn{P} for \eqn{i = 1, \dots, n}. We assume that
#'   \eqn{{\textrm{E}[X_i]} = {\mu_0}} is known and that \eqn{P} has a variance
#'   \eqn{\sigma_0^2}. Given a value of \eqn{\sigma}, the
#'   (profile) empirical likelihood ratio is defined by
#'   \deqn{R(\sigma) =
#'   \max_{p_i}\left\{\prod_{i = 1}^n np_i :
#'   \sum_{i = 1}^n p_i (X_i - \mu_0)^2 = \sigma^2,\
#'   p_i \geq 0,\
#'   \sum_{i = 1}^n p_i = 1
#'   \right\}.}
#'   [el_sd()] computes the empirical log-likelihood ratio statistic
#'   \eqn{-2\log R(\sigma)}, along with other values in \linkS4class{SD}.
#' @return An object of class \linkS4class{SD}.
#' @seealso \linkS4class{EL}, \linkS4class{SD}, [el_mean()], [elt()],
#'   [el_control()]
#' @examples
#' data("women")
#' x <- women$height
#' w <- women$weight
#' el_sd(x, mean = 65, sd = 5, weights = w)
#' @export
el_sd <- function(x, mean, sd, weights = NULL, control = el_control()) {
  nm <- names(x)
  x <- as.vector(x, mode = "numeric")
  stopifnot(
    "`x` must have at least two observations." = (length(x) >= 2L),
    "`x` must be a finite numeric vector" = (isTRUE(all(is.finite(x)))),
    "`mean` must be a finite single numeric." =
      (isTRUE(is.numeric(mean) && length(mean) == 1L && is.finite(mean))),
    "`sd` must be a finite single numeric." =
      (isTRUE(is.numeric(sd) && length(sd) == 1L && is.finite(sd))),
    "`sd` must be a positive single numeric." = (sd >= .Machine$double.eps),
    "Invalid `control` specified." = (is(control, "ControlEL"))
  )
  n <- length(x)
  mm <- (x - mean)^2L
  if (length(nm) == n) {
    names(mm) <- nm
  }
  w <- validate_weights(weights, n)
  if (length(w) == 0L) {
    est <- sqrt(sum((x - mean)^2L) / n)
  } else {
    est <- sqrt(sum((x - mean)^2L * w) / n)
    names(w) <- nm
  }
  out <- compute_EL("sd", sd, mm, control@maxit_l, control@tol_l, control@th, w)
  optim <- validate_optim(out$optim)
  names(optim$sd) <- names(sd)
  if (control@verbose) {
    message(
      "Convergence ",
      if (out$optim$convergence) "achieved." else "failed."
    )
  }
  new("SD",
    optim = optim, logp = setNames(out$logp, names(mm)), logl = out$logl,
    loglr = out$loglr, statistic = out$statistic, df = 1L,
    pval = pchisq(out$statistic, df = 1L, lower.tail = FALSE), nobs = n,
    npar = 1L, weights = w, coefficients = est, method = "sd",
    data = if (control@keep_data) mm else NULL
  )
}
