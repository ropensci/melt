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
#' @return An object of class \linkS4class{EL}.
#' @references Owen A (1990).
#'   “Empirical Likelihood Ratio Confidence Regions.”
#'   The Annals of Statistics, 18(1), 90–120.
#'   \doi{10.1214/aos/1176347494}.
#' @seealso [el_control()], [el_eval()], [elt()]
#' @examples
#' ## Scalar mean
#' set.seed(414)
#' x <- rnorm(100)
#' par <- 0
#' el_mean(x, par)
#'
#' ## Vector mean
#' x <- matrix(rnorm(100), ncol = 2)
#' par <- c(0, 0)
#' el_mean(x, par)
#'
#' ## Weighted data
#' x <- matrix(rnorm(100), ncol = 2)
#' par <- c(0, 0)
#' w <- rep(c(1, 2), each = 25)
#' el_mean(x, par, w)
#' @importFrom methods is new
#' @importFrom stats pchisq setNames
#' @export
#' @srrstats {G2.0, G2.0a} Assertions on lengths of inputs are clarified
#'   throughout the package documentation.
#' @srrstats {G2.16} All functions in the package strictly prohibit undefined
#'   values. They will trigger error messages in all cases.
#' @srrstats {G2.7} `el_mean()` accepts a numeric matrix (or an object that can
#'   be coerced to a numeric matrix by `as.matrix()`) for the argument `x`. This
#'   includes a data frame with only numeric variables. Unit tests use a data
#'   frame for the argument `x`. See `tests/testthat/test-confint.R`.
el_mean <- function(x,
                    par,
                    weights = NULL,
                    control = el_control()) {
  mm <- validate_x(x)
  stopifnot(
    "`par` must be a finite numeric vector." =
      (isTRUE(is.numeric(par) && all(is.finite(par)))),
    "Invalid `control` specified." = (is(control, "ControlEL"))
  )
  n <- nrow(mm)
  p <- ncol(mm)
  if (length(par) != p) {
    stop(gettextf("Length of `par` must be %d.", p, domain = NA))
  }
  w <- validate_weights(weights, n)
  if (length(w) == 0L) {
    est <- colMeans(mm)
  } else {
    est <- colSums(mm * w) / n
    names(w) <- rownames(mm)
  }
  out <- compute_EL(
    "mean", par, mm, control@maxit_l, control@tol_l, control@th, w
  )
  optim <- validate_optim(out$optim)
  names(optim$par) <- names(est)
  if (control@verbose) {
    message(
      "Convergence ",
      if (out$optim$convergence) "achieved." else "failed."
    )
  }
  new("EL",
    optim = optim, logp = setNames(out$logp, rownames(mm)), logl = out$logl,
    loglr = out$loglr, statistic = out$statistic, df = p,
    pval = pchisq(out$statistic, df = p, lower.tail = FALSE), nobs = n,
    npar = p, weights = w, data = if (control@keep_data) mm else NULL,
    coefficients = est, method = "mean"
  )
}
