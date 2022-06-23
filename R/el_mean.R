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
#'   \code{\link{el_control}}.
#' @param model A single logical. If `TRUE` the data matrix used for model
#'   fitting is returned.
#' @return An object of class \linkS4class{EL}.
#' @references Owen, Art. 1990. “Empirical Likelihood Ratio Confidence Regions.”
#'   The Annals of Statistics 18 (1): 90–120. \doi{10.1214/aos/1176347494}.
#' @seealso \link{el_control}, \link{el_eval}, \link{elt}
#' @examples
#' # scalar mean
#' x <- rnorm(100)
#' par <- 0
#' el_mean(x, par)
#'
#' # vector mean
#' x <- matrix(rnorm(100), ncol = 2)
#' par <- c(0, 0)
#' el_mean(x, par)
#'
#' # weighted data
#' x <- matrix(rnorm(100), ncol = 2)
#' par <- c(0, 0)
#' w <- rep(c(1, 2), each = 25)
#' el_mean(x, par, w)
#' @importFrom methods is new
#' @importFrom stats pchisq
#' @export
el_mean <- function(x,
                    par,
                    weights = NULL,
                    control = el_control(),
                    model = TRUE) {
  mm <- check_x_(x)
  n <- nrow(mm)
  p <- ncol(mm)
  stopifnot(
    "'par' must be a numeric vector" = (is.numeric(par)),
    "'par' must be a finite numeric vector" = (all(is.finite(par))),
    "'par' and 'x' have incompatible dimensions" = (length(par) == p),
    "invalid 'control' specified" = (is(control, "ControlEL"))
  )
  w <- check_weights_(weights, n)
  if (!is.null(weights)) {
    est <- colSums(mm * w) / n
  } else {
    est <- colMeans(mm)
  }
  el <- eval_("mean", par, mm, control@maxit_l, control@tol_l, control@th, w)
  new("EL",
    optim = el$optim, logp = el$logp, logl = el$logl, loglr = el$loglr,
    statistic = el$statistic, df = p,
    pval = pchisq(el$statistic, df = p, lower.tail = FALSE), npar = p,
    weights = w,
    data = if (model) mm else matrix(NA_real_, nrow = 0L, ncol = 0L),
    coefficients = est, method = "mean"
  )
}
