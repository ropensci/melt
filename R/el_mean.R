#' Empirical likelihood for the mean
#'
#' Computes empirical likelihood for the mean.
#'
#' @param par A numeric vector of parameter values to be tested.
#' @param x A numeric matrix, or an object that can be coerced to a numeric
#'   matrix. Each row corresponds to an observation.
#' @param weights An optional numeric vector of weights to be used in the
#'   fitting process. Defaults to \code{NULL}, corresponding to identical
#'   weights. If non-\code{NULL}, weighted empirical likelihood is computed.
#' @param control A list of control parameters set by \code{\link{el_control}}.
#' @param model A logical. If \code{TRUE} the data matrix used for model fitting
#'   is returned.
#' @return An object of class \linkS4class{EL}.
#' @references Owen, Art. 1990. “Empirical Likelihood Ratio Confidence Regions.”
#'   The Annals of Statistics 18 (1): 90–120. \doi{10.1214/aos/1176347494}.
#' @seealso \link{el_control}, \link{el_eval}, \link{lht}
#' @examples
#' # scalar mean
#' par <- 0
#' x <- rnorm(100)
#' el_mean(par, x)
#'
#' # vector mean
#' par <- c(0, 0)
#' x <- matrix(rnorm(100), ncol = 2)
#' el_mean(par, x)
#'
#' # weighted data
#' par <- c(0, 0)
#' x <- matrix(rnorm(100), ncol = 2)
#' w <- rep(c(1, 2), each = 25)
#' el_mean(par, x, w)
#' @importFrom methods is new
#' @importFrom stats pchisq
#' @export
el_mean <- function(par, x, weights = NULL, control = el_control(),
                    model = TRUE) {
  mm <- as.matrix(x)
  n <- nrow(mm)
  p <- ncol(mm)
  if (!is.numeric(mm) || !all(is.finite(mm))) {
    stop("'x' must be a finite numeric matrix")
  }
  if (n < 2L) {
    stop("not enough 'x' observations")
  }
  if (get_rank_(mm) != p || n <= p) {
    stop("'x' must have full column rank")
  }
  if (!is.numeric(par) || !all(is.finite(par))) {
    stop("'par' must be a finite numeric vector")
  }
  if (length(par) != p) {
    stop("'par' and 'x' have incompatible dimensions")
  }
  w <- check_weights(weights, n)
  if (!is.null(weights)) {
    est <- colSums(mm * w) / n
  } else {
    est <- colMeans(mm)
  }
  if (!is(control, "ControlEL")) {
    stop("invalid 'control' specified")
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
