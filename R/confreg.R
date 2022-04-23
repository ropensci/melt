#' Confidence region for model parameters
#'
#' Computes boundary points of a two-dimensional confidence region for model
#'   parameters.
#'
#' @param object A fitted \code{"el"} object.
#' @param parm A specification of which parameters are to be given a confidence
#'   region, either a vector of numbers or a vector of names. It should be a
#'   vector of length two of the form \code{c(x, y)}. If missing, the first two
#'   parameter in \code{object} are considered.
#' @param level A confidence level required. Defaults to \code{0.95}.
#' @param cv A critical value for calibration of empirical likelihood ratio
#'   statistic. Defaults to \code{qchisq(level, 2L)}.
#' @param npoints Number of boundary points to compute. Defaults to \code{50}.
#' @param control A list of control parameters set by \code{\link{control_el}}.
#' @importFrom stats qchisq
#' @return An S4 object of class \code{"\link{ConfregEL}"}.
#' @references Owen, Art. 1990. “Empirical Likelihood Ratio Confidence Regions.”
#'   The Annals of Statistics 18 (1): 90–120. \doi{10.1214/aos/1176347494}.
#' @seealso \link{control_el}, \link{confint.el}, \link{lht}
#' @examples
#' par <- c(0, 0, 0)
#' x <- matrix(rnorm(90), ncol = 3)
#' fit <- el_mean(par, x)
#' confreg(fit, c(1, 3))
#' @importFrom methods new
#' @export
confreg <- function(object, parm, level = 0.95, cv = qchisq(level, 2L),
                    npoints = 50L, control = control_el()) {
  if (!inherits(object, "el")) {
    stop("invalid 'object' supplied")
  }
  est <- object$coefficients
  if (length(est) == 0L) {
    stop("'confreg' method is not applicable to an empty model")
  } else if (length(est) == 1L) {
    stop("'confreg' method is not applicable to a model with one parameter")
  }
  pnames <- if (is.null(names(est))) seq(length(est)) else names(est)
  if (!missing(parm)) {
    if (length(parm) != 2L) {
      stop("length of 'parm' must be two")
    }
    if (is.numeric(parm) && all(is.finite(parm))) {
      est <- est[parm]
      pnames <- pnames[parm]
      if (length(unique(est)) != 2L) {
        stop("only one parameter specified by 'parm'")
      }
    } else if (is.character(parm)) {
      if (is.null(names(est))) {
        stop("'parm' is not recognized since 'object' has no named parameters")
      }
      idx <- match(parm, pnames)
      est <- est[idx]
      pnames <- pnames[idx]
      if (length(unique(est)) != 2L) {
        stop("only one parameter specified by 'parm'")
      }
    } else {
      stop("invalid 'parm' specified")
    }
  } else {
    est <- est[c(1L, 2L)]
    pnames <- pnames[c(1L, 2L)]
  }
  if (!missing(level) &&
    (length(level) != 1L || !is.finite(level) || level < 0 || level > 1)) {
    stop("'level' must be a number between 0 and 1")
  }
  if (isTRUE(all.equal(level, 0))) {
    return(new("ConfregEL",
      points = est, estimates = est, level = level, cv = cv, pnames = c("x", "y")
    ))
  } else if (isTRUE(all.equal(level, 1))) {
    stop("'level' must be a number between 0 and 1")
  }
  cv <- check_cv(cv, control$th)
  npoints <- as.integer(npoints)
  if (npoints <= 0) {
    stop("'npoints' must be a positive integer")
  }
  if (!inherits(control, "control_el") || !is.list(control)) {
    stop("invalid 'control' supplied")
  }
  w <- if (is.null(object$weights)) numeric(length = 0L) else object$weights
  ang <- seq(0, 2 * pi, length.out = npoints + 1L)[-(npoints + 1L)]
  circ <- rbind(cos(ang), sin(ang))
  cr <- confreg_(
    object$optim$method, est, object$data.matrix, cv, circ,
    control$maxit, control$maxit_l, control$tol, control$tol_l,
    control$step, control$th, control$nthreads, w
  )
  new("ConfregEL",
    points = t(circ) * cr + rep(est, each = ncol(circ)),
    estimates = est, level = level, cv = cv, pnames = as.character(pnames)
  )
}
