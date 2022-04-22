#' Confidence region for model parameters
#'
#' Computes confidence region for two model parameters in a fitted model.
#'
#' @param object A fitted \code{"el"} object.
#' @param parm A specification of which parameters are to be given confidence
#'   intervals, either a vector of numbers or a vector of names. If missing, all
#'   parameters are considered.
#' @param level A confidence level required.
#' @param cv A critical value for calibration of empirical likelihood ratio
#'   statistic. Defaults to \code{qchisq(level, 2L)}.
#' @param control A list of control parameters set by \code{\link{control_el}}.
#' @importFrom stats qchisq
#' @return A matrix with columns giving lower and upper confidence limits for
#'  each parameter. In contrast to other methods that rely on studentization,
#'  the lower and upper limits obtained from empirical likelihood do not
#'  correspond to the \code{(1 - level) / 2} and \code{1 - (1 - level) / 2} in
#'  \%, respectively.
#' @references Owen, Art. 1990. “Empirical Likelihood Ratio Confidence Regions.”
#'   The Annals of Statistics 18 (1): 90–120. \doi{10.1214/aos/1176347494}.
#' @seealso \link{control_el}, \link{confint.el}
#' @examples
#' fit <- el_lm(formula = mpg ~ wt, data = mtcars)
#' confint(fit)
#' @export
confreg <- function(object, parm, level = 0.95, cv = qchisq(level, 2L),
                    n = 50, control = control_el()) {
  est <- object$coefficients
  # no confidence region for empty model
  if (length(est) == 0L) {
    ci <- matrix(, nrow = 0L, ncol = 2L)
    colnames(ci) <- c("lower", "upper")
    return(ci)
  }
  # # index for tracking the parameters
  # idx <- seq(length(est))
  # # rownames of the confidence interval matrix
  # pnames <- if (is.null(names(est))) idx else names(est)
  # # if parm is supplied, modify idx and pnames accordingly
  # if (!missing(parm)) {
  #   if (is.numeric(parm) && all(is.finite(parm))) {
  #     pnames <- pnames[parm]
  #     # idx <- match(pnames, names(est))
  #     idx <- if (is.null(names(est))) {
  #       match(pnames, idx)
  #     } else {
  #       match(pnames, names(est))
  #     }
  #   } else if (is.character(parm)) {
  #     idx <- match(parm, pnames)
  #     pnames <- parm
  #   } else {
  #     stop("invalid 'parm' specified")
  #   }
  # }
  # # number of rows of the confidence interval matrix
  # p <- length(idx)
  # # check level
  # if (!missing(level) &&
  #     (length(level) != 1L || !is.finite(level) || level < 0 || level > 1)) {
  #   stop("'level' must be a single number between 0 and 1")
  # }
  # if (isTRUE(all.equal(level, 0))) {
  #   ci <- matrix(rep(est[idx], 2L), ncol = 2L)
  #   colnames(ci) <- c("lower", "upper")
  #   return(ci)
  # } else if (isTRUE(all.equal(level, 1))) {
  #   ci <- matrix(NA, nrow = p, ncol = 2L)
  #   ci[which(!is.na(idx)), ] <- c(-Inf, Inf)
  #   colnames(ci) <- c("lower", "upper")
  #   return(ci)
  # }

  # check control
  if (!inherits(control, "control_el") || !is.list(control)) {
    stop("invalid 'control' supplied")
  }
  w <- object$weights
  if (is.null(w)) {
    w <- numeric(length = 0L)
  }
  cv <- tryCatch(as.numeric(cv), warning = function(w) NA,
                 error = function(e) NA)
  if (any(length(cv) != 1L, is.na(cv), is.infinite(cv))) {
    stop("'cv' is not a number")
  }
  if (cv < .Machine$double.eps) {
    stop("'cv' is too small")
  }
  if (!is.null(control$th) && cv > 2 * control$th) {
    stop("'cv' is too large")
  }
  ang <- seq(0, 2 * pi, length.out = n + 1)[-(n + 1)]
  circ  <- rbind(cos(ang), sin(ang))
  cr <- confreg_(object$optim$method, est, object$data.matrix, cv, circ,
                 control$maxit, control$maxit_l, control$tol, control$tol_l,
                 control$step, control$th, control$nthreads, w)
  new("ConfregEL", points = t(circ) * cr + rep(est, each = ncol(circ)),
      estimates = est,
      level = level)
  # t(circ) * cr + rep(est, each = ncol(circ))
}
