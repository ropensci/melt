#' Confidence intervals for model parameters
#'
#' Computes confidence intervals for one or more parameters in a fitted model.
#' Package \strong{melt} adds a method for objects inheriting from class
#' \code{"el"}.
#'
#' @param object A fitted \code{"el"} object.
#' @param parm A specification of which parameters are to be given confidence
#'   intervals, either a vector of numbers or a vector of names. If missing, all
#'   parameters are considered.
#' @param level A confidence level required. Defaults to \code{0.95}.
#' @param cv A critical value for calibration of empirical likelihood ratio
#'   statistic. Defaults to \code{qchisq(level, 1L)}.
#' @param control A list of control parameters set by
#'   \code{\link{control_el}}.
#' @param ... Some methods for this generic function require extra arguments.
#'   None are used in this method.
#' @importFrom stats qchisq
#' @return A matrix with columns giving lower and upper confidence limits for
#'  each parameter. In contrast to other methods that rely on studentization,
#'  the lower and upper limits obtained from empirical likelihood do not
#'  correspond to the \code{(1 - level) / 2} and \code{1 - (1 - level) / 2} in
#'  \%, respectively.
#' @references Kim, E., MacEachern, S., and Peruggia, M., (2021),
#'   "Empirical Likelihood for the Analysis of Experimental Designs,"
#'   \href{https://arxiv.org/abs/2112.09206}{arxiv:2112.09206}.
#' @references Owen, Art. 1990. “Empirical Likelihood Ratio Confidence Regions.”
#'   The Annals of Statistics 18 (1): 90–120. \doi{10.1214/aos/1176347494}.
#' @seealso \link{confreg}, \link{control_el}, \link{lht}
#' @examples
#' fit <- el_lm(formula = mpg ~ wt, data = mtcars)
#' confint(fit)
#' @export
confint.el <- function(object, parm, level = 0.95, cv = qchisq(level, 1L),
                       control = control_el(), ...) {
  est <- object$coefficients
  # no confidence interval for empty model
  if (length(est) == 0L) {
    ci <- matrix(, nrow = 0L, ncol = 2L)
    colnames(ci) <- c("lower", "upper")
    return(ci)
  }
  # index for tracking the parameters
  idx <- seq(length(est))
  # rownames of the confidence interval matrix
  pnames <- if (is.null(names(est))) idx else names(est)
  # if parm is supplied, modify idx and pnames accordingly
  if (!missing(parm)) {
    if (is.numeric(parm) && all(is.finite(parm))) {
      pnames <- pnames[parm]
      # idx <- match(pnames, names(est))
      idx <- if (is.null(names(est))) {
        match(pnames, idx)
      } else {
        match(pnames, names(est))
      }
    } else if (is.character(parm)) {
      idx <- match(parm, pnames)
      pnames <- parm
    } else {
      stop("invalid 'parm' specified")
    }
  }
  # number of rows of the confidence interval matrix
  p <- length(idx)
  # check level
  if (!missing(level) &&
    (length(level) != 1L || !is.finite(level) || level < 0 || level > 1)) {
    stop("'level' must be a number between 0 and 1")
  }
  if (isTRUE(all.equal(level, 0))) {
    ci <- matrix(rep(est[idx], 2L), ncol = 2L)
    colnames(ci) <- c("lower", "upper")
    return(ci)
  } else if (isTRUE(all.equal(level, 1))) {
    ci <- matrix(NA, nrow = p, ncol = 2L)
    ci[which(!is.na(idx)), ] <- c(-Inf, Inf)
    colnames(ci) <- c("lower", "upper")
    return(ci)
  }

  # check control
  if (!inherits(control, "control_el") || !is.list(control)) {
    stop("invalid 'control' supplied")
  }
  method <- object$optim$method
  maxit <- control$maxit
  maxit_l <- control$maxit_l
  tol <- control$tol
  tol_l <- control$tol_l
  step <- control$step
  th <- control$th
  nthreads <- control$nthreads
  w <- if (is.null(object$weights)) numeric(length = 0L) else object$weights
  cv <- check_cv(cv, th)
  # compute the confidence interval matrix
  if (isTRUE(all.equal(level, 0))) {
    ci <- matrix(rep(est[idx], 2L), ncol = 2L)
  } else if (isTRUE(all.equal(level, 1))) {
    ci <- matrix(NA, nrow = p, ncol = 2L)
    ci[which(!is.na(idx)), ] <- c(-Inf, Inf)
  } else if (all(is.na(idx))) {
    ci <- matrix(NA, nrow = p, ncol = 2L)
  } else if (any(is.na(idx))) {
    idx_na <- which(is.na(idx))
    ci <- matrix(NA, nrow = p, ncol = 2L)
    ci[-idx_na, ] <- confint_(method, est, object$data.matrix, cv, idx[-idx_na],
                              maxit, tol, step, th, nthreads, w)
  } else {
    ci <- confint_(method, est, object$data.matrix, cv, idx, maxit, maxit_l, tol,
                   tol_l, step, th, nthreads, w)
  }
  dimnames(ci) <- list(pnames, c("lower", "upper"))
  ci
}

#' Empirical log-likelihood
#'
#' Extracts empirical log-likelihood of a model represented by
#'   \code{object} evaluated at the estimated coefficients. Package
#'   \strong{melt} adds a method for objects inheriting from class \code{"el"}.
#'
#' @param object A fitted \code{"el"} object.
#' @param ... Some methods for this generic function require extra arguments.
#'   None are used in this method.
#' @return An object of class \code{"logLik"} with an attribute \code{df} that
#'   gives the number of (estimated) parameters in the model.
#' @examples
#' fit <- el_lm(formula = mpg ~ wt, data = mtcars)
#' logLik(fit)
#' @export
logLik.el <- function(object, ...) {
  if (!missing(...)) {
    warning("extra arguments are not supported")
  }
  p <- object$npar
  rhs <- object$coefficients
  out <- lht(object, rhs = rhs)
  val <- out$loglik
  attr(val, "df") <- p
  class(val) <- "logLik"
  val
}

#' @noRd
#' @export
print.el <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nEmpirical Likelihood Test:", x$optim$method, "\n\n")
  out <- character()
  if (!is.null(x$statistic)) {
    out <- c(
      out, paste("Chisq:", format(x$statistic, digits = digits)),
      paste("df:", x$df),
      paste("p-value:", format.pval(x$p.value, digits = digits))
    )
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  if (!is.null(x$coefficients)) {
    cat("maximum EL estimates:\n")
    print(x$coefficients, digits = digits, ...)
  }
  cat("\n")
  invisible(x)
}

#' @noRd
#' @export
print.elt <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nEmpirical Likelihood Linear Hypothesis Test:", x$optim$method, "\n\n")
  out <- character()
  if (!is.null(x$statistic)) {
    out <- c(
      out, paste("Chisq:", format(x$statistic, digits = digits)),
      paste("df:", x$df),
      paste("p-value:", format.pval(x$p.value, digits = digits))
    )
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  if (!is.null(x$coefficients)) {
    cat("maximum EL estimates:\n")
    print(x$coefficients, digits = digits, ...)
  }
  cat("\n")
  invisible(x)
}
