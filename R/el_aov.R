#' Fit an one-way analysis of variance model with empirical likelihood
#'
#' \emph{This function is deprecated in favor of \link{el_lm} and will be
#'   removed in a future release.}
#'
#' @param formula A formula object. It must specify variables for response and
#'   treatment as 'response ~ treatment'.
#' @param data A data frame containing the variables in the formula.
#' @param maxit Maximum number of iterations for optimization. Defaults to
#'   10000.
#' @param abstol Absolute convergence tolerance for optimization. Defaults to
#'   1e-08.
#' @return A list with class \code{"el_aov"}.
#' @references Owen, Art. 1991. “Empirical Likelihood for Linear Models.”
#'   The Annals of Statistics 19 (4).
#'   \doi{10.1214/aos/1176348368}.
#' @seealso \link{el_lm}, \link{el_test}
#' @examples
#' \dontrun{
#' data("clothianidin")
#' el_aov(clo ~ trt, clothianidin)}
#' @importFrom stats .getXlevels setNames terms
#' @export
el_aov <- function(formula, data, maxit = 1e04, abstol = 1e-8) {
  .Deprecated("el_lm")
  ## check formula
  f <- attributes(terms(formula))
  if (any(
    # response required & no arbitrary manipulation on intercept
    f$response == 0, f$intercept == 0,
    length(f$variables) != 3,
    # no other formula
    typeof(f$variables[[3]]) != "symbol" ||
      length(f$variables[[3]]) != 1
  )) {
    stop("invalied model formula. specify formula as 'response ~ treatment'")
  }

  ## extract model frame
  mf <- cl <- match.call()
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(model.frame)
  mf <- eval(mf, parent.frame())

  ## type conversion
  # response
  mf[[1L]] <- as.numeric(mf[[1L]])
  # treatment
  mf[[2L]] <- as.factor(mf[[2L]])

  ## extract model terms & levels
  # model terms
  mt <- attr(mf, "terms")
  # number of levels
  p <- nlevels(mf[[2L]])
  if (p < 2L) {
    stop("contrasts can be applied only to factors with 2 or more levels")
  }
  # levels
  lv <- .getXlevels(mt, mf)
  # name for coefficients
  nm <- paste0(names(lv), lv[[1L]])

  ## construct a general block design
  # incidence matrix
  c <- unclass(table(
    factor(row.names(mf), levels = unique(row.names(mf))),
    mf[[2L]]
  ))
  # model matrix
  x <- mf[[1L]] * c

  ## specify hypothesis
  lhs <- matrix(0, nrow = p - 1, ncol = p)
  if (p == 2L) {
    lhs <- matrix(c(1, -1), nrow = 1L)
  } else {
    diag(lhs) <- 1
    diag(lhs[, -1L]) <- -1
  }
  rhs <- rep(0, p - 1)

  ## test hypothesis
  out <- ELtest(x, c, lhs, rhs, threshold = 500, maxit, abstol)
  out$coefficients <- setNames(out$coefficients, nm)
  out$xlevels <- lv
  out$call <- cl
  out$terms <- mt
  out$model <- list(model.matrix = x, incidence.matrix = c)
  if (!out$optim$convergence) {
    warning("convergence failed\n")
  }
  class(out) <- "el_aov"
  out
}

#' @noRd
#' @export
print.el_aov <- function(x, ...) {
  cat("\nCall:\n")
  dput(x$call, control = NULL)
  cat("\nminimizer:\n")
  cat(format(round(x$optim$par, 4), scientific = FALSE))
  cat("\n\n")
  cat("statistic:\n")
  cat(format(round(x$optim$n2logLR, 4), scientific = FALSE))
  cat("\n\n")
}
