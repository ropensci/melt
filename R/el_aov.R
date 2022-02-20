#' Fit an Analysis of Variance Model with Empirical Likelihood
#'
#' Fit an one-way analysis of variance model with empirical likelihood.
#'
#' @param formula A formula object. It must specify variables for response and treatment as 'response ~ treatment'.
#' @param data A data frame containing the variables in the formula.
#' @param maxit Maximum number of iterations for optimization. Defaults to 10000.
#' @param abstol Absolute convergence tolerance for optimization. Defaults to 1e-08.
#' @return A list with class \code{c("el_aov", "melt")}.
#' @references Owen, Art. 1991. “Empirical Likelihood for Linear Models.” The Annals of Statistics 19 (4). \doi{10.1214/aos/1176348368}.
#' @seealso \link{el_test}
#' @examples
#' data("clothianidin")
#' el_aov(clo ~ trt, clothianidin)
#' @importFrom stats .getXlevels setNames terms
#' @export
el_aov <- function(formula, data, maxit = 1e04, abstol = 1e-8) {
  ## check formula
  f <- attributes(terms(formula))
  if (any(
    # response required & no arbitrary manipulation on intercept
    f$response == 0, f$intercept == 0,
    length(f$variables) != 3,
    # no other formula
    typeof(f$variables[[3]]) != "symbol" ||
      length(f$variables[[3]]) != 1
  )
  ) {
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
  class(out) <- c("el_aov", oldClass(out))
  if (!out$optim$convergence) {
    warning("convergence failed\n")
  }
  out
}
