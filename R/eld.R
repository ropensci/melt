#' Empirical likelihood displacement
#'
#' Empirical likelihood displacement for objects inheriting from class
#'   \code{"el_test"}.
#'
#' @param object A fitted \code{"el_test"} object.
#' @param control A list of control parameters. See ‘Details’.
#' @export
eld <- function(object, control = list()) {
  if (!inherits(object, "el_test"))
    stop("invalid 'object' supplied")
  if (is.null(object$data.matrix))
    stop("'object' has no 'data.matrix'; fit the model with 'model' = TRUE")
  optcfg <- check_control(control)
  maxit <- optcfg$maxit
  tol <- optcfg$tol
  th <- optcfg$th
  wt <- object$weights
  out <- eld_(object$optim$method, object$coefficients, object$data.matrix,
              maxit, tol, th , wt)
  setNames(out, "eld")
  class(out) <- "eld"
  out
}

#' Plot for eld objects
#'
#' Takes a fitted \code{`eld`} object and plots the empirical likelihood
#'   displacement values versus the observation index.
#'
#' @param x An object of class \code{`eld`}.
#' @param ... Additional arguments to be passed to \code{\link[base]{plot}}.
#' @param main A title for the plot.
#' @param ylab A label for the y-axis.
#' @param pch A vector of plotting characters or symbols to use.
#' @seealso \link{eld}
#' @examples
#' data("clothianidin")
#' fit <- el_lm(clo ~ trt, clothianidin)
#' eld <- eld(fit)
#' plot(eld)
#' @export
plot.eld <- function(x, ..., main = "Empirical Likelihood Displacement",
                     ylab = "ELD", pch = 21) {
  plot(x[], ..., main = main, ylab = ylab, pch = pch)
}
