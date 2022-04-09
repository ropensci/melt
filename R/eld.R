#' Empirical likelihood displacement
#'
#' Computes empirical likelihood displacement for model diagnostics and outlier
#'   detection.
#'
#' @param object An object of class \code{`el`}.
#' @param control A list of control parameters. See ‘Details’ in
#'   \code{\link{el_eval}}.
#' @details Let \eqn{l(\theta)} be the empirical log-likelihood function based
#'   on the full sample with \eqn{n} observations. The maximum empirical
#'   likelihood estimate is denoted by \eqn{\hat{\theta}}. Consider a reduced
#'   sample with the \eqn{i}th observation deleted and the corresponding
#'   estimate \eqn{\hat{\theta}_{(i)}}. The empirical likelihood displacement is
#'   defined by
#'   \deqn{\textnormal{ELD}_i = 2\{l(\hat{\theta}) - l(\hat{\theta}_{(i)})\}.}
#'   If the value of \eqn{\textnormal{ELD}_i } is large, then the \eqn{i}th
#'   observation is an influential point and can be inspected as a possible
#'   outlier. \code{eld} computes \eqn{\textnormal{ELD}_i } for
#'   \eqn{i = 1, \dots, n }.
#' @return A numeric vector of class \code{"el"}.
#' @references Lazar, Nicole A. 2005. “Assessing the Effect of Individual Data
#'   Points on Inference From Empirical Likelihood.” Journal of Computational
#'   and Graphical Statistics 14 (3): 626–42.
#'   \doi{10.1198/106186005X59568}.
#' @references Zhu, H., J. G. Ibrahim, N. Tang, and H. Zhang. 2008. “Diagnostic
#'   Measures for Empirical Likelihood of General Estimating Equations.”
#'   Biometrika 95 (2): 489–507.
#'   \doi{10.1093/biomet/asm094}.
#' @seealso \link{plot.eld}
#' @examples
#' x <- rnorm(10)
#' y <- 10
#' fit <- el_mean(0, c(x, y))
#' eld(fit)
#' @export
eld <- function(object, control = list()) {
  if (!inherits(object, "el"))
    stop("invalid 'object' supplied")
  if (inherits(object, "elt"))
    stop("method not applicable for 'elt' object")
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
