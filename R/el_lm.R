#' Fit a Linear Model with Empirical Likelihood
#'
#' Fit a linear model with empirical likelihood.
#'
#' @param formula A formula object.
#' @param data A data frame containing the variables in the formula.
#' @param na.action A function which indicates what should happen when the data contain NAs.
#' @param model Logical. If TRUE the model frame is returned.
#' @param maxit Maximum number of iterations for optimization. Defaults to 10000.
#' @param abstol Absolute convergence tolerance for optimization. Defaults to 1e-08.
#' @return A list with class \code{c("el_lm", "melt")}.
#' @references Owen, Art. 1991. “Empirical Likelihood for Linear Models.” The Annals of Statistics 19 (4). \doi{10.1214/aos/1176348368}.
#' @seealso \link{el_aov}
#' @examples
#' model <- el_lm(clo ~ trt, clothianidin)
#' summary(model)
#' @importFrom stats .getXlevels is.empty.model model.matrix model.response setNames terms
#' @export
el_lm <- function(formula, data, na.action, model = TRUE, maxit = 1e04, abstol = 1e-08) {
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  if (is.empty.model(mt)) stop("empty model specified")

  y <- model.response(mf, "numeric")
  nm <- names(y)
  ny <- if (is.matrix(y)) nrow(y) else length(y)
  x <- model.matrix(mt, mf)

  out <- EL_lm(x, y, rep(0, ny), threshold = 500, maxit, abstol)
  out$coefficients <- setNames(out$coefficients, colnames(x))
  out$residuals <- setNames(out$residuals, nm)
  out$fitted.values <- setNames(out$fitted.values, nm)
  out$df.residual = nrow(x) - out$rank
  out$na.action <- attr(mf, "na.action")
  out$xlevels <- .getXlevels(mt, mf)
  out$call <- cl
  out$terms <- mt
  if (model)
    out$model <- mf
  return(out)
}
