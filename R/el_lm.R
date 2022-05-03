#' Empirical likelihood for linear models
#'
#' Fits a linear model with empirical likelihood.
#'
#' @param formula An object of class \code{"\link[stats]{formula}"} (or one that
#'   can be coerced to that class): a symbolic description of the model to be
#'   fitted.
#' @param data An optional data frame, list or environment (or object coercible
#'   by \code{\link[base]{as.data.frame}} to a data frame) containing the
#'   variables in the formula. If not found in data, the variables are taken
#'   from \code{environment(formula)}.
#' @param weights An optional numeric vector of weights to be used in the
#'   fitting process. Defaults to \code{NULL}, corresponding to identical
#'   weights. If non-\code{NULL}, weighted empirical likelihood is computed.
#' @param na.action A function which indicates what should happen when the data
#'   contain \code{NA}s. The default is set by the \code{na.action} setting of
#'   \code{\link[base]{options}}, and is \code{na.fail} if that is unset.
#' @param control A list of control parameters set by \code{\link{el_control}}.
#' @param model A logical. If \code{TRUE} the data matrix used for fitting is
#'   returned.
#' @param ... Additional arguments to be passed to the low level regression
#'   fitting functions. See ‘Details’.
#' @details Suppose that we observe \eqn{n} independent random variables
#'   \eqn{(X_i, Y_i)} from a common distribution, where \eqn{X_i} is the
#'   \eqn{p}-dimensional covariate (including the intercept if any) and
#'   \eqn{Y_i} is the response. We consider the following linear regression
#'   model:
#'   \deqn{Y_i = X_i^\top \beta + \epsilon_i,}
#'   where \eqn{\beta = (\beta_0, \dots, \beta_{p-1})} is an unknown
#'   \eqn{p}-dimensional parameter and the errors \eqn{\epsilon_i} are
#'   independent random variables that satisfy
#'   \eqn{\textnormal{E}(\epsilon_i | X_i)} = 0. We assume that the errors have
#'   finite conditional variance. Then the least square estimator of \eqn{\beta}
#'   solves the following estimating equation:
#'   \deqn{\sum_{i = 1}^n(Y_i - X_i^\top \beta)X_i = 0.}
#'   \code{\link{el_lm}} first computes the parameter estimates by calling
#'   \code{\link[stats]{lm.fit}} (with \code{...} if any) since the maximum
#'   empirical likelihood estimator is the same as the least square estimator in
#'   our model. Next, it performs hypothesis tests based on asymptotic
#'   chi-squared distribution of empirical likelihood ratio statistics. Included
#'   in the tests are the overall test with
#'   \deqn{H_0: \beta_1 = \beta_2 = \cdots = \beta_{p-1} = 0,}
#'   and the tests for each parameter with
#'   \deqn{H_{0j}: \beta_j = 0,\ j = 0, \dots, p-1.}
#'   The test results are returned as \code{optim} and \code{parTests},
#'   respectively.
#' @return An object of class of \linkS4class{LM}.
#' @references Owen, Art. 1991. “Empirical Likelihood for Linear Models.”
#'   The Annals of Statistics 19 (4): 1725–47. \doi{10.1214/aos/1176348368}.
#' @seealso \link{el_control}, \link{el_glm}, \link{lht}
#' @examples
#' fit <- el_lm(mpg ~ wt, mtcars)
#' summary(fit)
#' @importFrom stats .getXlevels is.empty.model lm.fit lm.wfit model.matrix
#'   model.response model.weights pchisq
#' @export
el_lm <- function(formula, data, weights = NULL, na.action,
                  control = el_control(), model = TRUE, ...) {
  cl <- match.call()
  if (missing(data)) {
    data <- environment(formula)
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  w <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w)) {
    stop("'weights' must be a numeric vector")
  }
  if (!is.null(w) && any(w < 0)) {
    stop("negative weights not allowed")
  }
  if (is.matrix(y)) {
    stop("'el_lm' does not support multiple responses")
  }
  if (is.empty.model(mt)) {
    x <- NULL
    mm <- cbind(y, x)
    return(new("LM",
      optim = list(
        method = "lm", par = numeric(), lambda = numeric(),
        iterations = integer(), convergence = logical()
      ),
      misc = list(
        call = cl, terms = mt, xlevels = .getXlevels(mt, mf),
        na.action = attr(mf, "na.action")
      )
    ))
  } else {
    x <- model.matrix(mt, mf, NULL)
    z <- if (is.null(w)) {
      lm.fit(x, y, offset = NULL, singular.ok = FALSE, ...)
    } else {
      lm.wfit(x, y, w, offset = NULL, singular.ok = FALSE, ...)
    }
  }
  intercept <- attr(mt, "intercept")
  mm <- cbind(y, x)
  p <- ncol(x)
  w <- check_weights(w, nrow(mm))
  if (!is(control, "ControlEL")) {
    stop("invalid 'control' specified")
  }
  el <- lm_(
    mm, z$coefficients, intercept, control@maxit, control@maxit_l,
    control@tol, control@tol_l, control@step, control@th,
    control@nthreads, w
  )
  df <- if (intercept && p > 1L) p - 1L else p
  pval <- pchisq(el$statistic, df = df, lower.tail = FALSE)
  new("LM",
    optim = el$optim, logp = el$logp, logl = el$logl, loglr = el$loglr,
    statistic = el$statistic, df = df, pval = pval, npar = p, weights = w,
    data = if (model) mm else matrix(NA_real_, nrow = 0L, ncol = 0L),
    coefficients = z$coefficients, parTests = el$parTests,
    misc = list(
      call = cl, terms = mt, xlevels = .getXlevels(mt, mf),
      na.action = attr(mf, "na.action")
    )
  )
}
