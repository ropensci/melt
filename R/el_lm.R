#' Empirical likelihood for linear models
#'
#' Fits a linear model with empirical likelihood.
#'
#' @param formula An object of class [`formula`] (or one that can be coerced to
#'   that class) for a symbolic description of the model to be fitted.
#' @param data An optional data frame, list or environment (or object coercible
#'   by [as.data.frame()] to a data frame) containing the variables in
#'   `formula`. If not found in data, the variables are taken from
#'   `environment(formula)`.
#' @param weights An optional numeric vector of weights to be used in the
#'   fitting process. Defaults to `NULL`, corresponding to identical weights. If
#'   non-`NULL`, weighted empirical likelihood is computed.
#' @param na.action A function which indicates what should happen when the data
#'   contain `NA`s. The default is set by the `na.action` setting of
#'   [`options`], and is `na.fail` if that is unset.
#' @param control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @param ... Additional arguments to be passed to the low level regression
#'   fitting functions. See ‘Details’.
#' @details Suppose that we observe \eqn{n} independent random variables
#'   \eqn{(X_i, Y_i)} from a common distribution, where \eqn{X_i} is the
#'   \eqn{p}-dimensional covariate (including the intercept if any) and
#'   \eqn{Y_i} is the response. We consider the following linear regression
#'   model:
#'   \deqn{Y_i = X_i^\top \theta + \epsilon_i,}
#'   where \eqn{\theta = (\theta_0, \dots, \theta_{p-1})} is an unknown
#'   \eqn{p}-dimensional parameter and the errors \eqn{\epsilon_i} are
#'   independent random variables that satisfy
#'   \eqn{\textnormal{E}(\epsilon_i | X_i)} = 0. We assume that the errors have
#'   finite conditional variance. Then the least square estimator of
#'   \eqn{\theta} solves the following estimating equation:
#'   \deqn{\sum_{i = 1}^n(Y_i - X_i^\top \theta)X_i = 0.}
#'   [el_lm()] first computes the parameter estimates by calling [lm.fit()]
#'   (with `...` if any) with the `model.frame` and `model.matrix` obtained from
#'   the `formula`. Note that the maximum empirical likelihood estimator is the
#'   same as the least square estimator in our model. Next, it performs
#'   hypothesis tests based on an asymptotic chi-squared distribution of
#'   empirical likelihood ratio statistics. Included in the tests are overall
#'   test with
#'   \deqn{H_0: \theta_1 = \theta_2 = \cdots = \theta_{p-1} = 0,}
#'   and significance tests for each parameter with
#'   \deqn{H_{0j}: \theta_j = 0,\ j = 0, \dots, p-1.}
#'   The test results are returned as `optim` and `parTests`, respectively.
#' @return An object of class of \linkS4class{LM}.
#' @references Owen A (1991).
#'   “Empirical Likelihood for Linear Models.”
#'   The Annals of Statistics, 19(4), 1725–1747.
#'   \doi{10.1214/aos/1176348368}.
#' @seealso [el_control()], [el_glm()], [elt()]
#' @examples
#' set.seed(5649)
#' df <- data.frame(y = rnorm(50), x = rnorm(50))
#' fit <- el_lm(y ~ x, df)
#' summary(fit)
#'
#' fit2 <- el_lm(y ~ x, df, weights = rep(c(1, 2), each = 25))
#' summary(fit2)
#'
#' df[1, 2] <- NA
#' fit3 <- el_lm(y ~ x, df, na.action = na.omit)
#' summary(fit3)
#' @importFrom stats .getXlevels is.empty.model lm.fit lm.wfit model.matrix
#'   model.response model.weights pchisq setNames
#' @export
#' @srrstats {G2.14, G2.14a, RE2.1, RE2.2} Missing values are handled by the
#'   `na.action` argument via `na.cation()`. `Inf` values are not allowed and
#'   produce an error.Partially missing values (missing response or missing
#'   predictors) are allowed unless a singular fit is encountered. Although
#'   singular fits can produce estimates and fitted values, there is no
#'   practical advantage of using the package for singular fits since further
#'   inference based on empirical likelihood is unavailable. Note that at least
#'   for the linear models the maximum empirical likelihood estimates (and thus
#'   the fitted values as well) are identical to the estimates returned by
#'   `lm.fit()`.
#' @srrstats {G5.8, G5.8d} Data with more fields than observations produces an
#'   error.
#' @srrstats {G2.11, G2.12} Given the `formula` and `data` input,
#'   `stats::model.frame()` (which calls `as.data.frame()`) checks whether the
#'   `data` input is valid or not. If not, an error is triggered. The
#'   requirement is documented and tested. See `test/testthat/el_lm.R` as well.
#' @srrstats {RE1.0, RE1.1} Formula interface is used to the `formula` argument,
#'   and how it is converted to a matrix input is documented as well.
#' @srrstats {RE1.2} The expected format for the argument `data` is documented.
#' @srrstats {RE1.3, RE1.3a} The transformation only applies to the `data`
#'   argument. The only attributes that are passed to the output are the row and
#'   column names, and these are preserved. See `test/testthat/test-el_lm.R` as
#'   well.
#' @srrstats {RE2.0} Documented the use of `model.frame` and `model.matrix`.
#'   Users should be well aware of the basic transformations done by `lm()`.
#' @srrstats {RE2.1} Missing values are handled by the `na.action` argument.
#'   `Inf` values are not allowed and produce an error.
#' @srrstats {RE2.4, RE2.4a, RE2.4b} Perfect collinearity is handled by
#'   `model.frame()`. Especially, perfect collinearity among predictor variables
#'   produces an error in `lm.fit()` since `singular.ok` is set to `FALSE`. This
#'   is because the underlying asymptotic empirical likelihood theory requires
#'   a full-rank covariance structure in order for a limiting argument to work.
#'   See `EL-class` documentation.
#' @srrstats {RE4.0} `el_lm()` returns an object of class `LM`.
el_lm <- function(formula,
                  data,
                  weights = NULL,
                  na.action,
                  control = el_control(),
                  ...) {
  stopifnot("Invalid `control` specified." = (is(control, "ControlEL")))
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
  stopifnot(
    "`weights` must be a numeric vector." =
      (isTRUE(is.null(w) || is.numeric(w))),
    "`weights` must be positive." = (isTRUE(is.null(w) || all(w > 0))),
    "`el_lm()` does not support multiple responses." = (isFALSE(is.matrix(y)))
  )
  if (is.empty.model(mt)) {
    x <- NULL
    mm <- cbind(y, x)
    return(new("LM",
      call = cl, terms = mt,
      misc = list(
        xlevels = .getXlevels(mt, mf),
        na.action = attr(mf, "na.action")
      ),
      optim = list(
        par = numeric(), lambda = numeric(), iterations = integer(),
        convergence = logical()
      )
    ))
  } else {
    x <- model.matrix(mt, mf, NULL)
    fit <- if (is.null(w)) {
      lm.fit(x, y, offset = NULL, singular.ok = FALSE, ...)
    } else {
      lm.wfit(x, y, w, offset = NULL, singular.ok = FALSE, ...)
    }
  }
  pnames <- names(fit$coefficients)
  intercept <- attr(mt, "intercept")
  mm <- cbind(y, x)
  n <- nrow(mm)
  p <- ncol(x)
  w <- validate_weights(w, n)
  names(w) <- if (length(w) != 0L) names(y) else NULL
  out <- test_LM(
    mm, fit$coefficients, intercept, control@maxit, control@maxit_l,
    control@tol, control@tol_l, control@step, control@th, control@nthreads, w
  )
  optim <- validate_optim(out$optim)
  names(optim$par) <- pnames
  df <- if (intercept && p > 1L) p - 1L else p
  pval <- pchisq(out$statistic, df = df, lower.tail = FALSE)
  if (control@verbose) {
    message(
      "Convergence ",
      if (out$optim$convergence) "achieved." else "failed."
    )
  }
  new("LM",
    parTests = lapply(out$par_tests, setNames, pnames), call = cl, terms = mt,
    misc = list(
      xlevels = .getXlevels(mt, mf),
      na.action = attr(mf, "na.action")
    ),
    optim = optim, logp = setNames(out$logp, names(y)), logl = out$logl,
    loglr = out$loglr, statistic = out$statistic, df = df, pval = pval,
    nobs = n, npar = p, weights = w, data = if (control@keep_data) mm else NULL,
    coefficients = fit$coefficients, method = "lm"
  )
}
