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
#' @param offset An optional expression for specifying an \emph{a priori} known
#'   component to be included in the linear predictor during fitting. This
#'   should be `NULL` or a numeric vector or matrix of extents matching those of
#'   the response. One or more [`offset`] terms can be included in the formula
#'   instead or as well, and if more than one are specified their sum is used.
#' @param control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @param ... Additional arguments to be passed to the low level regression
#'   fitting functions. See ‘Details’.
#' @details Suppose that we observe \eqn{n} independent random variables
#'   \eqn{{Z_i} \equiv {(X_i, Y_i)}} from a common distribution, where \eqn{X_i}
#'   is the \eqn{p}-dimensional covariate (including the intercept if any) and
#'   \eqn{Y_i} is the response. We consider the following linear model:
#'   \deqn{Y_i = X_i^\top \theta + \epsilon_i,}
#'   where \eqn{\theta = (\theta_0, \dots, \theta_{p-1})} is an unknown
#'   \eqn{p}-dimensional parameter and the errors \eqn{\epsilon_i} are
#'   independent random variables that satisfy
#'   \eqn{\textrm{E}(\epsilon_i | X_i)} = 0. We assume that the errors have
#'   finite conditional variance. Then the least square estimator of
#'   \eqn{\theta} solves the following estimating equations:
#'   \deqn{\sum_{i = 1}^n(Y_i - X_i^\top \theta)X_i = 0.}
#'   Given a value of \eqn{\theta}, let
#'   \eqn{{g(Z_i, \theta)} = {(Y_i - X_i^\top \theta)X_i}} and the (profile)
#'   empirical likelihood ratio is defined by
#'   \deqn{R(\theta) =
#'   \max_{p_i}\left\{\prod_{i = 1}^n np_i :
#'   \sum_{i = 1}^n p_i g(Z_i, \theta) = \theta,\
#'    p_i \geq 0,\
#'   \sum_{i = 1}^n p_i = 1
#'   \right\}.}
#'   [el_lm()] first computes the parameter estimates by calling [lm.fit()]
#'   (with `...` if any) with the `model.frame` and `model.matrix` obtained from
#'   the `formula`. Note that the maximum empirical likelihood estimator is the
#'   same as the the quasi-maximum likelihood estimator in our model. Next, it
#'   tests hypotheses based on asymptotic chi-square distributions of the
#'   empirical likelihood ratio statistics. Included in the tests are overall
#'   test with
#'   \deqn{H_0: \theta_1 = \theta_2 = \cdots = \theta_{p-1} = 0,}
#'   and significance tests for each parameter with
#'   \deqn{H_{0j}: \theta_j = 0,\ j = 0, \dots, p-1.}
#' @return An object of class of \linkS4class{LM}.
#' @references Owen A (1991).
#'   “Empirical Likelihood for Linear Models.”
#'   \emph{The Annals of Statistics}, 19(4), 1725--1747.
#'   \doi{10.1214/aos/1176348368}.
#' @seealso \linkS4class{EL}, \linkS4class{LM}, [el_glm()], [elt()],
#'   [el_control()]
#' @examples
#' ## Linear model
#' data("thiamethoxam")
#' fit <- el_lm(fruit ~ trt, data = thiamethoxam)
#' summary(fit)
#'
#' ## Weighted data
#' wfit <- el_lm(fruit ~ trt, data = thiamethoxam, weights = visit)
#' summary(wfit)
#'
#' ## Missing data
#' fit2 <- el_lm(fruit ~ trt + scb, data = thiamethoxam,
#'   na.action = na.omit, offset = NULL
#' )
#' summary(fit2)
#' @export
el_lm <- function(formula,
                  data,
                  weights = NULL,
                  na.action,
                  offset,
                  control = el_control(),
                  ...) {
  stopifnot("Invalid `control` specified." = (is(control, "ControlEL")))
  cl <- match.call()
  if (missing(data)) {
    data <- environment(formula)
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match(
    c("formula", "data", "weights", "na.action", "offset"), names(mf), 0L
  )
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
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != length(y)) {
      stop(gettextf(
        "Number of offsets is %d, should equal %d (number of observations).",
        length(offset), length(y)
      ), domain = NA)
    }
  }
  if (is.empty.model(mt)) {
    x <- matrix(numeric(0), length(y), 0L)
    return(new("LM",
      call = cl, terms = mt,
      misc = list(
        intercept = FALSE,
        xlevels = .getXlevels(mt, mf),
        na.action = attr(mf, "na.action")
      ),
      optim = list(
        par = numeric(), lambda = numeric(), iterations = integer(),
        convergence = logical(), cstr = logical()
      ), df = 0L, nobs = nrow(x), npar = 0L, method = NA_character_,
      control = control
    ))
  } else {
    x <- model.matrix(mt, mf, NULL)
    fit <- if (is.null(w)) {
      lm.fit(x, y, offset = offset, singular.ok = FALSE, ...)
    } else {
      lm.wfit(x, y, w, offset = offset, singular.ok = FALSE, ...)
    }
  }
  pnames <- names(fit$coefficients)
  intercept <- as.logical(attr(mt, "intercept"))
  s <- if (is.null(offset)) rep.int(0, length(y)) else offset
  mm <- cbind(s, y, x)
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
  optim$cstr <- intercept
  df <- if (intercept && p > 1L) p - 1L else p
  pval <- pchisq(out$statistic, df = df, lower.tail = FALSE)
  if (control@verbose) {
    message(
      "Convergence ",
      if (out$optim$convergence) "achieved." else "failed."
    )
  }
  new("LM",
    sigTests = lapply(out$sig_tests, setNames, pnames), call = cl, terms = mt,
    misc = list(
      intercept = intercept, xlevels = .getXlevels(mt, mf),
      na.action = attr(mf, "na.action"), offset = offset
    ),
    optim = optim, logp = setNames(out$logp, names(y)), logl = out$logl,
    loglr = out$loglr, statistic = out$statistic, df = df, pval = pval,
    nobs = n, npar = p, weights = w, coefficients = fit$coefficients,
    method = "lm", data = if (control@keep_data) mm else NULL, control = control
  )
}
