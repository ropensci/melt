#' Empirical likelihood for generalized linear models
#'
#' Fits a generalized linear model with empirical likelihood.
#'
#' @param formula An object of class [`formula`] (or one that can be coerced to
#'   that class): a symbolic description of the model to be fitted.
#' @param family A description of the error distribution and link function to be
#'   used in the model. Only the result of a call to a family function is
#'   supported. See ‘Details’.
#' @param data An optional data frame, list or environment (or object coercible
#'   by [as.data.frame()] to a data frame) containing the variables in the
#'   formula. If not found in data, the variables are taken from
#'   `environment(formula)`.
#' @param weights An optional numeric vector of weights to be used in the
#'   fitting process. Defaults to `NULL`, corresponding to identical weights. If
#'   non-`NULL`, weighted empirical likelihood is computed.
#' @param na.action A function which indicates what should happen when the data
#'   contain `NA`s. The default is set by the `na.action` setting of
#'   [`options`], and is `na.fail` if that is unset.
#' @param start Starting values for the parameters in the linear predictor.
#'   Defaults to `NULL` and is passed to [glm.fit()].
#' @param etastart Starting values for the linear predictor. Defaults to `NULL`
#'   and is passed to [glm.fit()].
#' @param mustart Starting values for the vector of means. Defaults to `NULL`
#'   and is passed to [glm.fit()].
#' @param offset An optional expression for specifying an \emph{a priori} known
#'   component to be included in the linear predictor during fitting. This
#'   should be `NULL` or a numeric vector or matrix of extents matching those of
#'   the response. One or more [`offset`] terms can be included in the formula
#'   instead or as well, and if more than one are specified their sum is used.
#' @param control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @param ... Additional arguments to be passed to [glm.control()].
#' @details Suppose that we observe \eqn{n} independent random variables
#'   \eqn{{Z_i} \equiv {(X_i, Y_i)}} from a common distribution, where \eqn{X_i}
#'   is the \eqn{p}-dimensional covariate (including the intercept if any) and
#'   \eqn{Y_i} is the response. A generalized linear model specifies that
#'   \eqn{{\textrm{E}(Y_i | X_i)} = {\mu_i}},
#'   \eqn{{G(\mu_i)} = {X_i^\top \theta}}, and
#'   \eqn{{\textrm{Var}(Y_i | X_i)} = {\phi V(\mu_i)}},
#'   where \eqn{\theta = (\theta_0, \dots, \theta_{p-1})} is an unknown
#'   \eqn{p}-dimensional parameter, \eqn{\phi} is an optional dispersion
#'   parameter, \eqn{G} is a known smooth link function, and \eqn{V} is a known
#'   variance function.
#'
#'   With \eqn{H} denoting the inverse link function, define the quasi-score
#'   \deqn{{g_1(Z_i, \theta)} =
#'   \left\{
#'   H^\prime(X_i^\top \theta) \left(Y_i - H(X_i^\top \theta)\right) /
#'   \left(\phi V\left(H(X_i^\top \theta)\right)\right)
#'   \right\}
#'   X_i.}
#'   Then we have the estimating equations
#'   \eqn{\sum_{i = 1}^n g_1(Z_i, \theta) = 0}.
#'   When \eqn{\phi} is known, the (profile) empirical likelihood ratio for a
#'   given \eqn{\theta} is defined by
#'   \deqn{R_1(\theta) =
#'   \max_{p_i}\left\{\prod_{i = 1}^n np_i :
#'   \sum_{i = 1}^n p_i g_1(Z_i, \theta) = 0,\
#'   p_i \geq 0,\
#'   \sum_{i = 1}^n p_i = 1
#'   \right\}.}
#'   With unknown \eqn{\phi}, we introduce another estimating function based on
#'   the squared residuals. Let \eqn{{\eta} = {(\theta, \phi)}} and
#'   \deqn{{g_2(Z_i, \eta)} =
#'   \left(Y_i - H(X_i^\top \theta)\right)^2 /
#'   \left(\phi^2 V\left(H(X_i^\top \theta)\right)\right) - 1 / \phi.}
#'   Now the empirical likelihood ratio is defined by
#'   \deqn{R_2(\eta) =
#'   \max_{p_i}\left\{\prod_{i = 1}^n np_i :
#'   \sum_{i = 1}^n p_i g_1(Z_i, \eta) = 0,\
#'   \sum_{i = 1}^n p_i g_2(Z_i, \eta) = 0,\
#'   p_i \geq 0,\
#'   \sum_{i = 1}^n p_i = 1
#'   \right\}.}
#'   [el_glm()] first computes the parameter estimates by calling [glm.fit()]
#'   (with `...` if any) with the `model.frame` and `model.matrix` obtained from
#'   the `formula`. Note that the maximum empirical likelihood estimator is the
#'   same as the the quasi-maximum likelihood estimator in our model. Next, it
#'   tests hypotheses based on asymptotic chi-square distributions of the
#'   empirical likelihood ratio statistics. Included in the tests are overall
#'   test with
#'   \deqn{H_0: \theta_1 = \theta_2 = \cdots = \theta_{p-1} = 0,}
#'   and significance tests for each parameter with
#'   \deqn{H_{0j}: \theta_j = 0,\ j = 0, \dots, p-1.}
#'
#'   The available families and link functions are as follows:
#'   * `gaussian`: `"identity"`, `"log"`, and `"inverse"`.
#'   * `binomial`: `"logit"`, `"probit"`, and `"log"`.
#'   * `poisson`: `"log"`, `"identity"`, and `"sqrt"`.
#'   * `quasipoisson`: `"log"`, `"identity"`, and `"sqrt"`.
#' @return An object of class of \linkS4class{GLM}.
#' @references Chen SX, Cui H (2003).
#'   “An Extended Empirical Likelihood for Generalized Linear Models.”
#'   \emph{Statistica Sinica}, 13(1), 69--81.
#' @references Kolaczyk ED (1994).
#'   “Empirical Likelihood for Generalized Linear Models.”
#'   \emph{Statistica Sinica}, 4(1), 199--218.
#' @seealso \linkS4class{EL}, \linkS4class{GLM}, [el_lm()], [elt()],
#'   [el_control()]
#' @examples
#' data("warpbreaks")
#' fit <- el_glm(wool ~ .,
#'   family = binomial, data = warpbreaks, weights = NULL, na.action = na.omit,
#'   start = NULL, etastart = NULL, mustart = NULL, offset = NULL
#' )
#' summary(fit)
#' @export
el_glm <- function(formula,
                   family = gaussian,
                   data,
                   weights = NULL,
                   na.action,
                   start = NULL,
                   etastart = NULL,
                   mustart = NULL,
                   offset,
                   control = el_control(),
                   ...) {
  stopifnot("Invalid `control` specified." = (is(control, "ControlEL")))
  cl <- match.call()
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (is.null(family$family)) {
    print(family)
    stop("`family` not recognized.")
  }
  if (missing(data)) {
    data <- environment(formula)
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c(
    "formula", "data", "weights", "na.action", "etastart", "mustart", "offset"
  ), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  glm_control <- do.call("glm.control", list(...))
  mt <- attr(mf, "terms")
  y <- model.response(mf, "any")
  if (length(dim(y)) == 1L) {
    nm <- rownames(y)
    dim(y) <- NULL
    if (!is.null(nm)) {
      names(y) <- nm
    }
  }
  stopifnot(
    "`el_glm()` does not support grouped data." = (isFALSE(is.matrix(y)))
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
    if (grepl("quasi", family$family)) {
      class <- "QGLM"
      npar <- 1L
    } else {
      class <- "GLM"
      npar <- 0L
    }
    return(new(class,
      family = family, call = cl, terms = mt,
      misc = list(
        formula = formula, offset = offset, control = glm_control,
        intercept = FALSE, method = "glm.fit", contrasts = attr(x, "contrasts"),
        xlevels = .getXlevels(mt, mf), na.action = attr(mf, "na.action")
      ),
      optim = list(
        par = numeric(), lambda = numeric(), iterations = integer(),
        convergence = logical(), cstr = logical()
      ), df = 0L, nobs = nrow(x), npar = npar, method = NA_character_,
      control = control
    ))
  } else {
    x <- model.matrix(mt, mf, NULL)
  }
  w <- as.vector(model.weights(mf))
  stopifnot(
    "`weights` must be a numeric vector." =
      (isTRUE(is.null(w) || is.numeric(w))),
    "`weights` must be positive." = (isTRUE(is.null(w) || all(w > 0)))
  )
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
  intercept <- as.logical(attr(mt, "intercept"))
  fit <- glm.fit(
    x = x, y = y, weights = w, start = start, etastart = etastart,
    mustart = mustart, offset = offset, family = family,
    control = glm_control, intercept = intercept,
    singular.ok = FALSE
  )
  pnames <- names(fit$coefficients)
  method <- validate_family(fit$family)
  s <- if (is.null(offset)) rep.int(0, length(y)) else offset
  mm <- cbind(s, fit$y, x)
  n <- nrow(mm)
  p <- ncol(x)
  w <- validate_weights(w, n)
  names(w) <- if (length(w) != 0L) names(y) else NULL
  if (fit$family$family %in% c("poisson", "binomial")) {
    dispersion <- 1L
  } else {
    yhat <- fitted(fit)
    if (length(w) == 0L) {
      dispersion <- sum((fit$y - yhat)^2L / fit$family$variance(yhat)) / n
    } else {
      dispersion <- sum(w * ((fit$y - yhat)^2L / fit$family$variance(yhat))) / n
    }
  }
  if (grepl("quasi", method)) {
    class <- "QGLM"
    npar <- p + 1L
    out <- test_QGLM(
      method, mm, c(fit$coefficients, dispersion), intercept, control@maxit,
      control@maxit_l, control@tol, control@tol_l, control@step, control@th,
      control@nthreads, w
    )
    optim <- validate_optim(out$optim)
    names(optim$par) <- c(pnames, "phi")
    optim$cstr <- intercept
  } else {
    class <- "GLM"
    npar <- p
    out <- test_GLM(
      method, mm, fit$coefficients, intercept, control@maxit, control@maxit_l,
      control@tol, control@tol_l, control@step, control@th, control@nthreads, w
    )
    optim <- validate_optim(out$optim)
    names(optim$par) <- pnames
    optim$cstr <- intercept
  }
  df <- if (intercept && p > 1L) p - 1L else p
  pval <- pchisq(out$statistic, df = df, lower.tail = FALSE)
  if (control@verbose) {
    message(
      "Convergence ",
      if (out$optim$convergence) "achieved." else "failed."
    )
  }
  new(class,
    family = fit$family, dispersion = dispersion,
    sigTests = lapply(out$sig_tests, setNames, pnames), call = cl, terms = mt,
    misc = list(
      iter = fit$iter, converged = fit$converged, boundary = fit$boundary,
      formula = formula, offset = offset, control = glm_control,
      intercept = intercept, method = "glm.fit",
      contrasts = attr(x, "contrasts"), xlevels = .getXlevels(mt, mf),
      na.action = attr(mf, "na.action")
    ),
    optim = optim, logp = setNames(out$logp, names(y)), logl = out$logl,
    loglr = out$loglr, statistic = out$statistic, df = df, pval = pval,
    nobs = n, npar = npar, weights = w, coefficients = fit$coefficients,
    method = method, data = if (control@keep_data) mm else NULL,
    control = control
  )
}
