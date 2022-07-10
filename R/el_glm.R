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
#' @param control An object of class \linkS4class{ControlEL} constructed by
#'   [el_control()].
#' @param start Starting values for the parameters in the linear predictor.
#'   Defaults to `NULL` and is passed to [glm.fit()].
#' @param etastart Starting values for the linear predictor. Defaults to `NULL`
#'   and is passed to [glm.fit()].
#' @param mustart Starting values for the vector of means. Defaults to `NULL`
#'   and is passed to [glm.fit()].
#' @param ... Additional arguments to be passed to [glm.control()].
#' @details The available families and link functions are as follows:
#'   \itemize{
#'   \item{`gaussian`}{: `identity`, `log`, and `inverse`.}
#'   \item{`bimomial`}{: `logit`, `probit`, and `log`.}
#'   \item{`poisson`}{: `log`, `identity`, and `sqrt`.}
#'   }
#'   Included in the tests are the overall test with
#'   \deqn{H_0: \beta_1 = \beta_2 = \cdots = \beta_{p-1} = 0,}
#'   and the tests for each parameter with
#'   \deqn{H_{0j}: \beta_j = 0,\ j = 0, \dots, p-1.}
#'   The test results are returned as `optim` and `parTests`, respectively.
#' @return An object of class of \linkS4class{GLM}.
#' @references Chen SX, Cui H (2003).
#'   “An Extended Empirical Likelihood for Generalized Linear Models.”
#'   Statistica Sinica, 13(1), 69–81.
#' @references Kolaczyk ED (1994).
#'   “Empirical Likelihood for Generalized Linear Models.”
#'   Statistica Sinica, 4(1), 199–218.
#' @seealso [el_control()], [el_lm()], [elt()]
#' @examples
#' set.seed(20010)
#' n <- 50
#' x <- rnorm(n)
#' x2 <- rnorm(n)
#' l <- -2 + 0.2 * x + 3 * x2
#' mu <- 1 / (1 + exp(-l))
#' y <- rbinom(n, 1, mu)
#' df <- data.frame(y, x, x2)
#' fit <- el_glm(y ~ x + x2,
#'   family = binomial, df, weights = NULL,
#'   na.action = na.omit, start = NULL, etastart = NULL, mustart = NULL
#' )
#' summary(fit)
#' @importFrom stats gaussian glm.fit model.extract model.weights pchisq
#'   setNames
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
#'   `glm.fit()`.
#' @srrstats {G5.8, G5.8d} Data with more fields than observations produces an
#'   error.
#' @srrstats {G2.11, G2.12} Given the `formula` and `data` input,
#'   `stats::model.frame()` (which calls `as.data.frame()`) checks whether the
#'   `data` input is valid or not. If not, an error is triggered. The
#'   requirement is documented and tested.
#' @srrstats {RE1.0, RE1.1} Formula interface is used to the `formula` argument,
#'   and how it is converted to a matrix input is documented as well.
#' @srrstats {RE1.2} The expected format for the argument `data` is documented.
#' @srrstats {RE1.3, RE1.3a} The transformation only applies to the `data`
#'   argument. The only attributes that are passed to the output are the row and
#'   column names, and these are preserved.
#' @srrstats {RE2.0} Documented the use of `model.frame` and `model.matrix`.
#'   Users should be well aware of the basic transformations done by `glm()`.
#' @srrstats {RE2.1} Missing values are handled by the `na.action` argument.
#'   `Inf` values are not allowed and produce an error.
#' @srrstats {RE2.4, RE2.4a, RE2.4b} Perfect collinearity is handled by
#'   `model.frame()`. Especially, perfect collinearity among predictor variables
#'   produces an error in `glm.fit()` since `singular.ok` is set to `FALSE`.
#'   This is because the underlying asymptotic empirical likelihood theory
#'   requires a full-rank covariance structure in order for a limiting argument
#'   to work. See `EL-class` documentation.
#' @srrstats {RE4.0} `el_glm()` returns an object of class `GLM`.
el_glm <- function(formula,
                   family = gaussian,
                   data,
                   weights = NULL,
                   na.action,
                   control = el_control(),
                   start = NULL,
                   etastart = NULL,
                   mustart = NULL,
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
    "formula", "data", "weights", "na.action", "etastart",
    "mustart"
  ), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  glm_control <- do.call("glm.control", list(...))
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm)) {
      names(Y) <- nm
    }
  }
  stopifnot(
    "`el_glm()` does not support grouped data." = (isFALSE(is.matrix(Y)))
  )
  if (is.empty.model(mt)) {
    X <- matrix(, NROW(Y), 0L)
    return(new("GLM",
      call = cl, terms = mt,
      misc = list(
        formula = formula, offset = NULL, control = glm_control,
        intercept = FALSE, method = "glm.fit", contrasts = attr(X, "contrasts"),
        xlevels = .getXlevels(mt, mf), na.action = attr(mf, "na.action")
      ),
      optim = list(
        par = numeric(), lambda = numeric(), iterations = integer(),
        convergence = logical()
      )
    ))
  } else {
    X <- model.matrix(mt, mf, NULL)
  }
  w <- as.vector(model.weights(mf))
  stopifnot(
    "`weights` must be a numeric vector." =
      (isTRUE(is.null(w) || is.numeric(w))),
    "`weights` must be positive." = (isTRUE(is.null(w) || all(w > 0)))
  )
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
  intercept <- attr(mt, "intercept") > 0L
  fit <- glm.fit(
    x = X, y = Y, weights = w, start = start, etastart = etastart,
    mustart = mustart, offset = NULL, family = family,
    control = glm_control, intercept = intercept,
    singular.ok = FALSE
  )
  pnames <- names(fit$coefficients)
  method <- validate_family(fit$family)
  mm <- cbind(fit$y, X)
  n <- nrow(mm)
  p <- ncol(X)
  w <- validate_weights(w, n)
  names(w) <- if (length(w) != 0L) names(Y) else NULL
  out <- test_GLM(
    method, mm, fit$coefficients, intercept, control@maxit, control@maxit_l,
    control@tol, control@tol_l,control@step, control@th, control@nthreads, w
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
  new("GLM",
    parTests = lapply(out$par_tests, setNames, pnames), call = cl, terms = mt,
    misc = list(
      family = fit$family, iter = fit$iter, converged = fit$converged,
      boundary = fit$boundary, formula = formula, offset = NULL,
      control = glm_control, intercept = intercept, method = "glm.fit",
      contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, mf),
      na.action = attr(mf, "na.action")
    ),
    optim = optim, logp = setNames(out$logp, names(Y)), logl = out$logl,
    loglr = out$loglr, statistic = out$statistic, df = df, pval = pval,
    nobs = n, npar = p, weights = w, data = if (control@keep_data) mm else NULL,
    coefficients = fit$coefficients, method = method
  )
}
