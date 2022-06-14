#' Empirical likelihood for generalized linear models
#'
#' Fits a generalized linear model with empirical likelihood.
#'
#' @param formula An object of class \code{"\link[stats]{formula}"} (or one that
#'   can be coerced to that class): a symbolic description of the model to be
#'   fitted.
#' @param family A description of the error distribution and link function to be
#'   used in the model. Only the result of a call to a family function is
#'   supported. See ‘Details’.
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
#' @param start Starting values for the parameters in the linear predictor.
#'   Defaults to \code{NULL} and is passed to \code{\link[stats]{glm.fit}}.
#' @param etastart Starting values for the linear predictor. Defaults to
#'   \code{NULL} and is passed to \code{\link[stats]{glm.fit}}.
#' @param mustart Starting values for the vector of means. Defaults to
#'   \code{NULL} and is passed to \code{\link[stats]{glm.fit}}.
#' @param ... Additional arguments to be passed to
#'   \code{\link[stats]{glm.control}}.
#' @details The available families and link functions are as follows:
#'   \itemize{
#'   \item{\code{gaussian}}{: \code{identity}, \code{log}, and \code{inverse}.}
#'   \item{\code{bimomial}}{: \code{logit}, \code{probit}, and \code{log}.}
#'   \item{\code{poisson}}{: \code{log}, \code{identity}, and \code{sqrt}.}
#'   }
#'   Included in the tests are the overall test with
#'   \deqn{H_0: \beta_1 = \beta_2 = \cdots = \beta_{p-1} = 0,}
#'   and the tests for each parameter with
#'   \deqn{H_{0j}: \beta_j = 0,\ j = 0, \dots, p-1.}
#'   The test results are returned as \code{optim} and \code{parTests},
#'   respectively.
#' @return An object of class of \linkS4class{GLM}.
#' @references Chen, Song Xi, and Hengjian Cui. 2003.
#'   “An Extended Empirical Likelihood for Generalized Linear Models.”
#'   Statistica Sinica 13: 69–81.
#' @seealso \link{el_control}, \link{el_lm}, \link{lht}
#' @examples
#' n <- 50
#' x <- rnorm(n)
#' x2 <- rnorm(n)
#' l <- -2 + 0.2 * x + 3 * x2
#' mu <- 1 / (1 + exp(-l))
#' y <- rbinom(n, 1, mu)
#' df <- data.frame(y, x, x2)
#' fit <- el_glm(y ~ x + x2, family = binomial, df)
#' summary(fit)
#' @importFrom stats gaussian glm.fit model.extract model.weights pchisq
#' @export
el_glm <- function(formula, family = gaussian, data, weights = NULL, na.action,
                   control = el_control(), model = TRUE, start = NULL,
                   etastart = NULL, mustart = NULL, ...) {
  cl <- match.call()
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
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
  if (is.matrix(Y)) {
    stop("'el_glm' does not support grouped data")
  }
  X <- if (!is.empty.model(mt)) {
    model.matrix(mt, mf, NULL)
  } else {
    matrix(, NROW(Y), 0L)
  }
  w <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w)) {
    stop("'weights' must be a numeric vector")
  }
  if (!is.null(w) && any(w < 0)) {
    stop("negative weights not allowed")
  }
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
  if (is.empty.model(mt)) {
    return(new("GLM",
      misc = list(
        call = cl, formula = formula, terms = mt, offset = NULL,
        control = glm_control, method = "glm.fit",
        contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, mf),
        na.action = attr(mf, "na.action")
      ),
      optim = list(
        par = numeric(), lambda = numeric(), iterations = integer(),
        convergence = logical()
      )
    ))
  }
  intercept <- attr(mt, "intercept") > 0L
  fit <- glm.fit(
    x = X, y = Y, weights = w, start = start, etastart = etastart,
    mustart = mustart, offset = NULL, family = family,
    control = glm_control, intercept = intercept,
    singular.ok = FALSE
  )
  method <- check_family(fit$family)
  mm <- cbind(fit$y, X)
  p <- ncol(X)
  w <- check_weights(w, nrow(mm))
  if (!is(control, "ControlEL")) {
    stop("invalid 'control' specified")
  }
  el <- glm_(
    method, mm, fit$coefficients, intercept,
    control@maxit, control@maxit_l, control@tol, control@tol_l,
    control@step, control@th, control@nthreads, w
  )
  df <- if (intercept && p > 1L) p - 1L else p
  pval <- pchisq(el$statistic, df = df, lower.tail = FALSE)
  new("GLM",
    parTests = el$parTests,
    misc = list(
      family = fit$family, iter = fit$iter, converged = fit$converged,
      boundary = fit$boundary, call = cl, formula = formula, terms = mt,
      offset = NULL, control = glm_control, method = "glm.fit",
      contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, mf),
      na.action = attr(mf, "na.action")
    ),
    optim = el$optim, logp = el$logp, logl = el$logl, loglr = el$loglr,
    statistic = el$statistic, df = df, pval = pval, npar = p, weights = w,
    data = if (model) mm else matrix(NA_real_, nrow = 0L, ncol = 0L),
    coefficients = fit$coefficients, method = method
  )
}











#' #' Empirical likelihood for generalized linear models
#' #'
#' #' Fits a generalized linear model with empirical likelihood.
#' #'
#' #' @param formula An object of class \code{"\link[stats]{formula}"} (or one that
#' #'   can be coerced to that class): a symbolic description of the model to be
#' #'   fitted.
#' #' @param family A description of the error distribution and link function to be
#' #'   used in the model. Only the result of a call to a family function is
#' #'   supported. See ‘Details’.
#' #' @param data An optional data frame, list or environment (or object coercible
#' #'   by \code{\link[base]{as.data.frame}} to a data frame) containing the
#' #'   variables in the formula. If not found in data, the variables are taken
#' #'   from \code{environment(formula)}.
#' #' @param weights An optional numeric vector of weights to be used in the
#' #'   fitting process. Defaults to \code{NULL}, corresponding to identical
#' #'   weights. If non-\code{NULL}, weighted empirical likelihood is computed.
#' #' @param na.action A function which indicates what should happen when the data
#' #'   contain \code{NA}s. The default is set by the \code{na.action} setting of
#' #'   \code{\link[base]{options}}, and is \code{na.fail} if that is unset.
#' #' @param control A list of control parameters set by \code{\link{el_control}}.
#' #' @param model A logical. If \code{TRUE} the data matrix used for fitting is
#' #'   returned.
#' #' @param start Starting values for the parameters in the linear predictor.
#' #'   Defaults to \code{NULL} and is passed to \code{\link[stats]{glm.fit}}.
#' #' @param etastart Starting values for the linear predictor. Defaults to
#' #'   \code{NULL} and is passed to \code{\link[stats]{glm.fit}}.
#' #' @param mustart Starting values for the vector of means. Defaults to
#' #'   \code{NULL} and is passed to \code{\link[stats]{glm.fit}}.
#' #' @param input Starting values for the vector of means. Defaults to
#' #'   \code{NULL} and is passed to \code{\link[stats]{glm.fit}}.
#' #' @param ... Additional arguments to be passed to
#' #'   \code{\link[stats]{glm.control}}.
#' #' @details The available families and link functions are as follows:
#' #'   \itemize{
#' #'   \item{\code{gaussian}}{: \code{identity}, \code{log}, and \code{inverse}.}
#' #'   \item{\code{bimomial}}{: \code{logit}, \code{probit}, and \code{log}.}
#' #'   \item{\code{poisson}}{: \code{log}, \code{identity}, and \code{sqrt}.}
#' #'   }
#' #'   Included in the tests are the overall test with
#' #'   \deqn{H_0: \beta_1 = \beta_2 = \cdots = \beta_{p-1} = 0,}
#' #'   and the tests for each parameter with
#' #'   \deqn{H_{0j}: \beta_j = 0,\ j = 0, \dots, p-1.}
#' #'   The test results are returned as \code{optim} and \code{parTests},
#' #'   respectively.
#' #' @return An object of class of \linkS4class{GLM}.
#' #' @references Chen, Song Xi, and Hengjian Cui. 2003.
#' #'   “An Extended Empirical Likelihood for Generalized Linear Models.”
#' #'   Statistica Sinica 13: 69–81.
#' #' @seealso \link{el_control}, \link{el_lm}, \link{lht}
#' #' @examples
#' #' n <- 50
#' #' x <- rnorm(n)
#' #' x2 <- rnorm(n)
#' #' l <- -2 + 0.2 * x + 3 * x2
#' #' mu <- 1 / (1 + exp(-l))
#' #' y <- rbinom(n, 1, mu)
#' #' df <- data.frame(y, x, x2)
#' #' fit <- el_glm(y ~ x + x2, family = binomial, df)
#' #' summary(fit)
#' #' @importFrom stats gaussian glm.fit model.extract model.weights pchisq
#' #' @export
#' el_glm2 <- function(formula, family = gaussian, data, weights = NULL, na.action,
#'                     control = el_control(), model = TRUE, start = NULL,
#'                     etastart = NULL, mustart = NULL, input = NULL, ...) {
#'   cl <- match.call()
#'   if (is.character(family)) {
#'     family <- get(family, mode = "function", envir = parent.frame())
#'   }
#'   if (is.function(family)) {
#'     family <- family()
#'   }
#'   if (is.null(family$family)) {
#'     print(family)
#'     stop("'family' not recognized")
#'   }
#'   if (missing(data)) {
#'     data <- environment(formula)
#'   }
#'   mf <- match.call(expand.dots = FALSE)
#'   m <- match(c(
#'     "formula", "data", "weights", "na.action", "etastart",
#'     "mustart"
#'   ), names(mf), 0L)
#'   mf <- mf[c(1L, m)]
#'   mf$drop.unused.levels <- TRUE
#'   mf[[1L]] <- quote(stats::model.frame)
#'   mf <- eval(mf, parent.frame())
#'   glm_control <- do.call("glm.control", list(...))
#'   mt <- attr(mf, "terms")
#'   Y <- model.response(mf, "any")
#'   if (length(dim(Y)) == 1L) {
#'     nm <- rownames(Y)
#'     dim(Y) <- NULL
#'     if (!is.null(nm)) {
#'       names(Y) <- nm
#'     }
#'   }
#'   if (is.matrix(Y)) {
#'     stop("'el_glm' does not support grouped data")
#'   }
#'   X <- if (!is.empty.model(mt)) {
#'     model.matrix(mt, mf, NULL)
#'   } else {
#'     matrix(, NROW(Y), 0L)
#'   }
#'   w <- as.vector(model.weights(mf))
#'   if (!is.null(w) && !is.numeric(w)) {
#'     stop("'weights' must be a numeric vector")
#'   }
#'   if (!is.null(w) && any(w < 0)) {
#'     stop("negative weights not allowed")
#'   }
#'   mustart <- model.extract(mf, "mustart")
#'   etastart <- model.extract(mf, "etastart")
#'   if (is.empty.model(mt)) {
#'     return(new("GLM",
#'                optim = list(
#'                  method = "GLM", par = numeric(), lambda = numeric(),
#'                  iterations = integer(), convergence = logical()
#'                ),
#'                misc = list(
#'                  call = cl, formula = formula, terms = mt,
#'                  offset = NULL, control = glm_control, method = "glm.fit",
#'                  contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, mf),
#'                  na.action = attr(mf, "na.action")
#'                )
#'     ))
#'   }
#'   intercept <- attr(mt, "intercept") > 0L
#'   fit <- glm.fit(
#'     x = X, y = Y, weights = w, start = start, etastart = etastart,
#'     mustart = mustart, offset = NULL, family = family,
#'     control = glm_control, intercept = intercept,
#'     singular.ok = FALSE
#'   )
#'   method <- check_family(fit$family)
#'   mm <- cbind(fit$y, X)
#'   p <- ncol(X)
#'   w <- check_weights(w, nrow(mm))
#'   if (!is(control, "ControlEL")) {
#'     stop("invalid 'control' specified")
#'   }
#'   disp <-
#'     mean((fit$y - fit$fitted.values)^2 / fit$family$variance(fit$fitted.values))
#'   tmp <- if (is.null(input)) fit$coefficients else input
#'   el <- glm2_(
#'     method$family, method$link, mm, c(tmp, disp), intercept,
#'     control@maxit, control@maxit_l, control@tol, control@tol_l,
#'     control@step, control@th, control@nthreads, w
#'   )
#'   el
#'   # df <- if (intercept && p > 1L) p - 1L else p
#'   # pval <- pchisq(el$statistic, df = df, lower.tail = FALSE)
#'   # new("GLM",
#'   #   optim = el$optim, logp = el$logp, logl = el$logl, loglr = el$loglr,
#'   #   statistic = el$statistic, df = df, pval = pval, npar = p, weights = w,
#'   #   data = if (model) mm else matrix(NA_real_, nrow = 0L, ncol = 0L),
#'   #   coefficients = fit$coefficients,
#'   #   parTests = el$parTests,
#'   #   misc = list(
#'   #     family = fit$family, iter = fit$iter, rdf = fit$df.residual,
#'   #     converged = fit$converged, boundary = fit$boundary, call = cl,
#'   #     formula = formula, terms = mt, offset = NULL, control = glm_control,
#'   #     method = "glm.fit", contrasts = attr(X, "contrasts"),
#'   #     xlevels = .getXlevels(mt, mf), na.action = attr(mf, "na.action")
#'   #   )
#'   # )
#' }
