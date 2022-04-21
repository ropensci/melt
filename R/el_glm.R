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
#' @param control A list of control parameters set by
#'   \code{\link{control_el}}.
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
#' @return A list of class \code{c("el_glm", "el_lm", "el")} with the following
#'   components:
#'   \item{optim}{A list with the following optimization results:
#'     \itemize{
#'       \item{\code{method } }{A character for method dispatch in internal
#'       functions.}
#'       \item{\code{par } }{The solution of the constrained minimization
#'       problem.}
#'       \item{\code{lambda } }{The Lagrange multiplier of dual problem.}
#'       \item{\code{logLR } }{The (weighted) empirical log-likelihood ratio
#'       value.}
#'       \item{\code{iterations } }{The number of iterations performed.}
#'       \item{\code{convergence } }{A logical vector. \code{TRUE} indicates
#'       convergence of the algorithm.}
#'     }
#'   }
#'   \item{par.tests}{A list with the test results for each parameter:
#'     \itemize{
#'       \item{\code{statistic } }{A numeric vector of chi-squared statistics.}
#'       \item{\code{convergence } }{A logical vector. \code{TRUE} indicates
#'       convergence of the algorithm.}
#'     }
#'   }
#'   \item{log.prob}{The log probabilities.}
#'   \item{loglik}{The log likelihood value.}
#'   \item{statistic}{The chi-square statistic.}
#'   \item{df}{The degrees of freedom of the statistic.}
#'   \item{p.value}{The \eqn{p}-value of the statistic.}
#'   \item{npar}{The number of parameters.}
#'   \item{weights}{The rescaled weights if non-\code{NULL} \code{weights} is
#'   supplied}
#'   \item{data.matrix}{The data matrix used for fitting if \code{model} is
#'   \code{TRUE}.}
#'   \item{coefficients}{The maximum empirical likelihood estimates of the
#'   parameters.}
#' @references Chen, Song Xi, and Hengjian Cui. 2003.
#'   “An Extended Empirical Likelihood for Generalized Linear Models.”
#'   Statistica Sinica 13: 69–81.
#' @seealso \link{el_lm}, \link{control_el}, \link{lht}
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
#' @importFrom stats gaussian glm.fit model.extract model.weights
#' @export
el_glm <- function(formula, family = gaussian, data, weights = NULL, na.action,
                   control = control_el(), model = TRUE, start = NULL,
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
    out <- list(
      optim = list(), log.prob = numeric(), loglik = numeric(),
      statistic = numeric(), df = 0L, p.value = numeric(), npar = 0L
    )
    out$na.action <- attr(mf, "na.action")
    return(structure(c(out, list(
      coefficients = numeric(), call = cl,
      formula = formula, terms = mt, offset = NULL,
      control = glm_control, method = "glm.fit",
      contrasts = attr(X, "contrasts"),
      xlevels = .getXlevels(mt, mf)
    )),
    class = c(out$class, c("el_glm", "el_lm", "el"))
    ))
  }

  if (!inherits(control, "control_el") || !is.list(control)) {
    stop("invalid 'control' supplied")
  }
  intercept <- attr(mt, "intercept") > 0L
  fit <- glm.fit(x = X, y = Y, weights = w, start = start, etastart = etastart,
                 mustart = mustart, offset = NULL, family = family,
                 control = glm_control, intercept = intercept,
                 singular.ok = FALSE)
  method <- check_family(fit$family)
  mm <- cbind(fit$y, X)
  p <- ncol(X)
  w <- check_weights(w, nrow(mm))
  out <- glm_(method$family, method$link, mm, fit$coefficients, intercept,
              control$maxit, control$maxit_l, control$tol, control$tol_l,
              control$step, control$th, control$nthreads, w)
  out$df <- if (intercept && p > 1L) p - 1L else p
  out$p.value <- pchisq(out$statistic, df = out$df, lower.tail = FALSE)
  out$npar <- p
  if (!is.null(weights)) {
    out$weights <- w
  }
  if (model) {
    out$data.matrix <- mm
  }
  out$na.action <- attr(mf, "na.action")
  structure(c(out, list(
    coefficients = fit$coefficients, family = fit$family,
    iter = fit$iter, converged = fit$converged,
    boundary = fit$boundary, call = cl, formula = formula,
    terms = mt, offset = NULL, control = glm_control,
    method = "glm.fit", contrasts = attr(X, "contrasts"),
    xlevels = .getXlevels(mt, mf)
  )),
  class = c(out$class, c("el_glm", "el_lm", "el"))
  )
}
