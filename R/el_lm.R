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
#' @param control A list of control parameters set by
#'   \code{\link{control_el}}.
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
#'   The test results are returned as \code{optim} and \code{par.tests},
#'   respectively.
#' @return A list of class \code{c("el_lm", "el")} with the following
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
#'   \item{par}{The value at which empricial likelihood is evaluated.}
#'   \item{npar}{The number of parameters.}
#'   \item{weights}{The rescaled weights if non-\code{NULL} \code{weights} is
#'   supplied}
#'   \item{data.matrix}{The data matrix used for fitting if \code{model} is
#'   \code{TRUE}.}
#'   \item{coefficients}{The maximum empirical likelihood estimates of the
#'   parameters.}
#' @references Owen, Art. 1991. “Empirical Likelihood for Linear Models.”
#'   The Annals of Statistics 19 (4): 1725–47. \doi{10.1214/aos/1176348368}.
#' @seealso \link{el_glm}, \link{control_el}, \link{lht}
#' @examples
#' fit <- el_lm(mpg ~ wt, mtcars)
#' summary(fit)
#' @importFrom stats .getXlevels is.empty.model lm.fit lm.wfit model.matrix
#'   model.response model.weights
#' @export
el_lm <- function(formula, data, weights = NULL, na.action,
                  control = control_el(), model = TRUE, ...) {
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
    out <- list(
      optim = list(), log.prob = numeric(), loglik = numeric(),
      statistic = numeric(), df = 0L, p.value = numeric(), npar = 0L,
      coefficients = numeric(), na.action = attr(mf, "na.action"),
      xlevels = .getXlevels(mt, mf), call = cl, terms = mt
    )
    class(out) <- c("el_lm", "el")
    return(out)
  } else {
    x <- model.matrix(mt, mf, NULL)
    z <- if (is.null(w)) {
      lm.fit(x, y, offset = NULL, singular.ok = FALSE, ...)
    } else {
      lm.wfit(x, y, w, offset = NULL, singular.ok = FALSE, ...)
    }
  }
  if (!inherits(control, "control_el") || !is.list(control)) {
    stop("invalid 'control' supplied")
  }
  intercept <- attr(mt, "intercept")
  mm <- cbind(y, x)
  p <- ncol(x)
  w <- check_weights(w, nrow(mm))
  out <- lm_(mm, z$coefficients, intercept, control$maxit, control$maxit_l,
             control$tol, control$tol_l, control$step, control$th,
             control$nthreads, w)
  out$df <- if (intercept && p > 1L) p - 1L else p
  out$p.value <- pchisq(out$statistic, df = out$df, lower.tail = FALSE)
  out$npar <- p
  if (!is.null(weights)) {
    out$weights <- w
  }
  if (model) {
    out$data.matrix <- mm
  }
  structure(c(out, list(
    coefficients = z$coefficients,
    na.action = attr(mf, "na.action"),
    xlevels = .getXlevels(mt, mf), call = cl, terms = mt
  )),
  class = c(out$class, c("el_lm", "el"))
  )
}

#' @exportS3Method print el_lm
print.el_lm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n",
      sep = "")
  if (length(x$coefficients) == 0L) {
    cat("No coefficients\n")
  } else {
    cat("Coefficients:\n")
    print.default(format(x$coefficients, digits = digits), print.gap = 2L,
                  quote = FALSE)
  }
  cat("\n")
  invisible(x)
}

#' @export
summary.el_lm <- function(object, ...) {
  if (!inherits(object, "el_lm")) {
    stop("invalid 'el_lm' object")
  }
  z <- object
  p <- z$npar
  if (p == 0L) {
    ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
    ans$coefficients <-
      matrix(NA_real_, 0L, 3L,
             dimnames = list(NULL, c("Estimate", "Chisq", "Pr(>Chisq)")))
    ans$aliased <- is.na(z$coefficients)
    class(ans) <- "summary.el_lm"
    return(ans)
  }
  if (is.null(z$terms)) {
    stop("invalid 'el_lm' object:  no 'terms' component")
  }
  ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
  ans$coefficients <- cbind(
    Estimate = z$coefficients,
    Chisq = z$par.tests$statistic,
    `Pr(>Chisq)` = pchisq(z$par.tests$statistic, df = 1L, lower.tail = FALSE)
  )
  ans$aliased <- is.na(z$coefficients)
  if (p != attr(z$terms, "intercept")) {
    # ans$chisq.statistic <- c(value = -2 * z$optim$logLR, df = z$df)
    ans$chisq.statistic <- c(value = z$statistic, df = z$df)
  }
  if (!is.null(z$na.action))
    ans$na.action <- z$na.action
  class(ans) <- "summary.el_lm"
  ans
}

#' @importFrom stats naprint pchisq
#' @exportS3Method print summary.el_lm
print.summary.el_lm <- function(x, digits = max(3L, getOption("digits") - 3L),
                                signif.stars = getOption("show.signif.stars"),
                                ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n",
    sep = ""
  )
  if (length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n")
  } else {
    cat("\nCoefficients:\n")
    coefs <- x$coefficients
    if (any(aliased <- x$aliased)) {
      cn <- names(aliased)
      coefs <-
        matrix(NA, length(aliased), 3L, dimnames = list(cn, colnames(coefs)))
      coefs[!aliased, ] <- x$coefficients
    }
    printCoefmat(coefs,
      digits = digits, signif.stars = signif.stars,
      P.values = TRUE, has.Pvalue = TRUE, na.print = "NA", ...
    )
  }
  cat("\n")
  if (nzchar(mess <- naprint(x$na.action))) {
    cat("  (", mess, ")\n", sep = "")
  }
  if (!is.null(x$chisq.statistic)) {
    out <- c(
      paste("Chisq:", format(x$chisq.statistic[1L], digits = digits)),
      paste("df:", x$chisq.statistic[2L]),
      paste("p-value:", format.pval(
        pchisq(x$chisq.statistic[1L], x$chisq.statistic[2L],
               lower.tail = FALSE), digits = digits)))
    cat(strwrap(paste(out, collapse = ", ")), "\n\n")
  }
  invisible(x)
}
