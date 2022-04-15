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
#'   See ‘Details’.
#' @param na.action A function which indicates what should happen when the data
#'   contain \code{NA}s. The default is set by the \code{na.action} setting of
#'   \code{\link[base]{options}}, and is \code{na.fail} if that is unset.
#' @param control A list of control parameters set by
#'   \code{\link{melt_control}}.
#' @param model A logical. If \code{TRUE} the model matrix used for fitting is
#'   returned.
#' @param ... Additional arguments to be passed to the low level regression
#'   fitting functions. See ‘Details’.
#' @return A list of class \code{c("el_lm", "el")}.
#' @references Owen, Art. 1991. “Empirical Likelihood for Linear Models.”
#'   The Annals of Statistics 19 (4): 1725–47. \doi{10.1214/aos/1176348368}.
#' @seealso \link{melt_control}, \link{lht}
#' @examples
#' fit <- el_lm(mpg ~ wt, mtcars)
#' summary(fit)
#' @importFrom stats .getXlevels is.empty.model lm.fit lm.wfit model.matrix
#'   model.response model.weights setNames
#' @export
el_lm <- function(formula, data, weights = NULL, na.action,
                  control = melt_control(), model = TRUE, ...) {
  cl <- match.call()
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  w <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w))
    stop("'weights' must be a numeric vector")
  if (!is.null(w) && any(w < 0))
    stop("negative weights not allowed")
  if (is.matrix(y))
    stop("'el_lm' does not support multiple responses")
  if (is.empty.model(mt)) {
    x <- NULL
    z <- list(optim = list(), log.prob = numeric(), loglik = numeric(),
              statistic = numeric(), df = 0L, p.value = numeric(), npar = 0L,
              coefficients = numeric(), na.action = attr(mf, "na.action"),
              xlevels = .getXlevels(mt, mf), call = cl, terms = mt)
    class(z) <- c("el_lm", "el")
    return(z)
  } else {
    x <- model.matrix(mt, mf, NULL)
    z <- if (is.null(w)) lm.fit(x, y, offset = NULL, singular.ok = FALSE, ...)
    else lm.wfit(x, y, w, offset = NULL, singular.ok = FALSE, ...)
  }

  if (!inherits(control, "melt_control") || !is.list(control))
    stop("invalid 'control' supplied")
  intercept <- attr(mt, "intercept")
  mm <- cbind(y, x)
  p <- ncol(x)
  w <- check_weights(w, nrow(mm))
  out <- lm_(mm, z$coefficients, intercept, control$maxit, control$maxit_l,
             control$tol, control$tol_l, control$th, control$nthreads, w)
  out$df <- if (intercept && p > 1L) p - 1L else p
  out$p.value <- pchisq(out$statistic, df = out$df, lower.tail = FALSE)
  out$npar <- p
  if (!is.null(weights))
    out$weights <- w
  if (model)
    out$data.matrix <- mm
  structure(c(out, list(coefficients = z$coefficients,
                        na.action = attr(mf, "na.action"),
                        xlevels = .getXlevels(mt, mf), call = cl, terms = mt)),
            class = c(out$class, c("el_lm", "el")))
}

#' @importFrom stats formula
#' @export
formula.el_lm <- function(x, ...) {
  form <- x$formula
  if (!is.null(form)) {
    form <- formula(x$terms)
    environment(form) <- environment(x$formula)
    form
  }
  else formula(x$terms)
}

#' @importFrom stats coef
#' @export
print.el_lm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  if (length(coef(x)) == 0L) {
    cat("No coefficients\n")
  } else {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2L,
                  quote = FALSE)
  }
  cat("\n")
  invisible(x)
}

#' @export
summary.el_lm <- function(object, ...) {
  z <- object
  p <- z$npar
  if (p == 0L) {
    r <- z$residuals
    n <- length(r)
    ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
    ans$coefficients <-
      matrix(NA_real_, 0L, 3L, dimnames =
               list(NULL, c("estimate", "chisq value", "p-value")))
    ans$aliased <- is.na(coef(object))
    ans$df <- c(0L, n, length(ans$aliased))
    class(ans) <- "summary.el_lm"
    return(ans)
  }
  if (!inherits(object, "el_lm"))
    stop("invalid 'el_lm' object")
  if (is.null(z$terms))
    stop("invalid 'el_lm' object:  no 'terms' component")
  ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
  ans$coefficients <- cbind(estimate = z$coefficients,
                            `chisq value` = z$par.tests$statistic,
                            # `Pr(>Chisq)` = z$optim$par.tests$p.value
                            `p-value` = pchisq(z$par.tests$statistic, df = 1L,
                                               lower.tail = FALSE))
  ans$aliased <- is.na(z$coefficients)
  if (p != attr(z$terms, "intercept")) {
    # df.int <- if (attr(z$terms, "intercept")) 1L else 0L
    # ans$chisq.statistic <- c(value = -2 * z$optim$logLR, df = p - df.int)
    ans$chisq.statistic <- c(value = -2 * z$optim$logLR, df = z$df)
  }
  if (!is.null(z$na.action)) ans$na.action <- z$na.action
  class(ans) <- "summary.el_lm"
  ans
}

#' @importFrom stats naprint pchisq
#' @export
print.summary.el_lm <- function(x, digits = max(3L, getOption("digits") - 3L),
                                signif.stars = getOption("show.signif.stars"),
                                ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
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
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                 P.values = TRUE, has.Pvalue = TRUE, na.print = "NA", ...)
  }
  cat("\n")
  if (nzchar(mess <- naprint(x$na.action)))
    cat("  (", mess, ")\n", sep = "")
  if (!is.null(x$chisq.statistic)) {
    out <- c(paste("Chisq:", format(x$chisq.statistic[1L], digits = digits)),
             paste("df:", x$chisq.statistic[2L]),
             paste("p-value:", format.pval(
               pchisq(x$chisq.statistic[1L], x$chisq.statistic[2L],
                      lower.tail = FALSE), digits = digits)))
    cat(strwrap(paste(out, collapse = ", ")), "\n\n")
  }
  invisible(x)
}
