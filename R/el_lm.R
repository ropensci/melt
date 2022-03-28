#' Fit a linear model with empirical likelihood
#'
#' Fit a linear model with empirical likelihood.
#'
#' @param formula A formula object.
#' @param data A data frame containing the variables in the formula.
#' @param na.action A function which indicates what should happen when the data
#'   contain \code{NA}s.
#' @param control A list of control parameters. See ‘Details’ in
#'   \code{\link{el_eval}}.
#' @param keep.data A logical. If \code{TRUE} the data matrix used in fitting is
#'   returned.
#' @inheritParams el_eval
#' @return A list with class \code{c("el_lm", "el_test")}.
#' @references Owen, Art. 1991. “Empirical Likelihood for Linear Models.”
#'   The Annals of Statistics 19 (4).
#'   \doi{10.1214/aos/1176348368}.
#' @seealso \link{el_eval}
#' @examples
#' fit <- el_lm(formula = mpg ~ wt, data = mtcars)
#' summary(fit)
#' @importFrom stats .getXlevels is.empty.model model.matrix model.response setNames
#' @export
el_lm <- function(formula, data, weights = NULL, na.action, control = list(),
                  keep.data = TRUE) {
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  action <- if (missing(na.action) && is.null(attr(mf, "na.action"))) list(NULL)
  else attr(mf, "na.action")
  if (is.null(action))
    action <- list(NULL)

  mt <- attr(mf, "terms")
  intercept <- attr(mt, "intercept")
  y <- model.response(mf, "numeric")
  nm <- names(y)
  if (is.matrix(y))
    stop("'el_lm' does not support multiple responses")
  x <- model.matrix(mt, mf)
  mm <- cbind(y, x)

  if (is.empty.model(mt)) {
    out <- list(coefficients = numeric(), residuals = y, rank = 0L,
                fitted.values = 0 * y, na.action = action,
                xlevels = .getXlevels(mt, mf), call = cl, terms = mt)
    if (keep.data)
      out$data.matrix <- mm
    class(out) <- c("el_lm", "el_test")
    return(out)
  }

  if (is.null(weights)) {
    w <- NULL
  } else {
    w <- check_weights(weights, NROW(mm))
  }

  optcfg <- check_control(control)
  out <- lm_(mm, intercept, optcfg$maxit, optcfg$tol, optcfg$th)
  out$coefficients <- setNames(out$coefficients, colnames(x))
  out$residuals <- setNames(out$residuals, nm)
  out$fitted.values <- setNames(out$fitted.values, nm)
  out$na.action <- action
  out$xlevels <- .getXlevels(mt, mf)
  out$call <- cl
  out$terms <- mt
  if (keep.data)
    out$data.matrix <- mm
  out
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

#' @importFrom stats formula
#' @export
logLik.el_lm <- function(object, ...) {
  mt <- object$terms
  mf <- object$model
  x <- model.matrix(mt, mf)
  y <- model.response(mf, "numeric")
  mele <- object$coefficients

  res <- object$residuals
  p <- object$rank
  N <- length(res)

  # no support for weights
  # if (is.null(w <- object$weights)) {
  #   w <- rep.int(1, N)
  # }
  # else {
  #   excl <- w == 0
  #   if (any(excl)) {
  #     res <- res[!excl]
  #     N <- length(res)
  #     w <- w[!excl]
  #   }
  # }

  N0 <- N
  val <- 0 - N * log(N)
  attr(val, "nall") <- N0
  attr(val, "nobs") <- N
  # df is p instead of p + 1 since EL does not estimate variance
  attr(val, "df") <- p
  class(val) <- "logLik"
  val
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
  p <- z$rank
  if (p == 0) {
    r <- z$residuals
    n <- length(r)
    ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
    ans$coefficients <-
      matrix(NA_real_, 0L, 3L, dimnames =
               list(NULL, c("estimate", "chisq-value", "p-value")))
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
                            `chisq-value` = z$optim$par.tests$statistic,
                            `p-value` = z$optim$par.tests$p.value)
  ans$aliased <- is.na(z$coefficients)
  if (p != attr(z$terms, "intercept")) {
    df.int <- if (attr(z$terms, "intercept")) 1L else 0L
    ans$chisq.statistic <- c(value = -2 * z$optim$logLR, df = p - df.int)
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
