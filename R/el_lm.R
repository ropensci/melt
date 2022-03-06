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
#' @seealso \link{el_aov}, \link{el_eval}
#' @examples
#' fit <- el_lm(formula = mpg ~ wt, data = mtcars)
#' summary(fit)
#' @importFrom stats .getXlevels is.empty.model model.matrix model.response setNames terms
#' @export
el_lm <- function(formula, data, weights = NULL, na.action, control = list(), keep.data = TRUE)
{
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  if (is.empty.model(mt))
    stop("empty model specified")
  intercept <- attr(mt, "intercept")
  y <- model.response(mf, "numeric")
  if (is.matrix(y))
    stop("'el_lm' does not support multiple responses")
  nm <- names(y)
  x <- model.matrix(mt, mf)
  mm <- cbind(y, x)

  if (is.null(weights)) {
    w <- NULL
  } else {
    if (!is.numeric(weights))
      stop("'weights' must be a numeric vector")
    w <- as.numeric(weights)
    if (any(!is.finite(w)))
      stop("'weights' must be a finite numeric vector")
    if (any(w < 0))
      stop("negative 'weights' are not allowed")
    if (length(w) != NROW(mm))
      stop("'data' and 'weights' have incompatible dimensions")
  }

  action <- if (missing(na.action) && is.null(attr(mf, "na.action"))) list(NULL)
  else attr(mf, "na.action")
  if (is.null(action))
    action <- list(NULL)

  optcfg <- check_control(control)
  out <- EL_lm(mm, intercept, optcfg$maxit, optcfg$abstol,
               optcfg$threshold)
  out$coefficients <- setNames(out$coefficients, colnames(x))
  out$residuals <- setNames(out$residuals, nm)
  out$fitted.values <- setNames(out$fitted.values, nm)
  out$df.residual = nrow(x) - out$df
  out$na.action <- action
  out$xlevels <- .getXlevels(mt, mf)
  out$call <- cl
  out$terms <- mt
  # out$tmp <- w
  if (keep.data)
    out$data.matrix <- mm
  out
}

# case.names.el_lm <- function(object, full = FALSE, ...) {
#    w <- weights(object)
#   dn <- names(residuals(object))
#   if (full || is.null(w))
#     dn
#   else dn[w != 0]
# }

# dummy.coef.el_lm <- function(object, use.na = FALSE, ...) {
#   xl <- object$xlevels
#   if (!length(xl))
#     return(as.list(coef(object)))
#   Terms <- terms(object)
#   tl <- attr(Terms, "term.labels")
#   int <- attr(Terms, "intercept")
#   facs <- attr(Terms, "factors")[-1, , drop = FALSE]
#   Terms <- delete.response(Terms)
#   mf <- object$model %||% model.frame(object)
#   vars <- dimnames(facs)[[1]]
#   xtlv <- lapply(mf[, vars, drop = FALSE], levels)
#   nxl <- pmax(lengths(xtlv), 1L)
#   lterms <- apply(facs, 2L, function(x) prod(nxl[x > 0]))
#   nl <- sum(lterms)
#   args <- sapply(vars, function(i) if (nxl[i] == 1)
#     rep.int(1, nl)
#     else factor(rep.int(xtlv[[i]][1L], nl), levels = xtlv[[i]]),
#     simplify = FALSE)
#   dummy <- do.call(data.frame, args)
#   names(dummy) <- vars
#   pos <- 0L
#   rn <- rep.int(tl, lterms)
#   rnn <- character(nl)
#   for (j in tl) {
#     i <- vars[facs[, j] > 0]
#     ifac <- i[nxl[i] > 1]
#     lt.j <- lterms[[j]]
#     if (length(ifac) == 0L) {
#       rnn[pos + 1L] <- j
#     }
#     else {
#       p.j <- pos + seq_len(lt.j)
#       if (length(ifac) == 1L) {
#         dummy[p.j, ifac] <- x.i <- xtlv[[ifac]]
#         rnn[p.j] <- as.character(x.i)
#       }
#       else {
#         tmp <- expand.grid(xtlv[ifac], KEEP.OUT.ATTRS = FALSE)
#         dummy[p.j, ifac] <- tmp
#         rnn[p.j] <- apply(as.matrix(tmp), 1L, paste,
#                           collapse = ":")
#       }
#     }
#     pos <- pos + lt.j
#   }
#   attr(dummy, "terms") <- attr(mf, "terms")
#   lcontr <- object$contrasts
#   lci <- vapply(dummy, is.factor, NA)
#   lcontr <- lcontr[names(lci)[lci]]
#   mm <- model.matrix(Terms, dummy, lcontr, xl)
#   if (anyNA(mm)) {
#     warning("some terms will have NAs due to the limits of the method")
#     mm[is.na(mm)] <- NA
#   }
#   coef <- object$coefficients
#   if (!use.na)
#     coef[is.na(coef)] <- 0
#   asgn <- attr(mm, "assign")
#   res <- setNames(vector("list", length(tl)), tl)
#   if (isM <- is.matrix(coef)) {
#     for (j in seq_along(tl)) {
#       keep <- which(asgn == j)
#       cf <- coef[keep, , drop = FALSE]
#       ij <- rn == tl[j]
#       cf <- if (any(na <- is.na(cf))) {
#         if (ncol(cf) >= 2)
#           stop("multivariate case with missing coefficients is not yet implemented")
#         rj <- t(mm[ij, keep[!na], drop = FALSE] %*% cf[!na])
#         rj[apply(mm[ij, keep[na], drop = FALSE] != 0,
#                  1L, any)] <- NA
#         rj
#       }
#       else t(mm[ij, keep, drop = FALSE] %*% cf)
#       dimnames(cf) <- list(colnames(coef), rnn[ij])
#       res[[j]] <- cf
#     }
#   }
#   else {
#     for (j in seq_along(tl)) {
#       keep <- which(asgn == j)
#       cf <- coef[keep]
#       ij <- rn == tl[j]
#       res[[j]] <- if (any(na <- is.na(cf))) {
#         rj <- setNames(drop(mm[ij, keep[!na], drop = FALSE] %*%
#                               cf[!na]), rnn[ij])
#         rj[apply(mm[ij, keep[na], drop = FALSE] != 0,
#                  1L, any)] <- NA
#         rj
#       }
#       else setNames(drop(mm[ij, keep, drop = FALSE] %*%
#                            cf), rnn[ij])
#     }
#   }
#   if (int > 0)
#     res <- c(list(`(Intercept)` = if (isM) coef[int, ] else coef[int]),
#              res)
#   structure(res, class = "dummy_coef", matrix = isM)
# }


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

#' @importFrom stats nobs
#' @export
nobs.el_lm <- function(object, ...) {
  if (!is.null(w <- object$weights)) sum(w != 0) else NROW(object$residuals)
}

#' @importFrom stats formula
#' @export
logLik.el_lm <- function(object, ...) {
  mt <- object$terms
  mf <- object$model
  x <- model.matrix(mt, mf)
  y <- model.response(mf, "numeric")
  mele <- object$coefficients
  # out <- EL_lm2(x, mele, threshold = 500, maxit = 1e4, abstol = 1e-8, ...)

  res <- object$residuals
  p <- object$df
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
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits),
                  print.gap = 2L, quote = FALSE)
  } else {cat("No coefficients\n")
  }
  cat("\n")
  invisible(x)
}

#' @importFrom stats var
#' @export
summary.el_lm <- function(object, ...) {
  z <- object
  if (!inherits(object, "el_lm")) stop("invalid 'el_lm' object")
  if (is.null(z$terms)) stop("invalid 'el_lm' object:  no 'terms' component")
  p <- z$df
  rdf <- z$df.residual
  n <- p + rdf
  if (is.na(z$df.residual) || n - p != z$df.residual)
    warning("residual degrees of freedom in object suggest this is not an \"el_lm\" fit")
  r <- z$residuals
  f <- z$fitted.values
  mss <- if (attr(z$terms, "intercept")) sum((f - mean(f))^2) else sum(f^2)
  rss <- sum(r^2)
  resvar <- rss/rdf
  if (is.finite(resvar) && resvar < (mean(f)^2 + var(c(f))) * 1e-30)
    warning("essentially perfect fit: summary may be unreliable")
  p1 <- 1L:p
  ans <- z[c("call", "terms")]
  ans$residuals <- r
  ans$coefficients <- cbind(estimate = z$coefficients,
                            `chisq-value` = z$optim$par.tests$statistic,
                            `p-value` = z$optim$par.tests$p.value)
  ans$aliased <- is.na(z$coefficients)
  ans$df <- c(p, rdf, p)
  if (p != attr(z$terms, "intercept")) {
    df.int <- if (attr(z$terms, "intercept")) 1L else 0L
    ans$r.squared <- mss/(mss + rss)
    ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - df.int)/rdf)
    ans$chisq.statistic <- c(value = -2 * z$optim$logLR, df = p - df.int)
  } else
    ans$r.squared <- ans$adj.r.squared <- 0
  if (!is.null(z$na.action)) ans$na.action <- z$na.action
  class(ans) <- "summary.el_lm"
  ans
}

#' @importFrom stats naprint pchisq quantile
#' @export
print.summary.el_lm <- function(
  x, digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  resid <- x$residuals
  df <- x$df
  rdf <- df[2L]
  if (rdf > 5L) {
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    rq <- if (length(dim(resid)) == 2L)
      structure(apply(t(resid), 1L, quantile),
                dimnames = list(nam, dimnames(resid)[[2L]]))
    else {
      zz <- zapsmall(quantile(resid), digits + 1L)
      structure(zz, names = nam)
    }
    print(rq, digits = digits, ...)
  }
  else if (rdf > 0L) {
    print(resid, digits = digits, ...)
  }
  else {
    cat("ALL", df[1L], "residuals are 0: no residual degrees of freedom!")
    cat("\n")
  }
  if (length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n")
  }
  else {
    if (nsingular <- df[3L] - df[1L])
      cat("\nCoefficients: (", nsingular,
          " not defined because of singularities)\n", sep = "")
    else cat("\nCoefficients:\n")
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
    cat("Multiple R-squared: ", formatC(x$r.squared, digits = digits))
    cat(",\tAdjusted R-squared: ", formatC(x$adj.r.squared, digits = digits),
        "\nChisq-statistic:", formatC(x$chisq.statistic[1L], digits = digits),
        "on", x$chisq.statistic[2L], "DF, p-value:",
        format.pval(pchisq(x$chisq.statistic[1L], x$chisq.statistic[2L],
                           lower.tail = FALSE),
                    digits = digits))
    cat("\n")
  }
  cat("\n")
  invisible(x)
}
