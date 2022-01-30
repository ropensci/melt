#' @export
print.el_aov <- function(x, ...) {
  stopifnot(inherits(x, "melt"))
  cat("Call:\n")
  dput(x$call, control = NULL)
  cat("\nminimizer:\n")
  cat(format(round(x$optim$par, 4), scientific = F))
  cat("\n\n")
  cat("statistic:\n")
  cat(format(round(x$optim$n2logLR, 4), scientific = F))
  cat("\n\n")
}

#' @importFrom stats coef
#' @export
print.el_lm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  stopifnot(inherits(x, "melt"))
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
  p <- z$rank
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

    ans$chisq.statistic <- c(value = z$optim$n2logLR, df = p)

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

#' @export
print.el_test <- function(x, ...) {
  stopifnot(inherits(x, "melt"))
  cat("\n")
  cat("Empirical Likelihood Hypothesis Testing\n\n")
  cat("minimizer:\n")
  cat(format(round(x$optim$par, 4), scientific = F))
  cat("\n\n")
  cat("statistic:\n")
  cat(format(round(x$optim$n2logLR, 4), scientific = F))
  cat("\n\n")
}

#' @importFrom stats printCoefmat
#' @export
print.pairwise <- function(x, ...) {
  stopifnot(inherits(x, "melt"))
  cat("\n")
  cat("Empirical Likelihood Multiple Hypothesis Testing\n\n")
  # set row names
  if (is.null(x$control)) {
    cat("Test: all pairwise comparisons\n\n")
    rname <- vector("character", length = 0)
    for (i in 1L:(length(x$trt) - 1L)) {
      for (j in (i + 1L):length(x$trt)) {
        rname <- c(rname, paste(x$trt[i], "-", x$trt[j]))
      }
    }
  } else {
    cat("Test: comparisons with control\n\n")
    diff <- setdiff(x$trt, x$control)
    rname <- vector("character", length = length(diff))
    for (i in 1L:length(diff)) {
      rname[i] <- paste(diff[i], "-", x$control)
    }
  }
  out <- data.frame(row.names = rname, estimate  = x$estimate,
                    statistic = x$statistic, lwr.ci = x$lower,
                    upr.ci = x$upper,
                    p.adj = round(x$p.adj, 4))
  printCoefmat(out, digits = min(4L, getOption("digits")), cs.ind = c(1, 3, 4),
               tst.ind = 2L, dig.tst = min(3L, getOption("digits")),
               P.values = T, has.Pvalue = T, eps.Pvalue = 1e-03)
  cat("\n")
  cat(paste(c("k", "level", "method", "cutoff"),
            c(x$k, x$level, x$method, round(x$cutoff, 4)),
            collapse = ", ", sep = ": "))
  cat("\n\n")
}

# check_hypotheses <- function(hypotheses, rhs, trt) {
#   # hypotheses must be either character vector or list of numeric matrices
#   match.arg(class(hypotheses), c("character", "list"))
#   # test for pairwise comparisons
#   test <- NULL
#   if (is.character(hypotheses)) {
#     if (length(hypotheses) == 1 && hypotheses == "pairwise") {
#       test <- 0
#       # all pairwise
#     } else if (length(hypotheses) == 2 && hypotheses[1] == "pairwise") {
#       if (hypotheses[2] %in% trt) {
#         test <- which(trt == hypotheses[2])
#       } else {
#         stop(paste("there is no treatment named", hypotheses[2]))
#       }
#     } else {
#       stop("invalid hypotheses specified")
#     }
#     return(test)
#   }
#
#   # check validity of hypotheses and rhs
#   if (!all(sapply(hypotheses, function(x) is.numeric(x) && is.matrix(x)))) {
#     stop("hypotheses must consists of numeric matrices")
#   } else if (!all(sapply(hypotheses, function(x) {
#     nrow(x) <= ncol(x) && ncol(x) == length(trt)}))) {
#     stop("invalid matrix dimensions for hypotheses")
#   }
#
#   if (!is.null(rhs)) {
#     if (!is.list(rhs)) {
#       stop("invalid rhs specified")
#     } else if (length(hypotheses) != length(rhs)) {
#       stop("number of hypotheses and rhs do not match")
#     } else if (!all(sapply(rhs, function(x) is.numeric(x) && is.vector(x)))) {
#       stop("rhs must consists of numeric vectors")
#     }
#   }
#   test
# }

# elmulttest_rownames <- function(m, trt, test) {
#   # if not pairwise comparison, return general hypotheses names
#   if (is.null(test)) {
#     return(paste0("H_", 1:m))
#   }
#
#   # pairwise comparisons
#   if (length(test) == 2) {
#     # comparisons with control
#     diff <- setdiff(trt, test[2])
#     row_names <- vector("character", length = length(diff))
#     for (i in 1:length(diff)) {
#       row_names[i] <- paste(diff[i], "-", test[2])
#     }
#     return(row_names)
#   } else {
#     # all pairwise comparisons
#     row_names <- vector("character", length = 0)
#     for (i in 1:(length(trt) - 1)) {
#       for (j in (i + 1):length(trt)) {
#         row_names <- c(row_names, paste(trt[i], "-", trt[j]))
#       }
#     }
#     row_names
#   }
# }

