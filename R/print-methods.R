#' @rdname print
setMethod("print", "EL", function(x,
                                  digits = max(3L, getOption("digits") - 3L),
                                  ...) {
  cat("\nEmpirical Likelihood:", getMethodEL(x), "\n\n")
  if (length(coef(x)) != 0L) {
    cat("Maximum EL estimates:\n")
    print.default(coef(x), digits = digits, ...)
  }
  cat("\n")

  out <- character()
  if (length(x@statistic) != 0L) {
    out <- c(
      out, paste("Chisq:", format.default(x@statistic, digits = digits)),
      paste("df:", x@df),
      paste("Pr(>Chisq):", format.pval(x@pval, digits = digits))
    )
  } else {
    out <- c("Empty model")
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n\n")

  if (length(x@statistic) != 0L) {
    cat(
      "\nEL evaluation:",
      if (x@optim$convergence) "converged" else "not converged", "\n"
    )
  }
  cat("\n")
  invisible(x)
})
setMethod("show", "EL", function(object) print(object))


#' @rdname print
#' @importFrom stats naprint pchisq printCoefmat
setMethod(
  "print", "SummaryLM", function(x,
                                 digits = max(3L, getOption("digits") - 3L),
                                 signif.stars = getOption("show.signif.stars"),
                                 ...) {
    cat("\nCall:\n", paste(deparse(x@call), sep = "\n", collapse = "\n"),
      "\n",
      sep = ""
    )
    if (length(x@aliased) == 0L) {
      cat("\nNo Coefficients\n")
    } else {
      cat("\nCoefficients:\n")
      coefs <- x@parMatrix
      if (any(aliased <- x@aliased)) {
        cn <- names(aliased)
        coefs <-
          matrix(NA, length(aliased), 3L, dimnames = list(cn, colnames(coefs)))
        coefs[!aliased, ] <- x@parMatrix
      }
      printCoefmat(coefs,
        digits = digits, signif.stars = signif.stars,
        P.values = TRUE, has.Pvalue = TRUE, na.print = "NA", ...
      )
    }
    cat("\n")

    mess <- naprint(x@na.action)
    if (isTRUE(nzchar(mess))) {
      cat("  (", mess, ")\n", sep = "", "\n")
    }

    out <- character()
    if (length(x@statistic) != 0L) {
      out <- c(
        out,
        paste("Chisq:", format(x@statistic, digits = digits)),
        paste("df:", x@df),
        paste("Pr(>Chisq):", format.pval(
          pchisq(x@statistic, x@df, lower.tail = FALSE),
          digits = digits
        ))
      )
    } else {
      out <- c("Empty model")
    }
    cat(strwrap(paste(out, collapse = ", ")), "\n\n")

    if (length(x@statistic) != 0L) {
      cat(
        "Constrained EL:",
        if (x@convergence) "converged" else "not converged", "\n\n"
      )
    }
    invisible(x)
  }
)
setMethod("show", "SummaryLM", function(object) print(object))


#' @rdname print
setMethod("print", "logLikEL", function(x, digits = getOption("digits"), ...) {
  cat("'Empirical log Lik.' ", paste(format(c(x@logLik), digits = digits),
    collapse = ", "
  ),
  " (df=", format(x@df), ")\n",
  sep = ""
  )
  invisible(x)
})
setMethod("show", "logLikEL", function(object) print(object))


#' @rdname print
setMethod("print", "ELT", function(x, digits = getOption("digits"), ...) {
  cat("\nEmpirical Likelihood Test\n\n")
  method <- switch(x@calibrate,
    "chisq" = "Chi-square",
    "boot" = "Bootstrap",
    "f" = "F"
  )
  out <- character()
  out <- c(
    out, paste("Significance level:", x@alpha),
    paste("Calibration:", method)
  )
  cat(strwrap(paste(out, collapse = ", ")), "\n\n")
  out2 <- character()
  out2 <- c(
    out2, paste("Statistic:", format.default(x@statistic, digits = digits)),
    paste("Critical value:", format.default(x@cv, digits = digits))
  )
  cat(strwrap(paste(out2, collapse = ", ")), "\n\n")
  cat("p-value:", format.pval(x@pval, digits = digits), "\n\n")
  invisible(x)
})
setMethod("show", "ELT", function(object) print(object))









#' @param x An object of class \linkS4class{ELMT}.
#' @rdname ELMT-class
#' @exportMethod print
setMethod("print", "ELMT", function(x) {print(x@cv)})






