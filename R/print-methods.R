#' @rdname print-method
setMethod(
  "print", "EL",
  function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("\nEmpirical Likelihood:", x@optim$method, "\n\n")
    if (length(x@coefficients) != 0L) {
      cat("Maximum EL estimates:\n")
      print.default(x@coefficients, digits = digits, ...)
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
  }
)
setMethod("show", "EL", function(object) print(object))

#' @rdname print-method
#' @importFrom stats naprint pchisq
setMethod(
  "print", "SummaryLM",
  function(x, digits = max(3L, getOption("digits") - 3L),
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

    if (nzchar(mess <- naprint(x@na.action))) {
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

#' @rdname print-method
setMethod(
  "print", "logLikEL",
  function(x, digits = getOption("digits"), ...) {
    cat("'Empirical log Lik.' ", paste(format(c(x@logLik), digits = digits),
      collapse = ", "
    ),
    " (df=", format(x@df), ")\n",
    sep = ""
    )
    invisible(x)
  }
)
setMethod("show", "logLikEL", function(object) print(object))
