#' @rdname print
setMethod("print", "EL", function(x,
                                  digits = max(3L, getOption("digits") - 3L),
                                  ...) {
  if (is.null(weights(x))) {
    cat("\n\tEmpirical Likelihood\n")
  } else {
    cat("\n\tWeighted Empirical Likelihood\n")
  }
  cat("\nModel:", getMethodEL(x), "\n\nMaximum EL estimates:\n")
  print.default(coef(x), digits = digits, ...)
  cat("\n")
  out <- c(
    paste("Chisq:", format.default(chisq(x), digits = digits, ...)),
    paste("df:", getDF(x)),
    paste("Pr(>Chisq):", format.pval(pVal(x), digits = digits, ...))
  )
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n\n")
  cat(show_convergence(x), "\n\n")
  invisible(x)
})
setMethod("show", "EL", function(object) print(object))

#' @rdname print
setMethod("print", "ELMT", function(x,
                                    digits = max(3L, getOption("digits") - 3L),
                                    ...) {
  method <- switch(x@calibrate,
    "mvchisq" = "Multivariate chi-square"
  )
  cat("\n\tEmpirical Likelihood Multiple Tests\n\n")
  cat("Overall significance level:", x@alpha, "\n\n")
  cat("Calibration:", method, "\n\n")
  cat("Hypotheses:\n")
  est <- getEstimates(x)
  if (isTRUE(all(vapply(est, FUN = length, FUN.VALUE = integer(1L)) == 1L))) {
    out <- cbind(Estimate = unlist(est), Chisq = chisq(x), Df = getDF(x))
    rownames(out) <- describe_hypothesis(x@rhs, x@lhs, colnames(x@lhs), digits)
    ind <- 2L
  } else {
    out <- cbind(Chisq = chisq(x), Df = getDF(x))
    rownames(out) <- seq_along(pVal(x))
    ind <- 1L
  }
  printCoefmat(out, digits = digits, tst.ind = ind)
  cat("\n")
  invisible(x)
})
setMethod("show", "ELMT", function(object) print(object))

#' @rdname print
setMethod("print", "ELT", function(x,
                                   digits = max(3L, getOption("digits") - 3L),
                                   ...) {
  cat("\n\tEmpirical Likelihood Test\n\n")
  cat("Hypothesis:", describe_hypothesis(
    x@rhs, x@lhs, names(getOptim(x)$par)[seq_len(ncol(x@lhs))], digits
  ), sep = "\n")
  method <- switch(x@calibrate,
    "ael" = "Adjusted EL",
    "boot" = "Bootstrap",
    "chisq" = "Chi-square",
    "f" = "F"
  )
  cat(paste0(
    "\nSignificance level: ", format.default(x@alpha, digits = digits, ...),
    ", Calibration: ", method
  ), "\n")
  cat(
    paste0(
      "\nStatistic: ", format.default(chisq(x), digits = digits, ...),
      ", Critical value: ", format.default(critVal(x), digits = digits, ...),
      "\np-value: ", format.pval(pVal(x), digits = digits)
    ), "\n"
  )
  cat(show_convergence(x), "\n\n")
  invisible(x)
})
setMethod("show", "ELT", function(object) print(object))

#' @rdname print
setMethod("print", "LM", function(x,
                                  digits = max(3L, getOption("digits") - 3L),
                                  ...) {
  if (is.null(weights(x))) {
    cat("\n\tEmpirical Likelihood\n\n")
  } else {
    cat("\n\tWeighted Empirical Likelihood\n\n")
  }
  method <- getMethodEL(x)
  if (isFALSE(is.na(method))) {
    if (is(x, "GLM")) {
      model <- unlist(strsplit(method, split = "_", fixed = TRUE))
      cat("Model: glm (", model[1L], " family with ", model[2L], " link)",
        "\n\n",
        sep = ""
      )
    } else {
      cat("Model: lm", "\n\n")
    }
  }
  if (length(coef(x)) != 0L) {
    cat("Maximum EL estimates:\n")
    print.default(coef(x), digits = digits, ...)
    cat("\n")
  }
  if (length(chisq(x)) != 0L) {
    out <- c(
      paste("Chisq:", format.default(chisq(x), digits = digits, ...)),
      paste("df:", getDF(x)),
      paste("Pr(>Chisq):", format.pval(pVal(x), digits = digits, ...))
    )
  } else {
    out <- c("Empty model")
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n\n")
  if (length(chisq(x)) != 0L) {
    cat(show_convergence(x), "\n")
  }
  cat("\n")
  invisible(x)
})
setMethod("show", "LM", function(object) print(object))

#' @rdname print
setMethod(
  "print", "SummaryEL",
  function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    if (x@weighted) {
      cat("\n\tWeighted Empirical Likelihood\n")
    } else {
      cat("\n\tEmpirical Likelihood\n")
    }
    cat(
      "\nModel:", getMethodEL(x),
      "\n\nNumber of observations:", nobs(x),
      "\nNumber of parameters:", getNumPar(x),
      "\n\nParameter values under the null hypothesis:\n"
    )
    print.default(getOptim(x)$par, digits = digits, ...)
    cat("\nLagrange multipliers:\n")
    print.default(getOptim(x)$lambda, digits = digits, ...)
    cat("\nMaximum EL estimates:\n")
    print.default(coef(x), digits = digits, ...)
    cat(
      paste0(
        "\nlogL: ", format.default(logL(x), digits = digits, ...),
        ", logLR: ", format.default(logLR(x), digits = digits, ...)
      ), "\n"
    )
    out <- c(
      paste("Chisq:", format.default(chisq(x), digits = digits), ...),
      paste("df:", getDF(x)),
      paste("Pr(>Chisq):", format.pval(pVal(x), digits = digits), ...)
    )
    cat(strwrap(paste(out, collapse = ", ")), sep = "\n\n")
    cat(show_convergence(x), "\n\n")
    invisible(x)
  }
)
setMethod("show", "SummaryEL", function(object) print(object))

#' @rdname print
setMethod(
  "print", "SummaryELMT",
  function(x,
           digits = max(3L, getOption("digits") - 3L),
           signif.stars = getOption("show.signif.stars"),
           ...) {
    method <- switch(x@calibrate,
      "mvchisq" = "Multivariate chi-square"
    )
    cat("\n\tEmpirical Likelihood Multiple Tests\n\n")
    cat("Overall significance level:", x@alpha, "\n\n")
    cat("Calibration:", method, "\n\n")
    cat("Hypotheses:\n")
    est <- getEstimates(x)
    if (isTRUE(all(vapply(est, FUN = length, FUN.VALUE = integer(1L)) == 1L))) {
      out <- cbind(
        Estimate = unlist(est), Chisq = chisq(x), Df = getDF(x), p.adj = pVal(x)
      )
      rownames(out) <- describe_hypothesis(
        x@rhs, x@lhs, colnames(x@lhs), digits
      )
      ind <- 2L
    } else {
      out <- cbind(Chisq = chisq(x), Df = getDF(x), p.adj = pVal(x))
      rownames(out) <- seq_along(pVal(x))
      ind <- 1L
    }
    printCoefmat(out,
      digits = digits, signif.stars = signif.stars, tst.ind = ind,
      P.values = TRUE, has.Pvalue = TRUE, eps.Pvalue = 1e-03
    )
    cat("\n")
    cat(paste(
      "Common critical value:", format.default(critVal(x), digits = digits, ...)
    ), "\n\n")
    invisible(x)
  }
)
setMethod("show", "SummaryELMT", function(object) print(object))

#' @rdname print
setMethod(
  "print", "SummaryELT",
  function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("\n\tEmpirical Likelihood Test\n\n")
    cat("Hypothesis:", describe_hypothesis(
      x@rhs, x@lhs, names(getOptim(x)$par)[seq_len(ncol(x@lhs))], digits
    ), sep = "\n")
    method <- switch(x@calibrate,
      "ael" = "Adjusted EL",
      "boot" = "Bootstrap",
      "chisq" = "Chi-square",
      "f" = "F"
    )
    cat(paste0(
      "\nSignificance level: ", format.default(x@alpha, digits = digits, ...),
      ", Calibration: ", method
    ), "\n")
    cat("\nParameter values under the null hypothesis:\n")
    print.default(getOptim(x)$par, digits = digits, ...)
    cat("\nLagrange multipliers:\n")
    print.default(getOptim(x)$lambda, digits = digits, ...)
    cat(
      paste0(
        "\nlogL: ", format.default(logL(x), digits = digits, ...),
        ", logLR: ", format.default(logLR(x), digits = digits, ...),
        "\nStatistic: ", format.default(chisq(x), digits = digits, ...),
        ", Critical value: ", format.default(critVal(x), digits = digits, ...),
        "\np-value: ", format.pval(pVal(x), digits = digits)
      ), "\n"
    )
    cat(show_convergence(x), "\n\n")
    invisible(x)
  }
)
setMethod("show", "SummaryELT", function(object) print(object))

#' @rdname print
setMethod(
  "print", "SummaryGLM", function(x,
                                  digits = max(3L, getOption("digits") - 3L),
                                  signif.stars = getOption("show.signif.stars"),
                                  ...) {
    if (x@weighted) {
      cat("\n\tWeighted Empirical Likelihood\n")
    } else {
      cat("\n\tEmpirical Likelihood\n")
    }
    cat("\nModel: glm (", x@family$family, " family with ", x@family$link,
      " link)", "\n",
      sep = ""
    )
    cat("\nCall:\n", paste(deparse(x@call, ...), sep = "\n", collapse = "\n"),
      "\n",
      sep = ""
    )
    if (length(chisq(x)) != 0L) {
      cat(
        "\nNumber of observations:", nobs(x),
        "\nNumber of parameters:", getNumPar(x), "\n"
      )
      cat("\nParameter values under the null hypothesis:\n")
      print.default(getOptim(x)$par, digits = digits, ...)
      cat("\nLagrange multipliers:\n")
      print.default(getOptim(x)$lambda, digits = digits, ...)
      cat("\nMaximum EL estimates:\n")
      print.default(coef(x)[, 1L], digits = digits, ...)
      cat(paste(
        "\nlogL:", format.default(logL(x), digits = digits, ...),
        ", logLR:", format.default(logLR(x), digits = digits, ...)
      ), "\n")
      out <- c(
        paste("Chisq:", format(chisq(x), digits = digits)),
        paste("df:", getDF(x)),
        paste("Pr(>Chisq):", format.pval(pVal(x), digits = digits))
      )
      cat(strwrap(paste(out, collapse = ", ")), sep = "\n\n")
      cat(show_convergence(x), "\n")
    } else {
      cat("\nEmpty model\n")
    }
    aliased <- x@aliased
    if (length(aliased) == 0L) {
      cat("\nNo Coefficients\n")
    } else {
      cat("\nCoefficients:\n")
      coefs <- coef(x)
      if (any(aliased)) {
        cn <- names(aliased)
        coefs <-
          matrix(NA, length(aliased), 3L, dimnames = list(cn, colnames(coefs)))
        coefs[!aliased, ] <- coef(x)
      }
      printCoefmat(coefs,
        digits = digits, signif.stars = signif.stars,
        P.values = TRUE, has.Pvalue = TRUE, na.print = "NA", ...
      )
      cat("\nDispersion for ", x@family$family, " family: ",
        format(x@dispersion), "\n",
        sep = ""
      )
    }
    mess <- naprint(x@na.action)
    if (isTRUE(nzchar(mess))) {
      cat("  (", mess, ")\n", sep = "")
    }
    cat("\n")
    invisible(x)
  }
)
setMethod("show", "SummaryGLM", function(object) print(object))

#' @rdname print
setMethod(
  "print", "SummaryLM", function(x,
                                 digits = max(3L, getOption("digits") - 3L),
                                 signif.stars = getOption("show.signif.stars"),
                                 ...) {
    if (x@weighted) {
      cat("\n\tWeighted Empirical Likelihood\n")
    } else {
      cat("\n\tEmpirical Likelihood\n")
    }
    cat("\nModel:", getMethodEL(x), "\n")
    cat("\nCall:\n", paste(deparse(x@call, ...), sep = "\n", collapse = "\n"),
      "\n",
      sep = ""
    )
    if (length(chisq(x)) != 0L) {
      cat(
        "\nNumber of observations:", nobs(x),
        "\nNumber of parameters:", getNumPar(x), "\n"
      )
      cat("\nParameter values under the null hypothesis:\n")
      print.default(getOptim(x)$par, digits = digits, ...)
      cat("\nLagrange multipliers:\n")
      print.default(getOptim(x)$lambda, digits = digits, ...)
      cat("\nMaximum EL estimates:\n")
      print.default(coef(x)[, 1L], digits = digits, ...)
      cat(paste(
        "\nlogL:", format.default(logL(x), digits = digits, ...),
        ", logLR:", format.default(logLR(x), digits = digits, ...)
      ), "\n")
      out <- c(
        paste("Chisq:", format(chisq(x), digits = digits)),
        paste("df:", getDF(x)),
        paste("Pr(>Chisq):", format.pval(pVal(x), digits = digits))
      )
      cat(strwrap(paste(out, collapse = ", ")), sep = "\n\n")
      cat(show_convergence(x), "\n")
    } else {
      cat("\nEmpty model\n")
    }
    aliased <- x@aliased
    if (length(aliased) == 0L) {
      cat("\nNo Coefficients\n")
    } else {
      cat("\nCoefficients:\n")
      coefs <- coef(x)
      if (any(aliased)) {
        cn <- names(aliased)
        coefs <-
          matrix(NA, length(aliased), 3L, dimnames = list(cn, colnames(coefs)))
        coefs[!aliased, ] <- coef(x)
      }
      printCoefmat(coefs,
        digits = digits, signif.stars = signif.stars,
        P.values = TRUE, has.Pvalue = TRUE, na.print = "NA", ...
      )
    }
    mess <- naprint(x@na.action)
    if (isTRUE(nzchar(mess))) {
      cat("  (", mess, ")\n", sep = "")
    }
    cat("\n")
    invisible(x)
  }
)
setMethod("show", "SummaryLM", function(object) print(object))
