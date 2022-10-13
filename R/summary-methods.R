#' @describeIn summary Summarizes the test results of the specified parameters.
setMethod("summary", "EL", function(object, ...) {
  z <- object
  new("SummaryEL",
    logl = logL(z), loglr = logLR(z), statistic = chisq(z), df = getDF(z),
    pval = pVal(z), nobs = nobs(z), npar = getNumPar(z), coefficients = coef(z),
    method = getMethodEL(z), weighted = !is.null(weights(z)),
    convergence = conv(z), par = getOptim(z)$par, lambda = getOptim(z)$lambda
  )
})

#' @describeIn summary Summarizes the results of the overall model test and the
#'   significance tests for coefficients.
setMethod("summary", "LM", function(object, ...) {
  z <- object
  p <- getNumPar(z)
  est <- coef(z)
  if (p == 0L) {
    return(new("SummaryLM",
      sigTests = matrix(NA_real_, 0L, 3L,
        dimnames = list(NULL, c("Estimate", "Chisq", "Pr(>Chisq)"))
      ),
      intercept = z@misc$intercept, na.action = z@misc$na.action, call = z@call,
      terms = z@terms, aliased = is.na(est), weighted = !is.null(weights(z)),
      convergence = conv(z), logl = logL(z), loglr = logLR(z),
      statistic = chisq(z), df = getDF(z), pval = pVal(z), nobs = nobs(z),
      npar = p, coefficients = est, method = getMethodEL(z)
    ))
  }
  new("SummaryLM",
    sigTests = cbind(
      Estimate = est,
      Chisq = sigTests(z)$statistic,
      `Pr(>Chisq)` = pchisq(sigTests(z)$statistic,
        df = 1L,
        lower.tail = FALSE
      )
    ),
    intercept = z@misc$intercept, na.action = z@misc$na.action, call = z@call,
    terms = z@terms, aliased = is.na(est), par = getOptim(z)$par,
    lambda = getOptim(z)$lambda, weighted = !is.null(weights(z)),
    convergence = conv(z), logl = logL(z), loglr = logLR(z),
    statistic = chisq(z), df = getDF(z), pval = pVal(z), nobs = nobs(z),
    npar = p, coefficients = est, method = getMethodEL(z)
  )
})

#' @describeIn summary Summarizes the results of the overall model test and the
#'   significance tests for coefficients. The dispersion parameter is extracted
#'   for display.
setMethod("summary", "GLM", function(object, ...) {
  z <- object
  p <- getNumPar(z)
  est <- coef(z)
  if (p == 0L) {
    return(new("SummaryGLM",
      family = z@family, dispersion = z@dispersion,
      sigTests = matrix(NA_real_, 0L, 3L,
        dimnames = list(NULL, c("Estimate", "Chisq", "Pr(>Chisq)"))
      ),
      intercept = z@misc$intercept, na.action = z@misc$na.action, call = z@call,
      terms = z@terms, aliased = is.na(est), weighted = !is.null(weights(z)),
      convergence = conv(z), logl = logL(z), loglr = logLR(z),
      statistic = chisq(z), df = getDF(z), pval = pVal(z), nobs = nobs(z),
      npar = p, coefficients = est, method = getMethodEL(z)
    ))
  }
  new("SummaryGLM",
    family = z@family, dispersion = z@dispersion,
    sigTests = cbind(
      Estimate = est,
      Chisq = sigTests(z)$statistic,
      `Pr(>Chisq)` = pchisq(sigTests(z)$statistic,
        df = 1L,
        lower.tail = FALSE
      )
    ),
    intercept = z@misc$intercept, na.action = z@misc$na.action, call = z@call,
    terms = z@terms, aliased = is.na(est), par = getOptim(z)$par,
    lambda = getOptim(z)$lambda, weighted = !is.null(weights(z)),
    convergence = conv(z), logl = logL(z), loglr = logLR(z),
    statistic = chisq(z), df = getDF(z), pval = pVal(z), nobs = nobs(z),
    npar = p, coefficients = est, method = getMethodEL(z)
  )
})

#' @describeIn summary Summarizes the results of the overall model test and the
#'   significance tests for coefficients. The estimated dispersion parameter is
#'   extracted for display.
setMethod("summary", "QGLM", function(object, ...) {
  z <- object
  p <- getNumPar(z) - 1L
  est <- coef(z)
  if (p == 0L) {
    return(new("SummaryQGLM",
      family = z@family, dispersion = z@dispersion,
      sigTests = matrix(NA_real_, 0L, 3L,
        dimnames = list(NULL, c("Estimate", "Chisq", "Pr(>Chisq)"))
      ),
      intercept = z@misc$intercept, na.action = z@misc$na.action, call = z@call,
      terms = z@terms, aliased = is.na(est), weighted = !is.null(weights(z)),
      convergence = conv(z), logl = logL(z), loglr = logLR(z),
      statistic = chisq(z), df = getDF(z), pval = pVal(z), nobs = nobs(z),
      npar = p, coefficients = est, method = getMethodEL(z)
    ))
  }
  new("SummaryQGLM",
    family = z@family, dispersion = z@dispersion,
    sigTests = cbind(
      Estimate = est,
      Chisq = sigTests(z)$statistic,
      `Pr(>Chisq)` = pchisq(sigTests(z)$statistic,
        df = 1L,
        lower.tail = FALSE
      )
    ),
    intercept = z@misc$intercept, na.action = z@misc$na.action, call = z@call,
    terms = z@terms, aliased = is.na(est), par = getOptim(z)$par,
    lambda = getOptim(z)$lambda, weighted = !is.null(weights(z)),
    convergence = conv(z), logl = logL(z), loglr = logLR(z),
    statistic = chisq(z), df = getDF(z), pval = pVal(z), nobs = nobs(z),
    npar = p, coefficients = est, method = getMethodEL(z)
  )
})

#' @describeIn summary Summarizes the.
setMethod("summary", "ELT", function(object, ...) {
  new("SummaryELT")
})

#' @describeIn summary Summarizes the.
setMethod("summary", "ELMT", function(object, ...) {
  new("SummaryELMT")
})
