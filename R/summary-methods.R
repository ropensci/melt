#' @describeIn summary Summarizes the test results of the specified parameters.
setMethod("summary", "EL", function(object, ...) {
  z <- object
  new("SummaryEL",
    optim = getOptim(z), logl = logL(z), loglr = logLR(z), statistic = chisq(z),
    df = getDF(z), pval = pVal(z), nobs = nobs(z), npar = getNumPar(z),
    weighted = !is.null(weights(z)), coefficients = coef(z),
    method = getMethodEL(z), control = object@control
  )
})

#' @describeIn summary Summarizes the multiple testing results.
setMethod("summary", "ELMT", function(object, ...) {
  z <- object
  new("SummaryELMT",
    estimates = getEstimates(z), statistic = chisq(z), df = getDF(z),
    pval = pVal(z), cv = z@cv, rhs = z@rhs, lhs = z@lhs, alpha = z@alpha,
    calibrate = z@calibrate
  )
})

#' @describeIn summary Summarizes the hypothesis test results.
setMethod("summary", "ELT", function(object, ...) {
  z <- object
  new("SummaryELT",
    optim = getOptim(z), logl = logL(z), loglr = logLR(z), statistic = chisq(z),
    df = getDF(z), pval = pVal(z), cv = z@cv, rhs = z@rhs, lhs = z@lhs,
    alpha = z@alpha, calibrate = z@calibrate, control = object@control
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
      coefficients = matrix(NA_real_, 0L, 3L,
        dimnames = list(NULL, c("Estimate", "Chisq", "Pr(>Chisq)"))
      ),
      intercept = z@misc$intercept, na.action = z@misc$na.action, call = z@call,
      terms = z@terms, aliased = is.na(est), optim = getOptim(z),
      logl = logL(z), loglr = logLR(z), statistic = chisq(z), df = getDF(z),
      pval = pVal(z), nobs = nobs(z), npar = p, weighted = !is.null(weights(z)),
      method = getMethodEL(z), control = object@control
    ))
  }
  new("SummaryGLM",
    family = z@family, dispersion = z@dispersion,
    coefficients = cbind(
      Estimate = est,
      Chisq = sigTests(z)$statistic,
      `Pr(>Chisq)` = pchisq(sigTests(z)$statistic,
        df = 1L,
        lower.tail = FALSE
      )
    ),
    intercept = z@misc$intercept, na.action = z@misc$na.action, call = z@call,
    terms = z@terms, aliased = is.na(est), optim = getOptim(z),
    logl = logL(z), loglr = logLR(z), statistic = chisq(z), df = getDF(z),
    pval = pVal(z), nobs = nobs(z), npar = p, weighted = !is.null(weights(z)),
    method = getMethodEL(z), control = object@control
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
      coefficients = matrix(NA_real_, 0L, 3L,
        dimnames = list(NULL, c("Estimate", "Chisq", "Pr(>Chisq)"))
      ),
      intercept = z@misc$intercept, na.action = z@misc$na.action, call = z@call,
      terms = z@terms, aliased = is.na(est), optim = getOptim(z),
      logl = logL(z), loglr = logLR(z), statistic = chisq(z), df = getDF(z),
      pval = pVal(z), nobs = nobs(z), npar = p, weighted = !is.null(weights(z)),
      method = getMethodEL(z), control = object@control
    ))
  }
  new("SummaryLM",
    coefficients = cbind(
      Estimate = est,
      Chisq = sigTests(z)$statistic,
      `Pr(>Chisq)` = pchisq(sigTests(z)$statistic,
        df = 1L,
        lower.tail = FALSE
      )
    ),
    intercept = z@misc$intercept, na.action = z@misc$na.action, call = z@call,
    terms = z@terms, aliased = is.na(est), optim = getOptim(z), logl = logL(z),
    loglr = logLR(z), statistic = chisq(z), df = getDF(z), pval = pVal(z),
    nobs = nobs(z), npar = p, weighted = !is.null(weights(z)),
    method = getMethodEL(z), control = object@control
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
      coefficients = matrix(NA_real_, 0L, 3L,
        dimnames = list(NULL, c("Estimate", "Chisq", "Pr(>Chisq)"))
      ),
      intercept = z@misc$intercept, na.action = z@misc$na.action, call = z@call,
      terms = z@terms, aliased = is.na(est), optim = getOptim(z),
      logl = logL(z), loglr = logLR(z), statistic = chisq(z), df = getDF(z),
      pval = pVal(z), nobs = nobs(z), npar = p, weighted = !is.null(weights(z)),
      method = getMethodEL(z), control = object@control
    ))
  }
  new("SummaryQGLM",
    family = z@family, dispersion = z@dispersion,
    coefficients = cbind(
      Estimate = est,
      Chisq = sigTests(z)$statistic,
      `Pr(>Chisq)` = pchisq(sigTests(z)$statistic,
        df = 1L,
        lower.tail = FALSE
      )
    ),
    intercept = z@misc$intercept, na.action = z@misc$na.action, call = z@call,
    terms = z@terms, aliased = is.na(est), optim = getOptim(z),
    logl = logL(z), loglr = logLR(z), statistic = chisq(z), df = getDF(z),
    pval = pVal(z), nobs = nobs(z), npar = p, weighted = !is.null(weights(z)),
    method = getMethodEL(z), control = object@control
  )
})
