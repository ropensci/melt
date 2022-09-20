#' @describeIn summary Summarizes the results of the overall model test and the
#'   significance tests for coefficients.
setMethod("summary", "LM", function(object, ...) {
  z <- object
  p <- getNumPar(z)
  if (p == 0L) {
    return(new("SummaryLM",
      statistic = chisq(z), df = getDF(z), convergence = conv(z),
      sigTests = matrix(NA_real_, 0L, 3L,
        dimnames = list(NULL, c("Estimate", "Chisq", "Pr(>Chisq)"))
      ),
      weighted = !is.null(weights(z)), intercept = z@misc$intercept,
      na.action = z@misc$na.action, call = z@call, terms = z@terms,
      aliased = is.na(coef(z))
    ))
  }
  new("SummaryLM",
    statistic = chisq(z), df = getDF(z), convergence = conv(z),
    sigTests = cbind(
      Estimate = coef(z),
      Chisq = sigTests(z)$statistic,
      `Pr(>Chisq)` = pchisq(sigTests(z)$statistic,
        df = 1L,
        lower.tail = FALSE
      )
    ),
    weighted = !is.null(weights(z)), intercept = z@misc$intercept,
    na.action = z@misc$na.action, call = z@call, terms = z@terms,
    aliased = is.na(coef(z))
  )
})

#' @describeIn summary Summarizes the results of the overall model test and the
#'   significance tests for coefficients. The dispersion parameter is extracted
#'   for display.
setMethod("summary", "GLM", function(object, ...) {
  z <- object
  p <- getNumPar(z)
  if (p == 0L) {
    return(new("SummaryGLM",
      family = z@family, dispersion = z@dispersion,
      statistic = chisq(z), df = getDF(z), convergence = conv(z),
      sigTests = matrix(NA_real_, 0L, 3L,
        dimnames = list(NULL, c("Estimate", "Chisq", "Pr(>Chisq)"))
      ),
      weighted = !is.null(weights(z)), intercept = z@misc$intercept,
      na.action = z@misc$na.action, call = z@call, terms = z@terms,
      aliased = is.na(coef(z))
    ))
  }
  new("SummaryGLM",
    family = z@family, dispersion = z@dispersion,
    statistic = chisq(z), df = getDF(z), convergence = conv(z),
    sigTests = cbind(
      Estimate = coef(z),
      Chisq = sigTests(z)$statistic,
      `Pr(>Chisq)` = pchisq(sigTests(z)$statistic,
        df = 1L,
        lower.tail = FALSE
      )
    ),
    weighted = !is.null(weights(z)), intercept = z@misc$intercept,
    na.action = z@misc$na.action, call = z@call, terms = z@terms,
    aliased = is.na(coef(z))
  )
})

#' @describeIn summary Summarizes the results of the overall model test and the
#'   significance tests for coefficients. The estimated dispersion parameter is
#'   extracted for display.
setMethod("summary", "QGLM", function(object, ...) {
  z <- object
  p <- getNumPar(z) - 1L
  if (p == 0L) {
    return(new("SummaryQGLM",
      family = z@family, dispersion = z@dispersion,
      statistic = chisq(z), df = getDF(z), convergence = conv(z),
      sigTests = matrix(NA_real_, 0L, 3L,
        dimnames = list(NULL, c("Estimate", "Chisq", "Pr(>Chisq)"))
      ),
      weighted = !is.null(weights(z)), intercept = z@misc$intercept,
      na.action = z@misc$na.action, call = z@call, terms = z@terms,
      aliased = is.na(coef(z))
    ))
  }
  new("SummaryQGLM",
    family = z@family, dispersion = z@dispersion,
    statistic = chisq(z), df = getDF(z), convergence = conv(z),
    sigTests = cbind(
      Estimate = coef(z),
      Chisq = sigTests(z)$statistic,
      `Pr(>Chisq)` = pchisq(sigTests(z)$statistic,
        df = 1L,
        lower.tail = FALSE
      )
    ),
    weighted = !is.null(weights(z)), intercept = z@misc$intercept,
    na.action = z@misc$na.action, call = z@call, terms = z@terms,
    aliased = is.na(coef(z))
  )
})
