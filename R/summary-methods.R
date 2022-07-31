#' @describeIn summary Summarizes the results of the overall model test and the
#'   significance tests for coefficients.
#' @importFrom stats pchisq
#' @srrstats {RE4.18} `summary` method is applicable to a `LM` object returned
#'   by `el_lm()` or `el_glm()`.
setMethod("summary", "LM", function(object, ...) {
  z <- object
  p <- getNumPar(z)
  if (p == 0L) {
    return(new("SummaryLM",
      statistic = chisq(z), df = getDF(z), convergence = conv(z),
      sigTests = matrix(NA_real_, 0L, 3L,
        dimnames = list(NULL, c("Estimate", "Chisq", "Pr(>Chis)"))
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
