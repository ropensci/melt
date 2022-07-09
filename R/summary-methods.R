#' @describeIn summary Summarizes the results of the overall test and the
#'   tests for each parameter.
#' @importFrom stats pchisq
#' @srrstats {RE4.18} `summary` method is applicable to a `LM` object returned
#'   by `el_lm()` or `el_glm()`.
setMethod("summary", "LM", function(object, ...) {
  z <- object
  p <- getNumPar(z)
  if (p == 0L) {
    return(new("SummaryLM",
      statistic = z@statistic, df = z@df, convergence = z@optim$convergence,
      parMatrix = matrix(NA_real_, 0L, 3L,
        dimnames = list(NULL, c("Estimate", "Chisq", "Pr(>Chis)"))
      ),
      weighted = length(z@weights) != 0L, na.action = z@misc$na.action,
      call = z@call, terms = z@terms,
      aliased = is.na(coef(z))
    ))
  }
  new("SummaryLM",
    statistic = z@statistic, df = z@df, convergence = z@optim$convergence,
    parMatrix = cbind(
      Estimate = coef(z),
      Chisq = z@parTests$statistic,
      `Pr(>Chisq)` = pchisq(z@parTests$statistic, df = 1L, lower.tail = FALSE)
    ),
    weighted = length(z@weights) != 0L, na.action = z@misc$na.action,
    call = z@call, terms = z@terms, aliased = is.na(coef(z))
  )
})
