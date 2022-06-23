#' @importFrom stats pf qf quantile
calibrate <- function(alpha, statistic, calibrate, p, par, object, control) {
  switch(calibrate,
    "chisq" = {
      c(
        cv = qchisq(p = 1 - alpha, df = p),
        pval = pchisq(q = statistic, df = p, lower.tail = FALSE)
      )
    },
    "boot" = {
      out <- boot_(
        control@B, control@seed, control@nthreads, getMethodEL(object),
        getDataMatrix(object), par, control@maxit_l, control@tol_l, control@th,
        getWeights(object)
      )
      c(
        cv = quantile(x = out, probs = 1 - alpha, names = FALSE),
        pval = mean.default(out > statistic)
      )
    },
    "f" = {
      n <- nrow(getDataMatrix(object))
      c(
        cv = qf(p = 1 - alpha, df1 = p, df2 = n - p) * (p * (n - 1)) / (n - p),
        pval = pf(
          q = statistic * (n - p) / (p * (n - 1)), df1 = p,
          df2 = n - p, lower.tail = FALSE
        )
      )
    }
  )
}
