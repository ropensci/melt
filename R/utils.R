#' @importFrom stats pf
tt_ <- function(calibrate, statistic, p, alpha, object) {
  switch(calibrate,
    "chisq" = {
      c(cv = 1, pval = pchisq(q = statistic, df = p, lower.tail = FALSE))
    },
    "boot" = {
      stop("not yet")
      c(cv = 2, pval = pchisq(q = statistic, df = p, lower.tail = FALSE))
    },
    "f" = {
      n <- nrow(getDataMatrix(object))
      c(cv = 3, pval = pf(
        q = statistic * (n - p) / (p * (n - 1)), df1 = p,
        df2 = n - p, lower.tail = FALSE
      ))
    }
  )
}
