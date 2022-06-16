#' @importFrom stats pf
calibrate_pval_ <- function(calibrate, statistic, p, object) {
  switch(calibrate,
    "chisq" = {
      pchisq(q = statistic, df = p, lower.tail = FALSE)
    },
    "boot" = {
      stop("not yet")
    },
    "f" = {
      n <- nrow(getDataMatrix(object))
      pf(
        q = statistic * (n - p) / (p * (n - 1)), df1 = p, df2 = n - p,
        lower.tail = FALSE
      )
    }
  )
}
