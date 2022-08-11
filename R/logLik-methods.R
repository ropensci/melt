#' @rdname logLik
#' @srrstats {RE4.11} `logLik()` method extracts the empirical log-likelihood
#'   value which can be interpreted as a measure of goodness-of-fit.
setMethod("logLik", "EL", function(object, ...) {
  if (!missing(...)) {
    warning("Extra arguments are not supported.")
  }
  stopifnot(
    "`object` has no `data`. Fit the model with `keep_data == TRUE`." =
      (isFALSE(is.null(getData(object))))
  )
  out <- elt(object, rhs = coef(object))
  new("logLikEL", .Data = logL(out), df = getNumPar(object))
})
