#' @rdname logLik
#' @srrstats {RE4.11} `logLik()` method extracts the empirical log-likelihood
#'   value which can be interpreted as a measure of goodness-of-fit.
setMethod("logLik", "EL", function(object, ...) {
  if (!missing(...)) {
    warning("Extra arguments are not supported.")
  }
  stopifnot(
    "`object` has no `data`. Fit the model with `keep_data == TRUE`." =
      (!is.null(getDataMatrix(object)))
  )
  out <- elt(object, rhs = coef(object))
  val <- out@logl
  new("logLikEL", logLik = val, df = getNumPar(object))
})
