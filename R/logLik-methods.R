#' @rdname logLik
setMethod("logLik", "EL", function(object, ...) {
  if (!missing(...)) {
    warning("Extra arguments are not supported.")
  }
  stopifnot(
    "`object` has no `data`. Fit the model with `keep_data == TRUE`." =
      (!is.null(getDataMatrix(object)))
  )
  p <- object@npar
  out <- elt(object, rhs = coef(object))
  val <- out@logl
  new("logLikEL", logLik = val, df = p)
})
