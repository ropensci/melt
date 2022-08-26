#' @rdname logLik
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
