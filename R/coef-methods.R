#' @rdname coef-method
setMethod(
  "coef", "EL",
  function(object, ...) {
    object@coefficients
  }
)
