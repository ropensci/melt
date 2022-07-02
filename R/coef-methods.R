#' @rdname coef
#' @srrstats {RE4.2} `coef()` method is exported.
setMethod("coef", "EL", function(object, ...) {
  object@coefficients
})
