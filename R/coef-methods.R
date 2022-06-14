#' @rdname coef
setMethod("coef", "EL", function(object, ...) {object@coefficients})
