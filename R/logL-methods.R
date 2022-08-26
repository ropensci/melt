#' @rdname logL
#' @srrstats {RE4.11} `logL()` method extracts the empirical log-likelihood
#'   value which can be interpreted as a measure of goodness-of-fit.
setMethod("logL", "EL", function(object, ...) {
  object@logl
})

#' @rdname logL
setMethod("logL", "ELT", function(object, ...) {
  object@logl
})
