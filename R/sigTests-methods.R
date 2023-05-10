#' @describeIn sigTests Extracts a list with the optimization results of
#'   significance tests.
setMethod("sigTests", "LM", function(object, ...) {
  object@sigTests
})
