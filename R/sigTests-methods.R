#' @describeIn sigTests Extracts a list with the optimization results of
#'   significance tests.
setMethod("sigTests", "LM", function(object, ...) {
  object@sigTests
})

#' @describeIn sigTests Extracts a matrix with the results of significance
#'   tests.
setMethod("sigTests", "SummaryLM", function(object, ...) {
  object@sigTests
})
