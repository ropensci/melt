#' @describeIn getSigTests Extracts a list with the optimization results of
#'   significance tests.
setMethod("getSigTests", "LM", function(object) {
  object@sigTests
})


#' @describeIn getSigTests Extracts a matrix with the results of significance
#'   tests.
setMethod("getSigTests", "SummaryLM", function(object) {
  object@sigTests
})
