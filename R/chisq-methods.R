#' @describeIn chisq Extracts the chi-square statistic.
setMethod("chisq", "EL", function(object, ...) {
  object@statistic
})

#' @describeIn chisq Extracts the vector of chi-square statistics.
setMethod("chisq", "ELMT", function(object, ...) {
  object@statistic
})

#' @describeIn chisq Extracts the chi-square statistic.
setMethod("chisq", "ELT", function(object, ...) {
  object@statistic
})

#' @describeIn chisq Extracts the chi-square statistic for the overall test of
#'   the model.
setMethod("chisq", "SummaryLM", function(object, ...) {
  object@statistic
})
