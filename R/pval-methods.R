#' @describeIn pval Extracts the \eqn{p}-value.
setMethod("pval", "EL", function(object, ...) {
  object@pval
})

#' @describeIn pval Extracts the \eqn{p}-value.
setMethod("pval", "ELT", function(object, ...) {
  object@pval
})

#' @describeIn pval Extracts the multiplicity adjusted \eqn{p}-values.
setMethod("pval", "ELMT", function(object, ...) {
  object@pval
})
