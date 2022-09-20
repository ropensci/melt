#' @describeIn LM Extracts the symbolic model formula used in [`el_lm()`] or
#'   [`el_glm()`].
#' @exportMethod formula
#' @usage NULL
setMethod("formula", "LM", function(x, ...) {
  formula(x@terms)
})
