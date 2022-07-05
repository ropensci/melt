#' @rdname eld
#' @importFrom methods is
setMethod("eld", "EL", function(object, control = el_control()) {
  stopifnot(
    "`eld` method is not applicable to a `GLM` object." = (!is(object, "GLM")),
    "`object` has no `data`; fit the model with `model` = TRUE." =
      (length(getDataMatrix(object)) > 1L),
    "Invalid `control` specified." = (is(control, "ControlEL"))
  )
  new("ELD", eld = compute_ELD(
    getMethodEL(object), coef(object), getDataMatrix(object), control@maxit_l,
    control@tol_l, control@th, control@nthreads, getWeights(object)
  ))
})

#' @rdname eld
#' @usage NULL
setMethod("eld", "missing", function(object, control = el_control()) {
  stopifnot(
    "Invalid `control` specified." = (is(control, "ControlEL"))
  )
  NULL
})
