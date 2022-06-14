#' @rdname eld
#' @importFrom methods is
setMethod(
  "eld", "EL",
  function(object, control = el_control()) {
    if (is(object, "GLM")) {
      stop("'eld' method is not applicable to a 'GLM' object")
    }
    if (length(object@data) == 0L) {
      stop("'object' has no 'data'; fit the model with 'model' = TRUE")
    }
    if (!is(control, "ControlEL")) {
      stop("invalid 'control' specified")
    }
    new("ELD", eld = eld_(
      getMethodEL(object), coef(object), object@data, control@maxit_l,
      control@tol_l, control@th, control@nthreads, object@weights
    ))
  }
)
