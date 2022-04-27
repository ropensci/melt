#' @rdname eld
#' @importFrom methods is
setMethod(
  "eld", "EL",
  function(object, control = el_control()) {
    if (is(object, "GLM")) {
      stop("'eld' method is not applicable to a 'GLM' object")
    }
    if (length(object@dataMatrix) == 0L) {
      stop("'object' has no 'dataMatrix'; fit the model with 'model' = TRUE")
    }
    if (!is(control, "ControlEL")) {
      stop("invalid 'control' specified")
    }
    new("ELD", eld = eld_(
      object@optim$method, object@coefficients, object@dataMatrix,
      control@maxit_l, control@tol_l, control@th, control@nthreads,
      object@weights
    ))
  }
)
