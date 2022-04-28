#' @rdname logLik
setMethod(
  "logLik", "EL",
  function(object, ...) {
    if (!missing(...)) {
      warning("extra arguments are not supported")
    }
    if (length(object@dataMatrix) == 0L) {
      stop("method is not applicable to an empty model")
    }
    p <- object@npar
    rhs <- object@coefficients
    out <- lht(object, rhs = rhs)
    val <- out@logl
    new("logLikEL", logLik = val, df = p)
  }
)
