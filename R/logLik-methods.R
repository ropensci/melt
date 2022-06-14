#' @rdname logLik
setMethod(
  "logLik", "EL",
  function(object, ...) {
    if (!missing(...)) {
      warning("extra arguments are not supported")
    }
    if (length(object@data) == 0L) {
      stop("method is not applicable to an empty model")
    }
    p <- object@npar
    out <- lht(object, rhs = coef(object))
    val <- out@logl
    new("logLikEL", logLik = val, df = p)
  }
)
