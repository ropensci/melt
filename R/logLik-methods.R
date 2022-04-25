# setMethod(
#   "logLik", "EL",
#   function(object, ...) {
#     if (!missing(...)) {
#       warning("extra arguments are not supported")
#     }
#     p <- object@npar
#     rhs <- object@coefficients
#     out <- lht(object, rhs = rhs)
#     val <- out$loglik
#     attr(val, "df") <- p
#     val
#   }
# )
