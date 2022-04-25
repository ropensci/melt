# setMethod(
#   "print", "EL",
#   function(y, digits = max(3L, getOption("digits") - 3L), ...) {
#     cat("\nEmpirical Likelihood Test:", x@optim$method, "\n\n")
#     out <- character()
#     if (!is.null(x@statistic)) {
#       out <- c(
#         out, paste("Chisq:", format(x@statistic, digits = digits)),
#         paste("df:", x@df),
#         paste("p-value:", format.pval(x@pval, digits = digits))
#       )
#     }
#     cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
#     if (!is.null(x@coefficients)) {
#       cat("maximum EL estimates:\n")
#       print(x@coefficients, digits = digits, ...)
#     }
#     cat("\n")
#     invisible(x)
#   }
# )
