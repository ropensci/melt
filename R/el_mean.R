#' Empirical likelihood test for mean
#'
#' Computes empirical likelihood for mean parameter.
#'
#' @param par Numeric vector of parameters to be tested.
#' @param x Numeric matrix or vector of data. For matrix \code{x}, each row corresponds to an observation.
#' @param maxit Maximum number of iterations for optimization. Defaults to 50.
#' @param abstol Absolute convergence tolerance for optimization. Defaults to 1e-08.
#'
#' @return A list with class \code{c("mean", "melt")}.
#'
#' @references Owen, Art. 1990. “Empirical Likelihood Ratio Confidence Regions.” The Annals of Statistics 18 (1). \doi{10.1214/aos/1176347494}.
#' @examples
#' ## scalar mean
#' par <- 0
#' x <- rnorm(100)
#' el_mean(par, x)
#'
#' ## vector mean
#' x <- matrix(rnorm(100), ncol = 2)
#' par <- c(0, 0)
#' el_mean(par, x)
#' @export
el_mean <- function(par, x, maxit = 1e04, abstol = 1e-08) {
  if (!is.numeric(par)) stop("'par' must be a numeric vector")
  out <- EL_mean(par, x, maxit, abstol)
  out$data.name <- deparse1(substitute(x))
  return(out)
}

#' @export
print.el_test <- function(x, digits = getOption("digits"), prefix = "\t", ...) {
  cat("\n")
  cat(strwrap(x$method, prefix = prefix), sep = "\n")
  cat("\n")
  cat("data: ", x$data.name, "\n", sep = "")
  out <- character()
  if (!is.null(x$statistic))
    out <- c(out, paste("Chisq", names(x$statistic), "=",
                        format(x$statistic, digits = max(1L, digits - 2L))))
  if (!is.null(x$df))
    out <- c(out, paste("df", "=", x$df))
  if (!is.null(x$p.value)) {
    fp <- format.pval(x$p.value, digits = max(1L, digits - 3L))
    out <- c(out, paste("p-value", if (startsWith(fp, "<")) fp else paste("=",
                                                                          fp)))
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  if (!is.null(x$alternative)) {
    cat("alternative hypothesis: ")
    if (!is.null(x$null.value)) {
      if (length(x$null.value) == 1L) {
        alt.char <- switch(x$alternative, two.sided = "not equal to",
                           less = "less than", greater = "greater than")
        cat("true ", names(x$null.value), " is ", alt.char,
            " ", x$null.value, "\n", sep = "")
      }
      else {
        cat(x$alternative, "\nnull values:\n", sep = "")
        print(x$null.value, digits = digits, ...)
      }
    }
    else cat(x$alternative, "\n", sep = "")
  }
  if (!is.null(x$estimate)) {
    cat("maximum EL estimates:\n")
    print(x$estimate, digits = digits, ...)
  }
  cat("\n")
  invisible(x)
}
