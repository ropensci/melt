#' Empirical Likelihood Multiple Hypothesis Testing
#'
#' @param df A data frame.
#' @param hypotheses hypotheses.
#' @param rhs Right hand side of hypotheses
#' @param k Integer k for k-FWER. Defaults to 1.
#' @param level Level of test.
#' @param method Method to be used; either `AMC' or `NB' is supported. Defaults to `AMC'.
#' @param B Number of replicates for AMC and NB.
#' @param approx Whether to use the approximation for NB. Defaults to FALSE.
#' @param nthread Number of cores(threads) to use for NB. Defaults to 1.
#' @param maxit Maximum number of iterations for optimization. Defaults to 1000.
#' @param abstol Absolute convergence tolerance. Defaults to 1e-08.
#'
#' @examples
#' sum(1:10)
#' sum(1:5, 6:10)
#'
#' @return A list with class \code{elmulttest}.
#'
#' @export
#' @importFrom stats reshape
ELmht <- function(df, hypotheses, rhs = NULL,
                  k = 1, level = 0.05,
                  method = c("AMC", "NB"), B = 1e4,
                  approx = F, nthread = 1, maxit = 1e4, abstol = 1e-8) {
  ## check method
  method <- match.arg(method)

  ## check df
  # take data.frame with 3 variables
  if (!is.data.frame(df) || length(df) != 3) {
    stop("error")
  }

  # 2 variables for treatments and blocks & 1 variable for observations
  type <- sapply(df, class)
  if (!setequal(type, c("factor", "factor", "numeric"))) {
    stop("error")
  }

  # infer variable types and name them as "block", "trt, and "obs"
  obs.loc <- unname(which(type == "numeric"))
  block.loc <- unname(which(type == "factor"))[1]
  trt.loc <- unname(which(type == "factor"))[2]
  if (length(levels(df[, block.loc])) == length(levels(df[, trt.loc]))) {
    stop("error")
  } else if (length(levels(df[, block.loc])) < length(levels(df[, trt.loc]))) {
    block.loc <- unname(which(type == "factor"))[2]
    trt.loc <- unname(which(type == "factor"))[1]
  }
  trt <- levels(df[, trt.loc])

  # incidence matrix
  c <- unclass(table(df[, block.loc], df[, trt.loc]))
  # data matrix
  x <- reshape(df[order(df[, trt.loc]),],
               idvar = names(df)[block.loc],
               timevar = names(df)[trt.loc],
               v.names = names(df)[obs.loc],
               direction = "wide")
  # remove block variable and convert to matrix
  x <- x[order(x[, block.loc]),]
  x[, names(df)[block.loc]] <- NULL
  x <- as.matrix(x)
  # replace NA with 0
  x[is.na(x)] <- 0
  dimnames(x) <- list(levels(df[, block.loc]), trt)


  ## check hypotheses and rhs
  test <- check_hypotheses(hypotheses, rhs, trt)

  ## checked everything and go to function call
  if (is.null(test)) {
    stop("general hypotheses not supported yet")
  } else {
    out <- el_pairwise(x, c, control = test, k, level, interval = T,
                       method, B, approx, nthread,
                       maxit, abstol)
    out$trt <- trt
    out$test <- hypotheses
  }
  out
}
