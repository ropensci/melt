#' Empirical Likelihood Multiple Hypothesis Testing
#'
#' Tests multiple hypotheses simultaneously for general block designs. Two single step \eqn{k}-FWER (generalized family-wise error rate) controlling procedures are available: asymptotic Monte Carlo (AMC) and nonparametric bootstrap (NB).
#'
#' @param data \code{data.frame} with three variables: blocks (\code{factor}), treatments (\code{factor}), and observations (\code{numeric}).
#' @param hypotheses List of numeric matrices specifying hypotheses to be tested. Each matrix denotes a linear hypothesis in terms of parameters.
#' @param rhs Optional list of numeric vectors specifying the right hand sides of \code{hypotheses}. The length of each vector must match the number of rows of the corresponding matrix. If not specified, they are all set to 0 vectors.
#' @param k Integer value \eqn{k} for \eqn{k}-FWER. Defaults to 1.
#' @param alpha Significance level of the test. Defaults to 0.05.
#' @param method Character value for the procedure to be used; either `AMC' or `NB' is supported. Defaults to `AMC'.
#' @param B Number of random samples for the AMC (number of bootstrap replicates for the NB).
#' @param approx If \code{TRUE}, approximations are used for computing bootstrap statistics. Only applied when the NB is chosen. Defaults to FALSE.
#' @param nthread Number of cores (threads) to be used for bootstrapping. Only applied when the NB is chosen. Defaults to 1.
#' @param maxit Maximum number of iterations for optimization. Defaults to 10000.
#' @param abstol Absolute convergence tolerance for optimization. Defaults to 1e-08.
#'
#' @return A list with class \code{elmulttest}.
#'
#' @examples
#' # General linear hypotheses
#' ELmht(df, c("pairwise"), k = 2, method = "AMC", B = 1e4)
#'
#' # All pairwise comparisons
#' ELmht(df, c("pairwise"), method = "NB", B = 1e3)
#'
#' # Comparisons with control
#' ELmht(df, c("pairwise", "1"), method = "NB", B = 1e3)
#'
#' @importFrom stats reshape
#' @export
ELmht <- function(data, hypotheses, rhs = NULL,
                  k = 1, alpha = 0.05,
                  method = c("AMC", "NB"), B,
                  approx = F, nthread = 1, maxit = 1e04, abstol = 1e-8) {
  ## check method
  method <- match.arg(method)

  ## check data
  # take data.frame with 3 variables
  if (!is.data.frame(data) || length(data) != 3) {
    stop("error")
  }

  # 2 variables for treatments and blocks & 1 variable for observations
  type <- sapply(data, class)
  if (!setequal(type, c("factor", "factor", "numeric"))) {
    stop("error")
  }

  # infer variable types and name them as "block", "trt, and "obs"
  obs.loc <- unname(which(type == "numeric"))
  block.loc <- unname(which(type == "factor"))[1]
  trt.loc <- unname(which(type == "factor"))[2]
  if (length(levels(data[, block.loc])) == length(levels(data[, trt.loc]))) {
    stop("error")
  } else if (length(levels(data[, block.loc])) <
             length(levels(data[, trt.loc]))) {
    block.loc <- unname(which(type == "factor"))[2]
    trt.loc <- unname(which(type == "factor"))[1]
  }
  trt <- levels(data[, trt.loc])

  # incidence matrix
  c <- unclass(table(data[, block.loc], data[, trt.loc]))
  # data matrix
  x <- reshape(data[order(data[, trt.loc]),],
               idvar = names(data)[block.loc],
               timevar = names(data)[trt.loc],
               v.names = names(data)[obs.loc],
               direction = "wide")
  # remove block variable and convert to matrix
  x <- x[order(x[, block.loc]),]
  x[, names(data)[block.loc]] <- NULL
  x <- as.matrix(x)
  # replace NA with 0
  x[is.na(x)] <- 0
  dimnames(x) <- list(levels(data[, block.loc]), trt)


  ## check hypotheses and rhs
  test <- check_hypotheses(hypotheses, rhs, trt)

  ## checked everything and go to function call
  if (is.null(test)) {
    stop("general hypotheses not supported yet")
  } else {
    out <- el_pairwise(x, c, control = test, k, alpha, interval = T,
                       method, B, approx, nthread,
                       maxit, abstol)
    out$trt <- trt
    out$test <- hypotheses
  }
  out
}
