#' Empirical Likelihood Hypothesis Testing
#'
#' Tests single hypothesis for general block designs.
#'
#' @param data \code{data.frame} with three variables: blocks (\code{factor}), treatments (\code{factor}), and observations (\code{numeric}).
#' @param lhs Numeric matrix specifying linear hypothesis in terms of parameters.
#' @param rhs Optional numeric vector for the right hand side of \code{lhs}. If not specified, it is set to 0 vector.
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
#' sum(1)
#'
#' @importFrom stats reshape
#' @export
el_test <- function(data, lhs, rhs = NULL,
                    alpha = 0.05,
                    method = c("AMC", "NB"), B,
                    approx = F, nthread = 1, maxit = 1e04, abstol = 1e-8) {
  ## check method
  method <- match.arg(method)
}
