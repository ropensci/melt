#' Linear hypothesis test
#'
#' Tests a linear hypothesis for objects inheriting from class \code{"el"}.
#'
#' @param object A fitted \code{"el"} object.
#' @param rhs A numeric vector for the right-hand-side of hypothesis, with as
#'   many entries as the rows in \code{lhs}. Defaults to \code{NULL}. See
#'   ‘Details’ in \code{\link{melt_control}}.
#' @param lhs A numeric matrix, or an object that can be coerced to a numeric
#'   matrix. It specifies the left-hand-side of hypothesis. Each row gives a
#'   linear combination of parameters. The number of columns should be equal to
#'   the number of parameters in \code{object}. Defaults to \code{NULL}.
#'   See ‘Details’ in \code{\link{melt_control}}.
#' @param control A list of control parameters set by
#'   \code{\link{melt_control}}.
#' @details \code{lht} performs the constrained minimization of
#'   \eqn{l(\theta)} using projected gradient descent method. It is required
#'   that \code{lhs} have full row rank \eqn{q \leq p} and \eqn{p} be equal to
#'   \code{object$npar}, the number of parameters.
#' @return A list of class \code{c("elt", "el")} with the following
#'   components:
#'   \describe{
#'   \item{optim}{A list with the following optimization results:
#'     \describe{
#'       \item{method}{The type of estimating function.}
#'       \item{lambda}{The Lagrange multiplier of dual problem.}
#'       \item{logLR}{The (weighted) empirical log-likelihood ratio value.}
#'       \item{iterations}{The number of iterations performed.}
#'       \item{convergence}{A logical vector. \code{TRUE} indicates
#'       convergence of the algorithm.}
#'     }
#'   }
#'   \item{npar}{The number of parameters.}
#'   \item{log.prob}{The log probabilities.}
#'   \item{loglik}{The log likelihood value evaluated at the estimated
#'   coefficients}
#'   \item{coefficients}{The solution of the optimization.}
#'   \item{statistic}{The chi-square statistic.}
#'   \item{df}{The degrees of freedom of the statistic.}
#'   \item{p.value}{The \eqn{p}-value of the statistic.}
#' }
#' @references Kim, E., MacEachern, S., and Peruggia, M., (2021),
#' "Empirical Likelihood for the Analysis of Experimental Designs,"
#' \href{https://arxiv.org/abs/2112.09206}{arxiv:2112.09206}.
#' @references Qin, Jing, and Jerry Lawless. 1995.
#'   “Estimating Equations, Empirical Likelihood and Constraints on Parameters.”
#'   Canadian Journal of Statistics 23 (2): 145–59.
#'   \doi{10.2307/3315441}.
#' @seealso \link{melt_control}
#' @examples
#' n <- 100L
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' y <- 1 + x1 + x2 + rnorm(n)
#' df <- data.frame(y, x1, x2)
#' fit <- el_lm(y ~ x1 + x2, df)
#' lhs <- matrix(c(0, 1, -1), nrow = 1L)
#' lht(fit, lhs = lhs)
#'
#' # test of no treatment effect
#' data("clothianidin")
#' lhs2 <- matrix(c(1, -1, 0, 0,
#'                  0, 1, -1, 0,
#'                  0, 0, 1, -1), byrow = TRUE, nrow = 3L)
#' fit2 <- el_lm(clo ~ -1 + trt, clothianidin)
#' lht(fit2, lhs = lhs2)
#' @export
lht <- function(object, rhs = NULL, lhs = NULL, control = melt_control()) {
  if (!inherits(object, "el"))
    stop("invalid 'object' supplied")
  if (is.null(object$data.matrix))
    stop("'object' has no 'data.matrix'; fit the model with 'model' = TRUE")
  h <- check_hypothesis(lhs, rhs, object$npar)
  if (!inherits(control, "melt_control") || !is.list(control))
    stop("invalid 'control' supplied")
  method <- object$optim$method
  maxit <- control$maxit
  tol <- control$tol
  th <- control$th
  w <- object$weights
  if (is.null(lhs)) {
    out <- eval_(method, h$r, object$data.matrix, maxit, tol, th, w)
    class(out) <- "el"
  } else {
    out <- lht_(method, object$coefficients, object$data.matrix, h$l, h$r,
                maxit, tol, th, w)
    class(out) <- c("elt", "el")
  }
  out
}
