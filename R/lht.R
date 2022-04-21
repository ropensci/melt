#' Linear hypothesis test
#'
#' Tests a linear hypothesis for objects inheriting from class \code{"el"}.
#'
#' @param object A fitted \code{"el"} object.
#' @param rhs A numeric vector for the right-hand-side of hypothesis, with as
#'   many entries as the rows in \code{lhs}. Defaults to \code{NULL}. See
#'   ‘Details’.
#' @param lhs A numeric matrix, or an object that can be coerced to a numeric
#'   matrix. It specifies the left-hand-side of hypothesis. Each row gives a
#'   linear combination of parameters. The number of columns should be equal to
#'   the number of parameters in \code{object}. Defaults to \code{NULL}.
#'   See ‘Details’.
#' @param control A list of control parameters set by
#'   \code{\link{control_el}}.
#' @details \code{\link{lht}} performs the constrained minimization of
#'   \eqn{l(\theta)} introduced in \code{\link{control_el}}.
#'   \code{rhs} and \code{lhs} cannot be both \code{NULL}. For non-\code{NULL}
#'   \code{lhs}, it is required that \code{lhs} have full row rank
#'   \eqn{q \leq p} and \eqn{p} be equal to \code{object$npar}, the number of
#'   parameters in the fitted model.
#'
#'   Depending on the specification of \code{rhs} and \code{lhs}, we have the following three cases:
#'   \enumerate{
#'   \item If both \code{rhs} and \code{lhs} are non-\code{NULL}, the
#'   constrained minimization is performed with the right-hand-side \eqn{r} and
#'   the left-hand-side \eqn{L} as
#'   \deqn{\min_{\theta: L\theta = r} l(\theta).}
#'   \item If \code{rhs} is \code{NULL}, \eqn{r} is set to the zero vector as
#'   \deqn{\min_{\theta: L\theta = 0} l(\theta).}
#'   \item If \code{lhs} is \code{NULL}, \eqn{L} is set to the identity matrix
#'   and the problem reduces to evaluating at \eqn{r} as
#'   \deqn{l(r).}
#'   }
#' @return A list of class \code{"elt"} with the following components:
#'   \item{optim}{A list with the following optimization results:
#'     \itemize{
#'       \item{\code{method } }{A character for method dispatch in internal
#'       functions.}
#'       \item{\code{par } }{The solution of the constrained minimization
#'       problem.}
#'       \item{\code{lambda } }{The Lagrange multiplier of dual problem.}
#'       \item{\code{logLR } }{The (weighted) empirical log-likelihood ratio
#'       value.}
#'       \item{\code{iterations } }{The number of iterations performed.}
#'       \item{\code{convergence } }{A logical vector. \code{TRUE} indicates
#'       convergence of the algorithm.}
#'     }
#'   }
#'   \item{log.prob}{The log probabilities.}
#'   \item{loglik}{The log likelihood value.}
#'   \item{statistic}{The chi-square statistic.}
#'   \item{df}{The degrees of freedom of the statistic.}
#'   \item{p.value}{The \eqn{p}-value of the statistic.}
#'   \item{npar}{The number of parameters.}
#' @references Kim, E., MacEachern, S., and Peruggia, M., (2021),
#' "Empirical Likelihood for the Analysis of Experimental Designs,"
#' \href{https://arxiv.org/abs/2112.09206}{arxiv:2112.09206}.
#' @references Qin, Jing, and Jerry Lawless. 1995.
#'   “Estimating Equations, Empirical Likelihood and Constraints on Parameters.”
#'   Canadian Journal of Statistics 23 (2): 145–59. \doi{10.2307/3315441}.
#' @seealso \link{control_el}
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
#' lhs2 <- matrix(c(
#'   1, -1, 0, 0,
#'   0, 1, -1, 0,
#'   0, 0, 1, -1
#' ), byrow = TRUE, nrow = 3L)
#' fit2 <- el_lm(clo ~ -1 + trt, clothianidin)
#' lht(fit2, lhs = lhs2)
#' @export
lht <- function(object, rhs = NULL, lhs = NULL, control = control_el()) {
  if (!inherits(object, "el")) {
    stop("invalid 'object' supplied")
  }
  if (is.null(object$data.matrix)) {
    stop("'object' has no 'data.matrix'; fit the model with 'model' = TRUE")
  }
  h <- check_hypothesis(lhs, rhs, object$npar)
  if (!inherits(control, "control_el") || !is.list(control)) {
    stop("invalid 'control' supplied")
  }
  method <- object$optim$method
  maxit <- control$maxit
  maxit_l <- control$maxit_l
  tol <- control$tol
  tol_l <- control$tol_l
  step <- control$step
  th <- control$th
  w <- object$weights
  if (is.null(w)) {
    w <- numeric(length = 0L)
  }
  if (is.null(lhs)) {
    out <- eval_(method, h$r, object$data.matrix, maxit_l, tol_l, th, w)
    out$df <- length(h$r)
    out$p.value <- pchisq(out$statistic, df = out$df, lower.tail = FALSE)
    out$npar <- length(h$r)
    # out$coefficients <- object$coefficients
    class(out) <- "elt"
  } else {
    out <- lht_(
      method, object$coefficients, object$data.matrix, h$l, h$r,
      maxit, maxit_l, tol, tol_l, step, th, w
    )
    out$df <- nrow(h$l)
    out$p.value <- pchisq(out$statistic, df = out$df, lower.tail = FALSE)
    out$npar <- ncol(h$l)
    class(out) <- "elt"
  }
  out
}
