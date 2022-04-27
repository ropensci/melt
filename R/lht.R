#' Linear hypothesis test
#'
#' Tests a linear hypothesis for objects that inherit from class
#'   \linkS4class{EL}.
#'
#' @param object Fitted \linkS4class{EL} object.
#' @param rhs Numeric vector for the right-hand-side of hypothesis, with as
#'   many entries as the rows in \code{lhs}. Defaults to \code{NULL}. See
#'   ‘Details’.
#' @param lhs Numeric matrix, or an object that can be coerced to a numeric
#'   matrix. It specifies the left-hand-side of hypothesis. Each row gives a
#'   linear combination of parameters. The number of columns should be equal to
#'   the number of parameters in \code{object}. Defaults to \code{NULL}.
#'   See ‘Details’.
#' @param control List of control parameters set by \code{\link{el_control}}.
#' @details \code{\link{lht}} performs the constrained minimization of
#'   \eqn{l(\theta)} introduced in \code{\link{el_control}}.
#'   \code{rhs} and \code{lhs} cannot be both \code{NULL}. For non-\code{NULL}
#'   \code{lhs}, it is required that \code{lhs} have full row rank
#'   \eqn{q \leq p} and \eqn{p} be equal to \code{object$npar}, the number of
#'   parameters in the fitted model.
#'
#'   Depending on the specification of \code{rhs} and \code{lhs}, we have the
#'   following three cases:
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
#' @return S4 object of class of class \linkS4class{EL} or \linkS4class{MinEL}.
#' @references Kim, E., MacEachern, S., and Peruggia, M., (2021),
#' "Empirical Likelihood for the Analysis of Experimental Designs,"
#' \href{https://arxiv.org/abs/2112.09206}{arxiv:2112.09206}.
#' @references Qin, Jing, and Jerry Lawless. 1995.
#'   “Estimating Equations, Empirical Likelihood and Constraints on Parameters.”
#'   Canadian Journal of Statistics 23 (2): 145–59. \doi{10.2307/3315441}.
#' @seealso \link{el_control}
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
#' @importFrom methods is
#' @importFrom stats pchisq
#' @export
lht <- function(object, rhs = NULL, lhs = NULL, control = el_control()) {
  if (all(!is(object, c("EL")), !is(object, c("MinEL")))) {
    stop("invalid 'object' supplied")
  }
  if (length(object@dataMatrix) == 0L) {
    stop("'object' has no 'dataMatrix'; fit the model with 'model' = TRUE")
  }
  h <- check_hypothesis(lhs, rhs, object@npar)
  if (!is(control, "ControlEL")) {
    stop("invalid 'control' specified")
  }
  method <- object@optim$method
  maxit <- control@maxit
  maxit_l <- control@maxit_l
  tol <- control@tol
  tol_l <- control@tol_l
  step <- control@step
  th <- control@th
  w <- object@weights
  if (is.null(lhs)) {
    el <- eval_(method, h$r, object@dataMatrix, maxit_l, tol_l, th, w)
    p <- length(h$r)
    return(new("EL",
      optim = el$optim, logp = el$logp, logl = el$logl, loglr = el$loglr,
      statistic = el$statistic, df = p,
      pval = pchisq(el$statistic, df = p, lower.tail = FALSE), npar = p,
      weights = w
    ))
  }
  el <- lht_(
    method, object@coefficients, object@dataMatrix, h$l, h$r,
    maxit, maxit_l, tol, tol_l, step, th, w
  )
  new("MinEL",
    optim = el$optim, logp = el$logp, logl = el$logl, loglr = el$loglr,
    statistic = el$statistic, df = nrow(h$l),
    pval = pchisq(el$statistic, df = nrow(h$l), lower.tail = FALSE),
    npar = ncol(h$l), weights = w
  )
}
