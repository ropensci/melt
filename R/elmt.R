#' Multiple tests with empirical likelihood
#'
#' Tests a linear hypothesis for objects that inherit from class
#'   \linkS4class{EL}.
#'
#' @param object A fitted \linkS4class{EL} object.
#' @param rhs A numeric vector for the right-hand-side of hypothesis, with as
#'   many entries as the rows in \code{lhs}. Defaults to \code{NULL}. See
#'   ‘Details’.
#' @param lhs A numeric matrix, or an object that can be coerced to a numeric
#'   matrix. It specifies the left-hand-side of hypothesis. Each row gives a
#'   linear combination of parameters. The number of columns should be equal to
#'   the number of parameters in \code{object}. Defaults to \code{NULL}.
#'   See ‘Details’.
#' @param alpha level.
#' @param control A list of control parameters set by \code{\link{el_control}}.
#' @details \code{\link{elt}} performs the constrained minimization of
#'   \eqn{l(\theta)} described in \linkS4class{CEL}. \code{rhs} and \code{lhs}
#'   cannot be both \code{NULL}. For non-\code{NULL} \code{lhs}, it is required
#'   that \code{lhs} have full row rank \eqn{q \leq p} and \eqn{p} be equal to
#'   \code{object$npar}, the number of parameters in the fitted model.
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
#' @return If \code{lhs} is \code{NULL}, an object of class \linkS4class{EL}
#' is returned. Otherwise, an object of class \linkS4class{CEL} is returned.
#' @references Adimari, Gianfranco, and Annamaria Guolo. 2010.
#'   “A Note on the Asymptotic Behaviour of Empirical Likelihood Statistics.”
#'   Statistical Methods & Applications 19 (4): 463–76.
#'   \doi{10.1007/s10260-010-0137-9}.
#' @references Qin, Jing, and Jerry Lawless. 1995.
#'   “Estimating Equations, Empirical Likelihood and Constraints on Parameters.”
#'   Canadian Journal of Statistics 23 (2): 145–59. \doi{10.2307/3315441}.
#' @seealso \link{el_control}
#' @examples
#' n <- 100L
#' @importFrom methods is
#' @importFrom stats pchisq
#' @export
elmt <- function(object,
                 rhs = NULL,
                 lhs = NULL,
                 alpha = 0.05,
                 control = el_control()) {
  stopifnot(
    "invalid 'object' supplied" = (is(object, "EL")),
    "'object' has no 'data'; fit the model with 'model' = TRUE" =
      (length(getDataMatrix(object)) > 1L),
    "invalid 'control' specified" = (is(control, "ControlEL"))
  )
  p <- object@npar
  h <- validate_hypotheses(rhs, lhs, p)
  # return(h)

  validate_alpha(alpha)
  method <- getMethodEL(object)
  maxit <- control@maxit
  maxit_l <- control@maxit_l
  tol <- control@tol
  tol_l <- control@tol_l
  step <- control@step
  th <- control@th


  elmt_statistic_(h$q, h$m, method, coef(object), getDataMatrix(object), h$r, h$l,
                  maxit, maxit_l, tol, tol_l, step, th, getWeights(object))
}
