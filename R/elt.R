#' Empirical likelihood test
#'
#' Tests a linear hypothesis for objects that inherit from class
#'   \linkS4class{EL}.
#'
#' @param object An object of class \linkS4class{EL}.
#' @param rhs A numeric vector or a column matrix for the right-hand-side of
#'   hypothesis, with as many entries as the rows in \code{lhs}. Defaults to
#'   \code{NULL}. See ‘Details’.
#' @param lhs A numeric matrix or a vector (treated as a row matrix) for the
#'   left-hand-side of hypothesis. Each row gives a linear combination of the
#'   parameters in \code{object}. The number of columns should be equal to the
#'   number of parameters. Defaults to \code{NULL}. See ‘Details’.
#' @param alpha A single numeric for the significance level. Defaults to \code{0.05}.
#' @param calibrate A single character for the calibration method. Defaults to
#' \code{chisq}. See ‘Details’.
#' @param control An object of class \linkS4class{ControlEL} constructed by
#'   \code{\link{el_control}}.
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
#'   \deqn{\inf_{\theta: L\theta = r} l(\theta).}
#'   \item If \code{rhs} is \code{NULL}, \eqn{r} is set to the zero vector as
#'   \deqn{\inf_{\theta: L\theta = 0} l(\theta).}
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
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' y <- 1 + x1 + x2 + rnorm(n)
#' df <- data.frame(y, x1, x2)
#' fit <- el_lm(y ~ x1 + x2, df)
#' elt(fit, lhs = c(0, 1, -1))
#' elt(fit, lhs = c(0, 1, 1), rhs = 2)
#'
#' # test of no treatment effect
#' data("clothianidin")
#' lhs <- matrix(c(
#'   1, -1, 0, 0,
#'   0, 1, -1, 0,
#'   0, 0, 1, -1
#' ), byrow = TRUE, nrow = 3)
#' fit2 <- el_lm(clo ~ -1 + trt, clothianidin)
#' elt(fit2, lhs = lhs)
#' @importFrom methods is
#' @importFrom stats pchisq
#' @export
elt <- function(object, rhs = NULL, lhs = NULL, alpha = 0.05,
                calibrate = c("chisq", "boot", "f"), control = el_control()) {
  stopifnot(
    "invalid 'object' supplied" = (is(object, "EL")),
    "'object' has no 'data'; fit the model with 'model' = TRUE" =
      (length(getDataMatrix(object)) > 1L),
    "invalid 'control' specified" = (is(control, "ControlEL"))
  )
  h <- check_hypothesis_(lhs, rhs, object@npar)
  alpha <- check_alpha_(alpha)
  calibrate <- match.arg(calibrate)
  method <- getMethodEL(object)
  maxit <- control@maxit
  maxit_l <- control@maxit_l
  tol <- control@tol
  tol_l <- control@tol_l
  step <- control@step
  th <- control@th
  w <- getWeights(object)
  if (is.null(lhs)) {
    if (calibrate == "f" && method != "mean") {
      stop("F calibration is applicable only to the mean")
    }
    par <- h$r
    el <- eval_(method, par, getDataMatrix(object), maxit_l, tol_l, th, w)
    p <- length(par)
    cal <- calibrate_(alpha, el$statistic, calibrate, p, par, object, control)
    return(new("ELT",
      optim = el$optim, alpha = alpha, logl = el$logl, statistic = el$statistic,
      cv = cal["cv"], pval = cal["pval"], calibrate = calibrate
    ))
  }
  # proceed with chi-square calibration for non-NULL 'lhs'
  stopifnot(
    "bootstrap calibration is applicable only when 'lhs' is NULL" =
      (calibrate != "boot"),
    "F calibration is applicable only when 'lhs' is NULL" = (calibrate != "f")
  )
  el <- elt_(
    method, coef(object), getDataMatrix(object), h$l, h$r,
    maxit, maxit_l, tol, tol_l, step, th, w
  )
  new("ELT",
    optim = el$optim, alpha = alpha, logl = el$logl, statistic = el$statistic,
    cv = qchisq(p = 1 - alpha, df = nrow(h$l)),
    pval = pchisq(el$statistic, df = nrow(h$l), lower.tail = FALSE),
    calibrate = calibrate
  )
}
