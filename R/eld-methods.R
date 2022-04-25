#' Empirical likelihood displacement
#'
#' Computes empirical likelihood displacement for model diagnostics and outlier
#'   detection.
#'
#' @param object Fitted \linkS4class{EL} object.
#' @param control List of control parameters set by \code{\link{control_el}}.
#' @details Let \eqn{L(\theta)} be the empirical log-likelihood function based
#'   on the full sample with \eqn{n} observations. The maximum empirical
#'   likelihood estimate is denoted by \eqn{\hat{\theta}}. Consider a reduced
#'   sample with the \eqn{i}th observation deleted and the corresponding
#'   estimate \eqn{\hat{\theta}_{(i)}}. The empirical likelihood displacement is
#'   defined by
#'   \deqn{\textnormal{ELD}_i = 2\{L(\hat{\theta}) - L(\hat{\theta}_{(i)})\}.}
#'   If \eqn{\textnormal{ELD}_i } is large, then the \eqn{i}th observation is an
#'   influential point and can be inspected as a possible outlier. \code{eld}
#'   computes \eqn{\textnormal{ELD}_i } for \eqn{i = 1, \dots, n }.
#' @return S4 object of class \linkS4class{ELD}.
#' @references Lazar, Nicole A. 2005. “Assessing the Effect of Individual Data
#'   Points on Inference From Empirical Likelihood.” Journal of Computational
#'   and Graphical Statistics 14 (3): 626–42.
#'   \doi{10.1198/106186005X59568}.
#' @references Zhu, H., J. G. Ibrahim, N. Tang, and H. Zhang. 2008. “Diagnostic
#'   Measures for Empirical Likelihood of General Estimating Equations.”
#'   Biometrika 95 (2): 489–507.
#'   \doi{10.1093/biomet/asm094}.
#' @seealso \link{el_eval}, \link{control_el}, \link{plot}
#' @examples
#' x <- rnorm(10L)
#' y <- 10
#' fit <- el_mean2(0, c(x, y))
#' eld(fit)
#' @aliases eld
setMethod(
  "eld", "EL",
  function(object, control = control_el()) {
    if (!inherits(control, "control_el") || !is.list(control)) {
      stop("invalid 'control' supplied")
    }
    if (inherits(object, "el_glm")) {
      stop("method is not applicable to 'el_glm' object")
    }
    if (length(object@dataMatrix) == 0L) {
      stop("'object' has no 'dataMatrix'; fit the model with 'model' = TRUE")
    }
    new("ELD", eld = eld_(
      object@optim$method, object@coefficients, object@dataMatrix,
      control$maxit_l, control$tol_l, control$th, control$nthreads,
      object@weights
    ))
  }
)
