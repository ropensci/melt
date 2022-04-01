#' Empirical likelihood hypothesis testing
#'
#' Tests single hypothesis for general block designs.
#'
#' @param formula A formula object. It must specify variables for response,
#'   treatment, and block as 'response ~ treatment | block'. Note that the use
#'   of vertical bar (|) separating treatment and block.
#' @param data A data frame containing the variables in the formula.
#' @param lhs Numeric matrix specifying linear hypothesis in terms of
#'   parameters.
#' @param rhs Optional numeric vector for the right hand side of \code{lhs}.
#'   If not specified, it is set to 0 vector.
#' @param maxit Maximum number of iterations for optimization.
#'   Defaults to 10000.
#' @param abstol Absolute convergence tolerance for optimization.
#'   Defaults to 1e-08.
#'
#' @return A list with class \code{c("el_test", "melt")}.
#' @references Kim, E., MacEachern, S., and Peruggia, M., (2021),
#' "Empirical Likelihood for the Analysis of Experimental Designs,"
#' \href{https://arxiv.org/abs/2112.09206}{arxiv:2112.09206}.
#'
#' @examples
#' # test of no treatment effect
#' data("clothianidin")
#' el_test(clo ~ trt | blk, clothianidin,
#'         lhs = matrix(c(1, -1, 0, 0,
#'                        0, 1, -1, 0,
#'                        0, 0, 1, -1), byrow = TRUE, nrow = 3))
#'
#' @importFrom stats reshape
#' @export
el_test <- function(formula, data, lhs, rhs = NULL, maxit = 1e04,
                    abstol = 1e-08) {
  ## check formula
  f <- attributes(terms(formula))
  if (any(
    # response required & no arbitrary manipulation on intercept
    f$response == 0, f$intercept == 0,
    length(f$variables) != 3L,
    # no other formula
    typeof(f$variables[[3L]]) != "language" ||
    length(f$variables[[3L]]) != 3L,
    # "|" operator needed
    f$variables[[3L]][[1L]] != "|",
    # no transformation of variables
    typeof(f$variables[[3L]][[2L]]) != "symbol" ||
    typeof(f$variables[[3L]][[3L]]) != "symbol",
    # distinct variables for treatment and block
    f$variables[[3L]][[2L]] == f$variables[[3L]][[3L]])
  ) {
    stop("specify formula as 'response ~ treatment | block'")
  }

  ## pseudo formula for model.frame
  l <- f$variables[[2L]]
  r <- c(f$variables[[3L]][[2L]], f$variables[[3L]][[3L]])
  pf <- formula(paste(l, paste(r, collapse = " + "), sep = " ~ "))

  ## extract model frame
  mf <- match.call()
  mf <- mf[c(1L, match(c("formula", "data"), names(mf), 0L))]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(model.frame)
  mf[[2L]] <- pf
  mf <- eval(mf, parent.frame())
  attributes(mf)$terms <- NULL

  ## type conversion
  # response
  mf[[1L]] <- as.numeric(mf[[1L]])
  # treatment
  mf[[2L]] <- as.factor(mf[[2L]])
  # block
  mf[[3L]] <- as.factor(mf[[3L]])
  if (nlevels(mf[[2L]]) >= nlevels(mf[[3L]])) {
    stop("number of blocks should be larger than number of treatments")
  }

  ## construct general block design
  # incidence matrix
  c <- unclass(table(mf[[3L]], mf[[2L]]))
  # model matrix
  x <- reshape(mf[order(mf[[2L]]), ],
               idvar = names(mf)[3L],
               timevar = names(mf)[2L],
               v.names = names(mf)[1L],
               direction = "wide")
  x <- x[order(x[[names(mf)[3L]]]), ]
  # replace NA with 0
  x[is.na(x)] <- 0
  # remove block variable and convert to matrix
  x[names(mf)[3L]] <- NULL
  x <- as.matrix(x)
  # name rows and columns
  dimnames(x) <- list(levels(mf[[3L]]), levels(mf[[2L]]))
  # general block design
  gbd <-
    list("model_matrix" = x, "incidence_matrix" = c, "trt" = levels(mf[[2L]]))
  class(gbd) <- c("gbd", "melt")

  ## test for lhs and rhs
  if (is.null(rhs)) {
    rhs <- rep(0, nrow(lhs))
  }

  ## test hypothesis
  out <- ELtest(gbd$model_matrix, gbd$incidence_matrix, lhs, rhs,
                threshold = nrow(lhs) * 500, maxit, abstol)
  out$trt <- gbd$trt
  out$model.matrix <- gbd$model_matrix
  out$incidence.matrix <- gbd$incidence_matrix
  class(out) <- c("el_test", oldClass(out))
  if (!out$optim$convergence) {
    warning("convergence failed\n")
  }
  out
}

#' Linear hypothesis test
#'
#' Tests a linear hypothesis for objects inheriting from class \code{"el_test"}.
#'
#' @param object A fitted \code{"el_test"} object.
#' @param rhs A numeric vector for the right-hand-side of hypothesis, with as
#'   many entries as the rows in \code{lhs}. Defaults to a vector of zeroes.
#' @param lhs A numeric matrix, or an object that can be coerced to a numeric
#'   matrix. It specifies the left-hand-side of hypothesis. Each row gives a
#'   linear combination of parameters. The number of columns should be equal to
#'   the number of parameters in \code{object}. See ‘Details’.
#' @param control A list of control parameters. See ‘Details’.
#' @details Consider a linear hypothesis of the form \deqn{L\theta = r,} where
#'   the left-hand-side \eqn{L} is a \eqn{q} by \eqn{p} matrix and the
#'   right-hand-side \eqn{r} is a \eqn{q}-dimensional vector. Let
#'   \eqn{l_n(\theta)} denote the minus twice the empirical log-likelihood ratio
#'   function. Under some regularity conditions, \eqn{l_n(\theta)} is
#'   asymptotically distributed as \eqn{\chi^2_q} under the constraint of
#'   hypothesis, i.e.,
#'   \deqn{\min_{\theta: L\theta = r} l_n(\theta) \to_d \chi^2_q .}
#'   \code{lht} solves the constrained optimization problem using projected
#'   gradient descent method. It is required that \code{lhs} have full row rank
#'   \eqn{q \leq p} and \eqn{p} be equal to \code{object$npar}, the number of
#'   parameters. \code{control} is a list that can supply any of the following
#'   components:
#'   \describe{
#'   \item{maxit}{The maximum number of iterations for the optimization.
#'   Defaults to \code{100}.}
#'   \item{tol}{The relative convergence tolerance, denoted by \eqn{\epsilon}.
#'   With the orthogonal projector matrix \eqn{P} and an initial value
#'   \eqn{\theta^{(0)}}, the iteration stops when
#'   \deqn{\|P \nabla l_n(\theta^{(k)})\| \leq
#'   \epsilon\|P \nabla l_n(\theta^{(0)})\| + \epsilon^2.}
#'   Defaults to \code{1e-06}.}
#'   \item{th}{The threshold for negative empirical log-likelihood ratio value.
#'   The iteration stops if the value exceeds the threshold.
#'   Defaults to \code{NULL} and sets the threshold to \eqn{200p}.}
#'   }
#' @return A list with class \code{c("el_lht", "el_test")} with the following
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
#' @seealso \link{el_eval}
#' @examples
#' n <- 100
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' y <- 1 + x1 + x2 + rnorm(n)
#' df <- data.frame(y, x1, x2)
#' fit <- el_lm(y ~ x1 + x2, df)
#' lhs <- matrix(c(0, 1, -1), nrow = 1)
#' lht(fit, lhs = lhs)
#'
#' # test of no treatment effect
#' data("clothianidin")
#' lhs2 <- matrix(c(1, -1, 0, 0,
#'                  0, 1, -1, 0,
#'                  0, 0, 1, -1), byrow = TRUE, nrow = 3)
#' fit2 <- el_lm(clo ~ -1 + trt, clothianidin)
#' lht(fit2, lhs = lhs2)
#' @export
lht <- function(object, rhs, lhs, control = list()) {
  if (!inherits(object, "el_test"))
    stop("invalid 'object' supplied")
  if (is.null(object$data.matrix))
    stop("'object' has no 'data.matrix'; fit the model with 'model' = TRUE")

  lhs <- as.matrix(lhs)
  q <- check_hypothesis(lhs, object$npar)

  if (missing(rhs)) {
    rhs <- rep(0, nrow(lhs))
  } else {
    rhs <- as.vector(rhs)
    if (!is.numeric(rhs) || !all(is.finite(rhs)))
      stop("'rhs' must be a finite numeric vector")
    if (length(rhs) != q)
      stop("'lhs' and 'rhs' have incompatible dimensions")
  }

  optcfg <- check_control(control)
  if (is.null(object$weights)) {
    out <- lht_(object$optim$method, object$coefficients, object$data.matrix,
                lhs, rhs, optcfg$maxit, optcfg$tol, optcfg$th)
  } else {
    out <- lht_w_(object$optim$method, object$coefficients, object$data.matrix,
                  object$weights, lhs, rhs, optcfg$maxit, optcfg$tol, optcfg$th)
  }
  class(out) <- c("el_lht", "el_test")
  out
}

#' Confidence intervals for model parameters
#'
#' Computes confidence intervals for one or more parameters in a fitted model.
#' Package \strong{melt} adds a method for objects inheriting from class
#' \code{"el_test"}.
#'
#' @param object A fitted \code{"el_test"} object.
#' @param parm A specification of which parameters are to be given confidence
#'   intervals, either a vector of numbers or a vector of names. If missing, all
#'   parameters are considered.
#' @param level A confidence level required.
#' @param control A list of control parameters. See ‘Details’ in
#'   \code{\link{lht}}.
#' @param ... Additional argument(s) for methods.
#' @importFrom stats complete.cases qchisq
#' @return A matrix with columns giving lower and upper confidence limits for
#'  each parameter. In contrast to other methods that rely on studentization,
#'  the lower and upper limits obatained from empirical likelihood do not
#'  correspond to the (1 - level) / 2 and 1 - (1 - level) / 2 in \%,
#'  respectively.
#' @references Kim, E., MacEachern, S., and Peruggia, M., (2021),
#' "Empirical Likelihood for the Analysis of Experimental Designs,"
#' \href{https://arxiv.org/abs/2112.09206}{arxiv:2112.09206}.
#' @seealso \link{lht}
#' @examples
#' fit <- el_lm(formula = mpg ~ wt, data = mtcars)
#' confint(fit)
#' @export
confint.el_test <- function(object, parm, level = 0.95, control = list(), ...) {
  # check level and control arguments
  if (!missing(level) &&
      (length(level) != 1L || !is.finite(level) || level < 0 || level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  optcfg <- check_control(control)
  # set cutoff and coefficients
  cutoff <- qchisq(level, 1L)
  cf <- coef(object)
  # no confidence interval for empty model
  if (length(cf) == 0L) {
    ci <- matrix(, nrow = 0L, ncol = 2L)
    colnames(ci) <- c("lower", "upper")
    return(ci)
  }
  # index for tracking the parameters
  idx <- seq(length(cf))
  # rownames of the confidence interval matrix
  pnames <- if (is.null(names(cf))) idx else names(cf)
  # if parm is supplied, modify idx and pnames accordingly
  if (!missing(parm)) {
    if (is.numeric(parm) && all(is.finite(parm))) {
      pnames <- pnames[parm]
      # idx <- match(pnames, names(cf))
      idx <- if (is.null(names(cf))) match(pnames, idx) else
        match(pnames, names(cf))
    } else if (is.character(parm)) {
      idx <- match(parm, pnames)
      pnames <- parm
    } else {
      stop("invalid 'parm' specified")
    }
  }
  # number of rows of the confidence interval matrix
  p <- length(idx)
  # compute the confidence interval matrix
  if (level == 0) {
    ci <- matrix(rep(cf[idx], 2L), ncol = 2L)
  } else if (level == 1) {
    ci <- matrix(NA, nrow = p, ncol = 2L)
    ci[which(!is.na(idx)), ] <- c(-Inf, Inf)
  } else if (all(is.na(idx))) {
    ci <- matrix(NA, nrow = p, ncol = 2L)
  } else if (any(is.na(idx))) {
    idx_na <- which(is.na(idx))
    ci <- matrix(NA, nrow = p, ncol = 2L)
    # ci[-idx_na, ] <- confint_(object$optim$method, cf, object$data.matrix,
    #                           cutoff, idx[-idx_na], optcfg$maxit, optcfg$tol,
    #                           optcfg$th)
    if (is.null(object$weights)) {
      ci[-idx_na, ] <- confint_(object$optim$method, cf, object$data.matrix,
                                cutoff, idx[-idx_na], optcfg$maxit, optcfg$tol,
                                optcfg$th)
    } else {
      ci[-idx_na, ] <- confint_w_(object$optim$method, cf, object$data.matrix,
                                  object$weights, cutoff, idx[-idx_na],
                                  optcfg$maxit, optcfg$tol, optcfg$th)
    }
  } else {
    # ci <- confint_(object$optim$method, cf, object$data.matrix, cutoff, idx,
    #                optcfg$maxit, optcfg$tol, optcfg$th)
    if (is.null(object$weights)) {
      ci <- confint_(object$optim$method, cf, object$data.matrix, cutoff, idx,
                     optcfg$maxit, optcfg$tol, optcfg$th)
    } else {
      ci <- confint_w_(object$optim$method, cf, object$data.matrix,
                       object$weights,  cutoff, idx, optcfg$maxit, optcfg$tol,
                       optcfg$th)
    }
  }
  dimnames(ci) <- list(pnames, c("lower", "upper"))
  ci
}

#' Empirical log-likelihood
#'
#' Computes the empirical log-likelihood value of the model represented by
#'   \code{object} evaluated at the estimated coefficients. Package
#'   \strong{melt} adds a method for objects inheriting from class
#'   \code{"el_test"}.
#'
#' @param object A fitted \code{"el_test"} object.
#' @param ... Some methods for this generic function require extra arguments.
#'   None are used in this method.
#' @return An object of class \code{"logLik"} with an attribute \code{df} that
#'   gives the number of (estimated) parameters in the model.
#' @examples
#' fit <- el_lm(formula = mpg ~ wt, data = mtcars)
#' logLik(fit)
#' @export
logLik.el_test <- function(object, ...) {
  if (!missing(...))
    warning("extra arguments are not supported")
  p <- object$npar
  if (inherits(object, "el_lht")) {
    val <- object$loglik
    # df is the number of estimated parameters
    attr(val, "df") <- p - object$df
  } else{
    rhs <- object$coefficients
    lhs <- diag(nrow = p, ncol = p)
    out <- lht(object, rhs = rhs, lhs = lhs)
    val <- out$loglik
    attr(val, "df") <- p
  }
  class(val) <- "logLik"
  val
}

#' @noRd
#' @export
print.el_test <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nEmpirical Likelihood Test:", x$optim$method, "\n\n")
  out <- character()
  if (!is.null(x$statistic)) {
    out <- c(out, paste("Chisq:", format(x$statistic, digits = digits)),
             paste("df:", x$df),
             paste("p-value:", format.pval(x$p.value, digits = digits)))
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  if (!is.null(x$coefficients)) {
    cat("maximum EL estimates:\n")
    print(x$coefficients, digits = digits, ...)
  }
  cat("\n")
  invisible(x)
}

#' @noRd
#' @export
print.el_lht <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nEmpirical Likelihood Linear Hypothesis Test:", x$optim$method, "\n\n")
  out <- character()
  if (!is.null(x$statistic)) {
    out <- c(out, paste("Chisq:", format(x$statistic, digits = digits)),
             paste("df:", x$df),
             paste("p-value:", format.pval(x$p.value, digits = digits)))
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  if (!is.null(x$coefficients)) {
    cat("maximum EL estimates:\n")
    print(x$coefficients, digits = digits, ...)
  }
  cat("\n")
  invisible(x)
}
