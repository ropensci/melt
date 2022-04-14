#' Fit an one-way analysis of variance model with empirical likelihood
#'
#' \emph{This function is deprecated in favor of \link{el_lm} and will be
#'   removed in a future release.}
#'
#' @param formula A formula object. It must specify variables for response and
#'   treatment as 'response ~ treatment'.
#' @param data A data frame containing the variables in the formula.
#' @param maxit Maximum number of iterations for optimization. Defaults to
#'   10000.
#' @param abstol Absolute convergence tolerance for optimization. Defaults to
#'   1e-08.
#' @return A list with class \code{"el_aov"}.
#' @references Owen, Art. 1991. “Empirical Likelihood for Linear Models.”
#'   The Annals of Statistics 19 (4): 1725–47. \doi{10.1214/aos/1176348368}.
#' @seealso \link{el_lm}, \link{el_test}
#' @examples
#' \dontrun{
#' data("clothianidin")
#' el_aov(clo ~ trt, clothianidin)}
#' @importFrom stats .getXlevels setNames terms
#' @export
el_aov <- function(formula, data, maxit = 1e04, abstol = 1e-8) {
  .Deprecated("el_lm")
  ## check formula
  f <- attributes(terms(formula))
  if (any(
    # response required & no arbitrary manipulation on intercept
    f$response == 0, f$intercept == 0,
    length(f$variables) != 3,
    # no other formula
    typeof(f$variables[[3]]) != "symbol" ||
      length(f$variables[[3]]) != 1
  )) {
    stop("invalied model formula. specify formula as 'response ~ treatment'")
  }

  ## extract model frame
  mf <- cl <- match.call()
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(model.frame)
  mf <- eval(mf, parent.frame())

  ## type conversion
  # response
  mf[[1L]] <- as.numeric(mf[[1L]])
  # treatment
  mf[[2L]] <- as.factor(mf[[2L]])

  ## extract model terms & levels
  # model terms
  mt <- attr(mf, "terms")
  # number of levels
  p <- nlevels(mf[[2L]])
  if (p < 2L) {
    stop("contrasts can be applied only to factors with 2 or more levels")
  }
  # levels
  lv <- .getXlevels(mt, mf)
  # name for coefficients
  nm <- paste0(names(lv), lv[[1L]])

  ## construct a general block design
  # incidence matrix
  c <- unclass(table(
    factor(row.names(mf), levels = unique(row.names(mf))),
    mf[[2L]]
  ))
  # model matrix
  x <- mf[[1L]] * c

  ## specify hypothesis
  lhs <- matrix(0, nrow = p - 1, ncol = p)
  if (p == 2L) {
    lhs <- matrix(c(1, -1), nrow = 1L)
  } else {
    diag(lhs) <- 1
    diag(lhs[, -1L]) <- -1
  }
  rhs <- rep(0, p - 1)

  ## test hypothesis
  out <- ELtest(x, c, lhs, rhs, threshold = 500, maxit, abstol)
  out$coefficients <- setNames(out$coefficients, nm)
  out$xlevels <- lv
  out$call <- cl
  out$terms <- mt
  out$model <- list(model.matrix = x, incidence.matrix = c)
  if (!out$optim$convergence) {
    warning("convergence failed\n")
  }
  class(out) <- "el_aov"
  out
}

#' Tests single hypothesis for general block designs
#'
#' \emph{This function is deprecated in favor of \link{lht} and will be
#'   removed in a future release.}
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
#' @return A list of class \code{c("el_test", "melt")}.
#' @references Kim, E., MacEachern, S., and Peruggia, M., (2021),
#' "Empirical Likelihood for the Analysis of Experimental Designs,"
#' \href{https://arxiv.org/abs/2112.09206}{arxiv:2112.09206}.
#'
#' @examples
#' \dontrun{
#' # test of no treatment effect
#' data("clothianidin")
#' el_test(clo ~ trt | blk, clothianidin,
#'         lhs = matrix(c(1, -1, 0, 0,
#'                        0, 1, -1, 0,
#'                        0, 0, 1, -1), byrow = TRUE, nrow = 3))}
#' @importFrom stats reshape
#' @export
el_test <- function(formula, data, lhs, rhs = NULL, maxit = 1e04,
                    abstol = 1e-08) {
  .Deprecated("lht")
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


#' @noRd
#' @export
print.el_aov <- function(x, ...) {
  cat("\nCall:\n")
  dput(x$call, control = NULL)
  cat("\nminimizer:\n")
  cat(format(round(x$optim$par, 4), scientific = FALSE))
  cat("\n\n")
  cat("statistic:\n")
  cat(format(round(x$optim$n2logLR, 4), scientific = FALSE))
  cat("\n\n")
}
