#' Empirical likelihood hypothesis testing
#'
#' Tests single hypothesis for general block designs.
#'
#' @param formula A formula object. It must specify variables for response, treatment, and block as 'response ~ treatment | block'. Note that the use of vertical bar (|) separating treatment and block.
#' @param data A data frame containing the variables in the formula.
#' @param lhs Numeric matrix specifying linear hypothesis in terms of parameters.
#' @param rhs Optional numeric vector for the right hand side of \code{lhs}. If not specified, it is set to 0 vector.
#' @param maxit Maximum number of iterations for optimization. Defaults to 10000.
#' @param abstol Absolute convergence tolerance for optimization. Defaults to 1e-08.
#'
#' @return A list with class \code{c("el_test", "melt")}.
#' @references Kim, E., MacEachern, S., and Peruggia, M., (2021),
#' "Empirical Likelihood for the Analysis of Experimental Designs,"
#' \href{https://arxiv.org/abs/2112.09206}{arxiv:2112.09206}.
#'
#' @examples
#' ## test for equal means
#' el_test(clo ~ trt | blk, clothianidin,
#'         lhs = matrix(c(1, -1, 0, 0,
#'                        0, 1, -1, 0,
#'                        0, 0, 1, -1), byrow = TRUE, nrow = 3))
#'
#' @importFrom stats reshape
#' @export
el_test <- function(formula, data, lhs, rhs = NULL, maxit = 1e04, abstol = 1e-8) {
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
    stop("invalied model formula. specify formula as 'response ~ treatment | block'")
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

#' @export
confint.el_test <- function(object, parm, conf.level = 0.95) {
  cf <- coef(object)
  pnames <- if (is.null(names(cf))) seq(length(cf)) else names(cf)
  idx <- seq(length(cf))
  if (!missing(parm)) {
    if (is.numeric(parm)) {
      idx <- parm
      pnames <- pnames[parm]
    } else if (is.character(parm)) {
      idx <- match(parm, pnames)
      pnames <- pnames[idx]
    } else {
      stop("invalid 'parm' specified")
    }
  }
  if (!missing(conf.level) &&
      (length(conf.level) != 1L || !is.finite(conf.level) ||
       conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  if (conf.level == 0) {
    ci <- matrix(rep(cf, 2L), ncol = 2L)
  } else if (conf.level == 1) {
    p <- length(pnames)
    ci <- matrix(c(rep(-Inf, p), rep(Inf, p)), ncol = 2L)
  } else {
    cutoff <- qchisq(conf.level, 1L)
    ci <- EL_confint(object$data, object$optim$type, cf, cutoff, idx,
                     maxit = 100, abstol = 1e-8)
  }
  a <- (1 - conf.level)/2
  a <- c(a, 1 - a)
  pct <- paste(round(100 * a, 1L), "%")
  dimnames(ci) <- list(pnames, pct)
  ci
}

#' @export
print.el_test <- function(object, digits = getOption("digits"), prefix = "\t", ...) {
  cat("\n")
  cat(strwrap(object$method, prefix = prefix), sep = "\n")
  cat("\n")
  cat("data: ", object$data.name, "\n", sep = "")
  out <- character()
  if (!is.null(object$statistic))
    out <- c(out, paste("Chisq", names(object$statistic), "=",
                        format(object$statistic, digits = max(1L, digits - 2L))))
  if (!is.null(object$df))
    out <- c(out, paste("df", "=", object$df))
  if (!is.null(object$p.value)) {
    fp <- format.pval(object$p.value, digits = max(1L, digits - 3L))
    out <- c(out, paste("p-value", if (startsWith(fp, "<")) fp else paste("=",
                                                                          fp)))
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  if (!is.null(object$alternative)) {
    cat("alternative hypothesis: ")
    if (!is.null(object$null.value)) {
      if (length(object$null.value) == 1L) {
        alt.char <- switch(object$alternative, two.sided = "not equal to",
                           less = "less than", greater = "greater than")
        cat("true ", names(object$null.value), " is ", alt.char,
            " ", object$null.value, "\n", sep = "")
      }
      else {
        cat(object$alternative, "\nnull values:\n", sep = "")
        print(object$null.value, digits = digits, ...)
      }
    }
    else cat(object$alternative, "\n", sep = "")
  }
  if (!is.null(object$coefficients)) {
    cat("maximum EL estimates:\n")
    print(object$coefficients, digits = digits, ...)
  }
  cat("\n")
  invisible(x)
}
