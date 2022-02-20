#' Empirical likelihood pairwise comparisons
#'
#' Tests all pairwise comparisons or comparisons with control for general block designs. Two single step asymptotic \eqn{k}-FWER (generalized family-wise error rate) controlling procedures are available: asymptotic Monte Carlo (AMC) and nonparametric bootstrap (NB).
#'
#' @param formula A formula object. It must specify variables for response, treatment, and block as 'response ~ treatment | block'. Note that the use of vertical bar (|) separating treatment and block.
#' @param data A data frame containing the variables in the formula.
#' @param control Optional character specifying the treatment for comparisons with control.
#' @param k Single integer for \eqn{k} in \eqn{k}-FWER. Defaults to 1.
#' @param alpha Level of the test. Defaults to 0.05.
#' @param method Single character for the procedure to be used; either 'AMC' or 'NB' is supported. Defaults to 'AMC'.
#' @param B Number of Monte Carlo samples for the AMC (number of bootstrap replicates for the NB).
#' @param nthread Number of threads to be used with OpenMP. Only applied when the NB is chosen as the method. Defaults to 1.
#' @param progress If TRUE, will show progress status.
#' @param maxit Maximum number of iterations for optimization. Defaults to 10000.
#' @param abstol Absolute convergence tolerance for optimization. Defaults to 1e-08.
#' @return A list with class \code{c("pairwise", "melt")}.
#' @references Kim, E., MacEachern, S., and Peruggia, M., (2021),
#' "Empirical Likelihood for the Analysis of Experimental Designs,"
#' \href{https://arxiv.org/abs/2112.09206}{arxiv:2112.09206}.
#' @examples
#' ## all pairwise comparisons
#' data("clothianidin")
#' el_pairwise(clo ~ trt | blk, clothianidin, B = 10000)
#'
#' ## comparisons with control
#' el_pairwise(clo ~ trt | blk, clothianidin, control = "Naked", method = "NB", B = 1000)
#'
#' @importFrom stats terms
#' @export
el_pairwise <- function(formula, data, control = NULL, k = 1, alpha = 0.05,
                        method = c("AMC", "NB"), B, nthread = 1, progress = TRUE,
                        maxit = 10000, abstol = 1e-08) {
  ## check method
  method <- match.arg(method)

  ## check formula
  f <- attributes(terms(formula))
  if (any(
    # response required & no arbitrary manipulation on intercept
    f$response == 0, f$intercept == 0,
    length(f$variables) != 3,
    # no other formula
    typeof(f$variables[[3]]) != "language" ||
    length(f$variables[[3]]) != 3,
    # "|" operator needed
    f$variables[[3]][[1]] != "|",
    # no transformation of variables
    typeof(f$variables[[3]][[2]]) != "symbol" ||
    typeof(f$variables[[3]][[3]]) != "symbol",
    # distinct variables for treatment and block
    f$variables[[3]][[2]] == f$variables[[3]][[3]])
  ) {
    stop("invalied model formula. specify formula as 'response ~ treatment | block'.")
  }

  ## pseudo formula for model.frame
  lhs <- f$variables[[2]]
  rhs <- c(f$variables[[3]][[2]], f$variables[[3]][[3]])
  pf <- formula(paste(lhs, paste(rhs, collapse = " + "), sep = " ~ "))

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
    stop("number of blocks should be larger than number of treatments.")
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

  ## check whether all pairwise comparisons or comparisons to control
  match.arg(control, gbd$trt)
  if (is.null(control)) {
    ctrl <- 0
  } else {
    ctrl <- which(control == gbd$trt)
  }

  ## pairwise comparisons
  out <- pairwise(gbd$model_matrix, gbd$incidence_matrix,
                  control = ctrl, k, alpha, interval = TRUE,
                  method, B, nthread, progress,
                  threshold = 50, maxit, abstol)
  out$trt <- gbd$trt
  out$control <- control
  out$model.matrix <- gbd$model_matrix
  out$incidence.matrix <- gbd$incidence_matrix
  if (!all(out$convergence)) {
    if (method == "NB") {
      warning("convergence failed. switched to AMC for confidence intervals.\n")
    } else {
      warning("convergence failed.\n")
    }
  }
  out
}
