#' Control parameters for computation
#'
#' Specifies details of computation of (constrained) empirical likelihood.
#'
#' @param maxit Maximum number of iterations for constrained minimization of
#'   empirical likelihood. Defaults to \code{200}.
#' @param maxit_l Maximum number of iterations of evaluation of empirical
#'   likelihood. Defaults to \code{25}.
#' @param tol Convergence tolerance for the constrained minimization.
#'   Defaults to \code{1e-06}.
#' @param tol_l Relative convergence tolerance for the evaluation. Defaults
#'   to \code{1e-06}.
#' @param step Step size for projected gradient method. Defaults to \code{NULL}
#'   and set to the reciprocal of sample size.
#' @param th Threshold for negative empirical log-likelihood ratio value.
#'   The iteration stops if the value exceeds the threshold. Defaults to
#'   \code{NULL}.
#' @param nthreads Number of threads for parallel computation via OpenMP (if
#'     available).
#' @return An object of class of \linkS4class{ControlEL}.
#' @seealso \link{el_eval}, \link{lht}
#' @examples
#' optcfg <- el_control(maxit = 300L, th = 200, nthreads = 1L)
#' @export
el_control <- function(maxit = 200L, maxit_l = 25L, tol = 1e-06, tol_l = 1e-06,
                       step = NULL, th = NULL, nthreads) {
  # maxit: integer (positive)
  maxit <- tryCatch(as.integer(maxit),
    warning = function(w) NA,
    error = function(e) NA
  )
  if (any(length(maxit) != 1L, is.na(maxit))) {
    stop("'maxit' is not an integer")
  }
  if (maxit < 1) {
    stop("'maxit' is not a positive integer")
  }
  # maxit_l: integer (positive)
  maxit_l <- tryCatch(as.integer(maxit_l),
    warning = function(w) NA,
    error = function(e) NA
  )
  if (any(length(maxit_l) != 1L, is.na(maxit_l))) {
    stop("'maxit_l' is not an integer")
  }
  if (maxit_l < 1) {
    stop("'maxit_l' is not a positive integer")
  }
  # tol: numeric (positive, finite)
  tol <- tryCatch(as.numeric(tol),
    warning = function(w) NA,
    error = function(e) NA
  )
  if (any(length(tol) != 1L, is.na(tol), is.infinite(tol))) {
    stop("'tol' is not a number")
  }
  if (tol < .Machine$double.eps) {
    stop("'tol' is too small")
  }
  # tol_l: numeric (positive, finite)
  tol_l <- tryCatch(as.numeric(tol_l),
    warning = function(w) NA,
    error = function(e) NA
  )
  if (any(length(tol_l) != 1L, is.na(tol_l), is.infinite(tol_l))) {
    stop("'tol' is not a number")
  }
  if (tol_l < .Machine$double.eps) {
    stop("'tol' is too small")
  }
  # step: numeric (positive, finite)
  if (!is.null(step)) {
    step <- tryCatch(as.numeric(step),
      warning = function(w) NA,
      error = function(e) NA
    )
    if (any(length(step) != 1L, is.na(step), is.infinite(step))) {
      stop("'step' is not a number")
    }
    if (step < .Machine$double.eps) {
      stop("'step' is too small")
    }
  }
  # th: numeric (positive, finite)
  if (!is.null(th)) {
    th <- tryCatch(as.numeric(th),
      warning = function(w) NA,
      error = function(e) NA
    )
    if (any(length(th) != 1L, is.na(th), is.infinite(th))) {
      stop("'th' is not a number")
    }
    if (th < .Machine$double.eps) {
      stop("'th' is too small")
    }
  }
  # nthreads: integer (positive)
  max_threads <- max_threads_()
  if (missing(nthreads)) {
    nthreads <- as.integer(max(1L, max_threads / 2L))
  } else {
    nthreads <- tryCatch(as.integer(nthreads),
      warning = function(w) NA,
      error = function(e) NA
    )
    if (any(length(nthreads) != 1L, is.na(nthreads))) {
      stop("'nthreads' is not an integer")
    }
    if (nthreads < 1) {
      warning("'nthreads' is set to 1")
      nthreads <- 1L
    }
    if (nthreads > max_threads) {
      warning("'nthreads' is set to the maximum number of threads available")
      nthreads <- max_threads
    }
  }
  new("ControlEL",
    maxit = maxit, maxit_l = maxit_l, tol = tol, tol_l = tol_l, step = step,
    th = th, nthreads = nthreads
  )
}
