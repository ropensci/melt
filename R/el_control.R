#' Control parameters for computation
#'
#' Specifies computational details of (constrained) empirical likelihood.
#'
#' @param maxit A single integer for the maximum number of iterations for
#'   constrained minimization of empirical likelihood. Defaults to \code{200}.
#' @param maxit_l A single integer for the maximum number of iterations for
#'   evaluation of empirical likelihood. Defaults to \code{25}.
#' @param tol A single numeric for the convergence tolerance for the constrained
#'   minimization. Defaults to \code{1e-06}.
#' @param tol_l A single numeric for the relative convergence tolerance for the
#'   evaluation. Defaults to \code{1e-06}.
#' @param step A single numeric for the step size for projected gradient descent
#'   method. Defaults to \code{NULL} and set to the reciprocal of sample size.
#' @param th A single numeric for the threshold for the negative empirical
#'   log-likelihood ratio. The iteration stops if the value exceeds the
#'   threshold. Defaults to \code{NULL} and sets the threshold to
#'   \code{200 * d}, where \code{d} corresponds to the degrees of freedom of the
#'   limiting chi-squared distribution of the statistic.
#' @param nthreads A single integer for the number of threads for parallel
#'   computation via OpenMP (if available). Defaults to the half of the
#'   available threads. For better performance, it is generally recommended to
#'   limit the number of threads to the number of physical cores. Note that it
#'   only applies to the following functions that involve multiple evaluations
#'   or minimizations: \code{\link{confint}}, \code{\link{confreg}},
#'   \code{\link{el_lm}}, \code{\link{el_glm}}, \code{\link{eld}}, and
#'   \code{\link{elt}}.
#' @param seed A single integer for the seed for random number generation. It
#'   only applies to \code{\link{elt}} when \code{calibrate} is set to
#'   \code{"boot"}. Defaults to a random integer generated from 1 to the maximum
#'   integer supported by \R on the machine, which is determined by
#'   \code{\link[base]{set.seed}}. Only one seed is needed even when multiple
#'   threads are used with \code{nthreads}. Each thread is given a separate seed
#'   to produce a non-overlapping but reproducible sequence of random numbers.
#'   The \code{xoshiro256+} pseudo-random number generator is used internally to
#'   work with OpenMP.
#' @param B A single integer for the number of bootstrap replicates. It only
#'   applies to \code{\link{elt}} when \code{calibrate} is set to \code{"boot"}.
#'   Defaults to \code{10000L}.
#' @return An object of class of \linkS4class{ControlEL}.
#' @seealso \link{el_eval}, \link{elt}
#' @examples
#' optcfg <- el_control(maxit = 300, th = 200, nthreads = 1)
#' @export
el_control <- function(maxit = 200L,
                       maxit_l = 25L,
                       tol = 1e-06,
                       tol_l = 1e-06,
                       step = NULL,
                       th = NULL,
                       nthreads,
                       seed = sample.int(.Machine$integer.max, 1L),
                       B = 10000L) {
  # maxit: single integer (positive)
  maxit <- tryCatch(as.integer(maxit),
    warning = function(w) NA,
    error = function(e) NA
  )
  stopifnot(
    "'maxit' must be a single integer" = (isTRUE(!is.na(maxit))),
    "'maxit' must be a positive single integer" = (maxit > 0L)
  )
  # maxit_l: single integer (positive)
  maxit_l <- tryCatch(as.integer(maxit_l),
    warning = function(w) NA,
    error = function(e) NA
  )
  stopifnot(
    "'maxit_l' must be a single integer" = (isTRUE(!is.na(maxit_l))),
    "'maxit_l' must be a positive single integer" = (maxit_l > 0L)
  )
  # tol: single numeric (positive, finite)
  tol <- tryCatch(as.numeric(tol),
    warning = function(w) NA,
    error = function(e) NA
  )
  stopifnot(
    "'tol' must be a single numeric" = (isTRUE(!is.na(tol))),
    "'tol' must be a finite single numeric" = (is.finite(tol)),
    "'tol' is too small" = (tol >= .Machine$double.eps)
  )
  # tol_l: single numeric (positive, finite)
  tol_l <- tryCatch(as.numeric(tol_l),
    warning = function(w) NA,
    error = function(e) NA
  )
  stopifnot(
    "'tol_l' must be a single numeric" = (isTRUE(!is.na(tol_l))),
    "'tol_l' must be a finite single numeric" = (is.finite(tol_l)),
    "'tol_l' is too small" = (tol_l >= .Machine$double.eps)
  )
  # step: single numeric (positive, finite)
  if (!is.null(step)) {
    step <- tryCatch(as.numeric(step),
      warning = function(w) NA,
      error = function(e) NA
    )
    stopifnot(
      "'step' must be a single numeric" = (isTRUE(!is.na(step))),
      "'step' must be a finite single numeric" = (is.finite(step)),
      "'step' is too small" = (step >= .Machine$double.eps)
    )
  }
  # th: single numeric (positive, finite)
  if (!is.null(th)) {
    th <- tryCatch(as.numeric(th),
      warning = function(w) NA,
      error = function(e) NA
    )
    stopifnot(
      "'th' must be a single numeric" = (isTRUE(!is.na(th))),
      "'th' must be a finite single numeric" = (is.finite(th)),
      "'th' is too small" = (th >= .Machine$double.eps)
    )
  }
  # nthreads: single integer (positive, finite)
  max_threads <- max_threads_()
  if (missing(nthreads)) {
    nthreads <- as.integer(max(1L, max_threads / 2L))
  } else {
    nthreads <- tryCatch(as.integer(nthreads),
      warning = function(w) NA,
      error = function(e) NA
    )
    stopifnot(
      "'nthreads' must be a single integer" = (isTRUE(!is.na(nthreads)))
    )
    if (nthreads < 1) {
      warning("'nthreads' is set to 1")
      nthreads <- 1L
    }
    if (nthreads > max_threads) {
      warning("'nthreads' is set to the maximum number of threads available")
      nthreads <- max_threads
    }
  }
  # seed: single integer (finite)
  seed <- tryCatch(as.integer(seed),
    warning = function(w) NA,
    error = function(e) NA
  )
  stopifnot("'seed' must be a single integer" = (isTRUE(!is.na(seed))))
  # B: single integer (positive, finite)
  B <- tryCatch(as.integer(B), warning = function(w) NA, error = function(e) NA)
  stopifnot(
    "'B' must be a single integer" = (isTRUE(!is.na(B))),
    "'B' must be a positive single integer" = (B > 0L)
  )
  new("ControlEL",
    maxit = maxit, maxit_l = maxit_l, tol = tol, tol_l = tol_l, step = step,
    th = th, nthreads = nthreads, seed = seed, B = B
  )
}
