#' Control parameters for computation
#'
#' Specifies computational details of (constrained) empirical likelihood.
#'
#' @param maxit A single integer for the maximum number of iterations for
#'   constrained minimization of empirical likelihood. Defaults to `200`.
#' @param maxit_l A single integer for the maximum number of iterations for
#'   evaluation of empirical likelihood. Defaults to `25`.
#' @param tol A single numeric for the convergence tolerance for the constrained
#'   minimization. Defaults to `1e-06`.
#' @param tol_l A single numeric for the relative convergence tolerance for the
#'   evaluation. Defaults to `1e-06`.
#' @param step A single numeric for the step size for projected gradient descent
#'   method. Defaults to `NULL` and sets the step size to the reciprocal of the
#'   sample size.
#' @param th A single numeric for the threshold for the negative empirical
#'   log-likelihood ratio. The iteration stops if the value exceeds the
#'   threshold. Defaults to `NULL` and sets the threshold to `200 * d`, where
#'   `d` corresponds to the degrees of freedom of the limiting chi-squared
#'   distribution of the statistic.
#' @param verbose A single logical. If `TRUE`, a message on the convergence
#'   status is printed when fitting objects that inherit from class
#'   \linkS4class{EL}.
#'   Defaults to `FALSE`.
#' @param keep_data A single logical. If `TRUE`, the data used for fitting
#'   objects that inherit from class \linkS4class{EL} are stored for later use
#'   with other methods. Defaults to `TRUE`.
#' @param nthreads A single integer for the number of threads for parallel
#'   computation via OpenMP (if available). Defaults to half the available
#'   threads. For better performance, it is generally recommended in most
#'   platforms to limit the number of threads to the number of physical cores.
#'   Note that it applies to the following functions that involve multiple
#'   evaluations or optimizations: [confint()], [confreg()], [el_lm()],
#'   [el_glm()], [eld()], and [elt()].
#' @param seed A single integer for the seed for random number generation. It
#'   only applies to [elt()] when `calibrate` is set to `"boot"`. Defaults to
#'   `NULL`. In this case, a seed is set to a random integer generated from 1 to
#'   the maximum integer supported by \R on the machine, which is determined by
#'   [set.seed()]. Only one seed is needed even when multiple threads are used
#'   with `nthreads`. Each thread is given a separate seed to produce a
#'   non-overlapping but reproducible sequence of random numbers. The
#'   Xoshiro256+ pseudo-random number generator is used internally to work with
#'   OpenMP.
#' @param b A single integer for the number of bootstrap replicates. It only
#'   applies to [elt()] when `calibrate` is set to `"boot"`. Defaults to
#'   `10000`.
#' @param m A single integer for the number of Monte Carlo samples. It only
#'   applies to [elmt()]. Defaults to `1e+06`.
#' @return An object of class of \linkS4class{ControlEL}.
#' @seealso [el_eval()], [elt()]
#' @examples
#' optcfg <- el_control(maxit = 300, step = 0.01, th = 200, nthreads = 1)
#' @export
el_control <- function(maxit = 200L,
                       maxit_l = 25L,
                       tol = 1e-06,
                       tol_l = 1e-06,
                       step = NULL,
                       th = NULL,
                       verbose = FALSE,
                       keep_data = TRUE,
                       nthreads,
                       seed = NULL,
                       b = 10000L,
                       m = 1e+06L) {
  maxit <- assert_int(maxit, lower = 1L, coerce = TRUE)
  maxit_l <- assert_int(maxit_l, lower = 1L, coerce = TRUE)
  tol <- assert_number(tol, lower = .Machine$double.eps, finite = TRUE)
  tol_l <- assert_number(tol_l, lower = .Machine$double.eps, finite = TRUE)
  if (isFALSE(is.null(step))) {
    step <- assert_number(step, lower = .Machine$double.eps, finite = TRUE)
  }
  if (isFALSE(is.null(th))) {
    th <- assert_number(th, lower = .Machine$double.eps, finite = TRUE)
  }
  verbose <- assert_logical(verbose,
    any.missing = FALSE, all.missing = FALSE, len = 1L, typed.missing = TRUE
  )
  keep_data <- assert_logical(keep_data,
    any.missing = FALSE, all.missing = FALSE, len = 1L, typed.missing = TRUE
  )
  max_threads <- get_max_threads()
  if (missing(nthreads)) {
    nthreads <- as.integer(max(1L, max_threads / 2L))
  } else {
    nthreads <- assert_int(nthreads, lower = 1L, coerce = TRUE)
    if (nthreads > max_threads) {
      warning("`nthreads` is set to the maximum number of threads available.")
      nthreads <- max_threads
    }
  }
  if (isFALSE(is.null(seed))) {
    seed <- assert_int(seed, coerce = TRUE)
  }
  b <- assert_int(b, lower = 1L, coerce = TRUE)
  m <- assert_int(m, lower = 1L, coerce = TRUE)
  new("ControlEL",
    maxit = maxit, maxit_l = maxit_l, tol = tol, tol_l = tol_l, step = step,
    th = th, verbose = verbose, keep_data = keep_data, nthreads = nthreads,
    seed = seed, b = b, m = m
  )
}
