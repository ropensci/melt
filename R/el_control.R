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
#'   method. Defaults to `NULL` and set to the reciprocal of sample size.
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
#'   only applies to [elt()] when `calibrate` is set to `"boot"`. Defaults to a
#'   random integer generated from 1 to the maximum integer supported by \R on
#'   the machine, which is determined by [set.seed()]. Only one seed is
#'   needed even when multiple threads are used with `nthreads`. Each thread is
#'   given a separate seed to produce a non-overlapping but reproducible
#'   sequence of random numbers. The Xoshiro256+ pseudo-random number generator
#'   is used internally to work with OpenMP.
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
#' @srrstats {G2.1, G2.1a} Assertions on types of inputs are clarified
#'   throughout the package documentation.
#' @srrstats {RE3.2, RE3.3} Convergence thresholds can be set by the `tol` and
#'   `tol_l` arguments with the default values documented.
#' @srrstats {RE3.0, RE3.1} The `verbose` argument controls whether to print a
#'   message on the convergence status when fitting objects. The status is
#'   printed regardless of the `verbose` when `print()` method is used to the
#'   fitted objects. Alternatively, `conv()` method also extracts the status.
el_control <- function(maxit = 200L,
                       maxit_l = 25L,
                       tol = 1e-06,
                       tol_l = 1e-06,
                       step = NULL,
                       th = NULL,
                       verbose = FALSE,
                       keep_data = TRUE,
                       nthreads,
                       seed = sample.int(.Machine$integer.max, 1L),
                       b = 10000L,
                       m = 1e+06L) {
  maxit <- validate_maxit(maxit)
  maxit_l <- validate_maxit_l(maxit_l)
  tol <- validate_tol(tol)
  tol_l <- validate_tol_l(tol_l)
  if (!is.null(step)) {
    step <- validate_step(step)
  }
  if (!is.null(th)) {
    th <- validate_th(th)
  }
  verbose <- validate_verbose(verbose)
  keep_data <- validate_keep_data(keep_data)
  max_threads <- get_max_threads()
  if (missing(nthreads)) {
    nthreads <- as.integer(max(1L, max_threads / 2L))
  } else {
    nthreads <- validate_nthreads(nthreads, max_threads)
  }
  seed <- validate_seed(seed)
  b <- validate_b(b)
  m <- validate_m(m)
  new("ControlEL",
    maxit = maxit, maxit_l = maxit_l, tol = tol, tol_l = tol_l, step = step,
    th = th, verbose = verbose, keep_data = keep_data, nthreads = nthreads,
    seed = seed, b = b, m = m
  )
}
