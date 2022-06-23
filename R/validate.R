#' Validate `maxit`
#'
#' Validate `maxit` in [el_control()].
#'
#' @param maxit A single integer.
#' @return A single integer.
#' @noRd
validate_maxit <- function(maxit) {
  maxit <- tryCatch(as.integer(maxit),
    warning = function(w) NA,
    error = function(e) NA
  )
  stopifnot(
    "'maxit' must be a single integer" = (isTRUE(!is.na(maxit))),
    "'maxit' must be a positive single integer" = (maxit > 0L)
  )
  maxit
}

#' Validate `maxit_l`
#'
#' Validate `maxit_l` in [el_control()].
#'
#' @param maxit_l A single integer.
#' @return A single integer.
#' @noRd
validate_maxit_l <- function(maxit_l) {
  maxit_l <- tryCatch(as.integer(maxit_l),
    warning = function(w) NA,
    error = function(e) NA
  )
  stopifnot(
    "'maxit_l' must be a single integer" = (isTRUE(!is.na(maxit_l))),
    "'maxit_l' must be a positive single integer" = (maxit_l > 0L)
  )
  maxit_l
}

#' Validate `tol`
#'
#' Validate `tol` in [el_control()].
#'
#' @param tol A single numeric.
#' @return A single numeric.
#' @noRd
validate_tol <- function(tol) {
  tol <- tryCatch(as.numeric(tol),
    warning = function(w) NA,
    error = function(e) NA
  )
  stopifnot(
    "'tol' must be a single numeric" = (isTRUE(!is.na(tol))),
    "'tol' must be a finite single numeric" = (is.finite(tol)),
    "'tol' is too small" = (tol >= .Machine$double.eps)
  )
  tol
}

#' Validate `tol_l`
#'
#' Validate `tol_l` in [el_control()].
#'
#' @param tol_l A single numeric.
#' @return A single numeric.
#' @noRd
validate_tol_l <- function(tol_l) {
  tol_l <- tryCatch(as.numeric(tol_l),
    warning = function(w) NA,
    error = function(e) NA
  )
  stopifnot(
    "'tol_l' must be a single numeric" = (isTRUE(!is.na(tol_l))),
    "'tol_l' must be a finite single numeric" = (is.finite(tol_l)),
    "'tol_l' is too small" = (tol_l >= .Machine$double.eps)
  )
  tol_l
}

#' Validate `step`
#'
#' Validate `step` in [el_control()].
#'
#' @param step A single numeric.
#' @return A single numeric.
#' @noRd
validate_step <- function(step) {
  step <- tryCatch(as.numeric(step),
    warning = function(w) NA,
    error = function(e) NA
  )
  stopifnot(
    "'step' must be a single numeric" = (isTRUE(!is.na(step))),
    "'step' must be a finite single numeric" = (is.finite(step)),
    "'step' is too small" = (step >= .Machine$double.eps)
  )
  step
}

#' Validate `th`
#'
#' Validate `th` in [el_control()].
#'
#' @param th A single numeric.
#' @return A single numeric.
#' @noRd
validate_th <- function(th) {
  th <- tryCatch(as.numeric(th),
    warning = function(w) NA,
    error = function(e) NA
  )
  stopifnot(
    "'th' must be a single numeric" = (isTRUE(!is.na(th))),
    "'th' must be a finite single numeric" = (is.finite(th)),
    "'th' is too small" = (th >= .Machine$double.eps)
  )
  th
}

#' Validate `nthreads`
#'
#' Validate `nthreads` in [el_control()].
#'
#' @param nthreads A single integer.
#' @param max_threads A single integer.
#' @return A single integer.
#' @noRd
validate_nthreads <- function(nthreads, max_threads) {
  nthreads <- tryCatch(as.integer(nthreads),
    warning = function(w) NA,
    error = function(e) NA
  )
  stopifnot("'nthreads' must be a single integer" = (isTRUE(!is.na(nthreads))))
  if (nthreads < 1) {
    warning("'nthreads' is set to 1")
    nthreads <- 1L
  }
  if (nthreads > max_threads) {
    warning("'nthreads' is set to the maximum number of threads available")
    nthreads <- max_threads
  }
  nthreads
}

#' Validate `seed`
#'
#' Validate `seed` in [el_control()].
#'
#' @param seed A single integer.
#' @return A single integer.
#' @noRd
validate_seed <- function(seed) {
  seed <- tryCatch(as.integer(seed),
    warning = function(w) NA,
    error = function(e) NA
  )
  stopifnot("'seed' must be a single integer" = (isTRUE(!is.na(seed))))
  seed
}

#' Validate `B`
#'
#' Validate `B` in [el_control()].
#'
#' @param B A single integer.
#' @return A single integer.
#' @noRd
validate_B <- function(B) {
  B <- tryCatch(as.integer(B), warning = function(w) NA, error = function(e) NA)
  stopifnot(
    "'B' must be a single integer" = (isTRUE(!is.na(B))),
    "'B' must be a positive single integer" = (B > 0L)
  )
  B
}








#' @noRd
validate_model <- function(model) {
  stopifnot(
    "'model' must be a single logical" = (is.logical(model)),
    "'model' must be a single logical" = (length(model) == 1L)
  )
}




#' @noRd
check_x_ <- function(x) {
  x <- as.matrix(x)
  stopifnot(
    "'x' must have at least two observations" = (nrow(x) >= 2L),
    "'x' must must have larger number of rows than columns" =
      (nrow(x) > ncol(x)),
    "'x' must be a numeric matrix" = (is.numeric(x)),
    "'x' must be a finite numeric matrix" = (all(is.finite(x))),
    "'x' must have full column rank" = (get_rank_(x) == ncol(x))
  )
  x
}

#' @noRd
check_weights_ <- function(weights, nw) {
  if (is.null(weights)) {
    return(numeric(length = 0L))
  }
  if (!is.numeric(weights)) {
    stop("'weights' is not a numeric vector")
  }
  w <- as.vector(weights, mode = "numeric")
  if (!all(is.finite(w))) {
    stop("'weights is not a finite numeric vector")
  }
  if (any(w <= 0)) {
    stop("negative 'weights' not allowed")
  }
  if (length(w) != nw) {
    stop("length of 'weights' is incompatible with data")
  }
  w <- (nw / sum(w)) * w
  w
}











#' Validate `rhs`
#'
#' Validate `rhs` in [elt()].
#'
#' @param rhs A numeric vector or a column matrix.
#' @param p A single integer.
#' @noRd
validate_rhs <- function(rhs, p) {
  UseMethod("validate_rhs", rhs)
}

#' Validate `rhs`
#'
#' Validate `rhs` in [elt()].
#'
#' @param rhs A numeric vector.
#' @param p A single integer.
#' @noRd
validate_rhs.numeric <- function(rhs, p) {
  stopifnot("'rhs' must be a finite numeric vector" = (all(is.finite(rhs))))
  if (length(rhs) != p) {
    stop(gettextf("length of 'rhs' should be %d", p, domain = NA))
  }
}

#' Validate `rhs`
#'
#' Validate `rhs` in [elt()].
#'
#' @param rhs A numeric matrix.
#' @param p A single integer.
#' @noRd
validate_rhs.matrix <- function(rhs, p) {
  stopifnot(
    "'rhs' must be a vector" = (ncol(rhs) == 1L),
    "'rhs' must be a finite numeric vector" = (all(is.finite(rhs)))
  )
  if (nrow(rhs) != p) {
    stop(gettextf("length of 'rhs' must be %d", p, domain = NA))
  }
}

#' Validate `lhs`
#'
#' Validate `lhs` in [elt()].
#'
#' @param lhs A numeric matrix or a vector (treated as a row matrix).
#' @param p A single integer.
#' @return A numeric matrix.
#' @noRd
validate_lhs <- function(lhs, p) {
  UseMethod("validate_lhs", lhs)
}

#' Validate `lhs`
#'
#' Validate `lhs` in [elt()].
#'
#' @param lhs A numeric vector.
#' @param p A single integer.
#' @return A numeric matrix.
#' @noRd
validate_lhs.numeric <- function(lhs, p) {
  stopifnot(
    "'lhs' must be a finite numeric vector" = (all(is.finite(lhs))),
    "'object' and 'lhs' have incompatible dimensions" = (length(lhs) == p),
    "'lhs' must have full row rank" = (get_rank_(lhs) == 1L)
  )
  matrix(lhs, nrow = 1L)
}

#' Validate `lhs`
#'
#' Validate `lhs` in [elt()].
#'
#' @param lhs A numeric matrix.
#' @param p A single integer.
#' @return A numeric matrix.
#' @noRd
validate_lhs.matrix <- function(lhs, p) {
  q <- nrow(lhs)
  stopifnot(
    "'lhs' must be a finite numeric matrix" = (all(is.finite(lhs))),
    "'object' and 'lhs' have incompatible dimensions" =
      (isTRUE(q > 0L && q <= p && ncol(lhs) == p)),
    "'lhs' must have full row rank" = (get_rank_(lhs) == q)
  )
  lhs
}

#' Validate `rhs` and `lhs`
#'
#' Validate `rhs` and `lhs` in [elt()].
#'
#' @param rhs A numeric vector or a column matrix.
#' @param lhs A numeric matrix or a vector (treated as a row matrix).
#' @param p A single integer.
#' @return A list.
#' @noRd
validate_hypothesis <- function(rhs, lhs, p) {
  if (is.null(rhs) && is.null(lhs)) {
    stop("either 'rhs' or 'lhs' must be provided")
  } else if (is.null(lhs)) {
    validate_rhs(rhs, p)
  } else if (is.null(rhs)) {
    lhs <- validate_lhs(lhs, p)
    rhs <- rep(0, nrow(lhs))
  } else {
    lhs <- validate_lhs(lhs, p)
    validate_rhs(rhs, nrow(lhs))
  }
  list(l = lhs, r = rhs)
}

#' Validate `alpha`
#'
#' Validate `alpha` in [elt()].
#'
#' @param alpha A single numeric.
#' @noRd
validate_alpha <- function(alpha) {
  stopifnot(
    "'alpha' must be a finite single numeric" = (is.numeric(alpha)),
    "'alpha' must be a finite single numeric" = (length(alpha) == 1L),
    "'alpha' must be a finite single numeric" = (is.finite(alpha)),
    "'alpha' must be between 0 and 1" = (isTRUE(alpha > 0 && alpha < 1))
  )
}





#' Validate `level`
#'
#' Validate `level` in [confint()] and [confreg()].
#'
#' @param level A single numeric.
#' @noRd
validate_level <- function(level) {
  stopifnot(
    "'level' must be a finite single numeric" = (is.numeric(level)),
    "'level' must be a finite single numeric" = (length(level) == 1L),
    "'level' must be a finite single numeric" = (is.finite(level)),
    "'level' must be between 0 and 1" = (isTRUE(level >= 0 && level <= 1))
  )
}

#' Validate `cv`
#'
#' Validate `cv` in [confint()] and [confreg()].
#'
#' @param cv A single numeric.
#' @param th A single numeric.
#' @return A single numeric.
#' @noRd
validate_cv <- function(cv, th) {
  stopifnot(
    "'cv' must be a finite single numeric" = (is.numeric(cv)),
    "'cv' must be a finite single numeric" = (length(cv) == 1L),
    "'cv' must be a finite single numeric" = (is.finite(cv)),
    "'cv' is too small" = (cv >= .Machine$double.eps)
  )
  if (is.null(th)) {
    if (cv > 400) {
      stop("'cv' is too large compared to 'th'")
    }
  } else {
    if (cv > 2 * th) {
      stop("'cv' is too large compared to 'th'")
    }
  }
  cv
}




#' Validate `family`
#'
#' Validate `family` in [el_glm()].
#'
#' @param family An object of class [`family`].
#' @return A single character.
#' @noRd
validate_family <- function(family) {
  f <- family$family
  l <- family$link
  switch(f,
    "gaussian" = {
      if (!any(l == c("identity", "log", "inverse"))) {
        stop(gettextf(
          "%s family with %s link not supported by 'el_glm'",
          sQuote(f), sQuote(l)
        ), domain = NA)
      }
    },
    "binomial" = {
      if (!any(l == c("logit", "probit", "log"))) {
        stop(gettextf(
          "%s family with %s link not supported by 'el_glm'",
          sQuote(f), sQuote(l)
        ), domain = NA)
      }
    },
    # "quasibinomial" = {
    #   if (!any(l == c("logit"))) {
    #     stop(gettextf(
    #       "%s family with %s link not supported by 'el_glm'",
    #       sQuote(f), sQuote(l)
    #     ), domain = NA)
    #   }
    # },
    "poisson" = {
      if (!any(l == c("log", "identity", "sqrt"))) {
        stop(gettextf(
          "%s family with %s link not supported by 'el_glm'",
          sQuote(f), sQuote(l)
        ), domain = NA)
      }
    },
    stop(gettextf("%s family not supported by 'el_glm'", sQuote(f)),
      domain = NA
    )
  )
  paste(f, l, sep = "_")
}














#' Validate `rhs` and `lhs`
#'
#' Validate `rhs` and `lhs` in [elt()].
#'
#' @param rhs A xxx xxx
#' @param lhs A xxx xxx
#' @param p A single integer.
#' @return A list.
#' @noRd
validate_hypotheses <- function(rhs, lhs, p) {
  if (is.null(rhs) && is.null(lhs)) {
    stop("either 'rhs' or 'lhs' must be provided")
  } else if (is.null(lhs)) {
    rhs <- check_rhs2_(rhs, p)
    lhs <- matrix(rep(diag(1, nrow = p, ncol = p), attr(rhs, "m")),
      ncol = p,
      byrow = TRUE
    )
    q <- attr(rhs, "q")
    m <- attr(rhs, "m")
  } else if (is.null(rhs)) {
    lhs <- check_lhs2_(lhs, p)
    rhs <- rep(0, nrow(lhs))
    q <- attr(lhs, "q")
    m <- attr(lhs, "m")
  } else {
    rhs <- check_rhs2_(rhs, p)
    lhs <- check_lhs2_(lhs, p)
    q <- attr(lhs, "q")
    m <- attr(lhs, "m")
  }
  if (any(
    length(rhs) != nrow(lhs), isFALSE(all.equal(attr(rhs, "q"), q)),
    attr(rhs, "m") != m
  )) {
    stop("'rhs' and 'lhs' are incompatible")
  }
  list(r = rhs, l = lhs, q = q, m = m)
}

check_rhs2_ <- function(rhs, p) {
  UseMethod("check_rhs2_", rhs)
}

check_rhs2_.numeric <- function(rhs, p) {
  if (isFALSE(all(is.finite(rhs)))) {
    stop("'rhs' is not a finite numeric vector")
  }
  m <- length(rhs)
  attr(rhs, "q") <- c(0L, cumsum(rep(1L, m)))
  attr(rhs, "m") <- m
  rhs
}

check_rhs2_.list <- function(rhs, p) {
  if (isFALSE(all(vapply(rhs, is.vector, TRUE)))) {
    stop("'rhs' is not a list of vectors")
  }
  m <- length(rhs)
  if (m < 2L) {
    stop("'rhs' does not specify multiple hypotheses")
  }
  if (any(vapply(rhs, \(x) {
    !is.numeric(x) || !all(is.finite(x))
  }, FALSE))) {
    stop("'rhs' is not a list of finite numeric vectors")
  }
  if (any(vapply(rhs, \(x) {
    length(x) != p
  }, FALSE))) {
    stop(gettextf("length of every vector in 'rhs' should be %s", sQuote(p),
      domain = NA
    ))
  }
  out <- do.call(c, rhs)
  attr(out, "q") <- c(0L, cumsum(vapply(rhs, length, 1L)))
  attr(out, "m") <- m
  out
}

check_lhs2_ <- function(lhs, p) {
  UseMethod("check_lhs2_", lhs)
}

check_lhs2_.matrix <- function(lhs, p) {
  if (!is.numeric(lhs) || !all(is.finite(lhs))) {
    stop("'lhs' is not a finite numeric matrix")
  }
  if (ncol(lhs) != p) {
    stop("'object' and 'lhs' have incompatible dimensions")
  }
  m <- nrow(lhs)
  if (m < 2L) {
    stop("'lhs' does not specify multiple hypotheses")
  }
  if (isFALSE(all(apply(lhs, 1L, get_rank_)))) {
    stop("every row of 'lhs' must be a nonzero vector")
  }
  attr(lhs, "q") <- c(0L, cumsum(rep(1L, m)))
  attr(lhs, "m") <- m
  lhs
}

check_lhs2_.list <- function(lhs, p) {
  if (isFALSE(all(vapply(lhs, is.matrix, TRUE)))) {
    stop("'lhs' is not a list of matrices")
  }
  m <- length(lhs)
  if (m < 2L) {
    stop("'lhs' does not specify multiple hypotheses")
  }
  if (any(vapply(lhs, \(x) {
    !is.numeric(x) || !all(is.finite(x))
  }, FALSE))) {
    stop("'lhs' is not a list of finite numeric matrices")
  }
  if (any(vapply(lhs, \(x) {
    ncol(x) != p
  }, FALSE))) {
    stop("'object' and 'lhs' have incompatible dimensions")
  }
  if (any(vapply(lhs, \(x) {
    nrow(x) == 0L || nrow(x) > p
  }, FALSE))) {
    stop("'object' and 'lhs' have incompatible dimensions")
  }
  if (any(vapply(lhs, \(x) {
    get_rank_(x) != nrow(x)
  }, FALSE))) {
    stop("every matrix in 'lhs' must have full row rank")
  }
  out <- do.call(rbind, lhs)
  attr(out, "q") <- c(0L, cumsum(vapply(lhs, nrow, 1L)))
  attr(out, "m") <- m
  out
}
