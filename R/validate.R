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

#' Validate `model`
#'
#' Validate `model` in [el_glm()], [el_lm()], and [el_mean()].
#'
#' @param model A single logical.
#' @return A single logical.
#' @noRd
validate_model <- function(model) {
  stopifnot(
    "'model' must be a single logical" =
      (isTRUE(is.logical(model) && length(model) == 1L))
  )
  model
}

#' Validate `x`
#'
#' Validate `x` in [el_mean()].
#'
#' @param x A numeric matrix, or an object that can be coerced to a numeric
#'   matrix.
#' @return A numeric matrix.
#' @noRd
validate_x <- function(x) {
  x <- as.matrix(x)
  stopifnot(
    "'x' must have at least two observations" = (nrow(x) >= 2L),
    "'x' must must have larger number of rows than columns" =
      (nrow(x) > ncol(x)),
    "'x' must be a finite numeric matrix" =
      (isTRUE(is.numeric(x) && all(is.finite(x)))),
    "'x' must have full column rank" = (getRank(x) == ncol(x))
  )
  x
}

#' Validate `weights`
#'
#' Validate `weights` in [el_eval()], [el_glm()], [el_lm()], and [el_mean()].
#'
#' @param weights An optional numeric vector.
#' @param nw A single integer.
#' @return A numeric vector.
#' @noRd
validate_weights <- function(weights, n) {
  if (is.null(weights)) {
    return(numeric(length = 0L))
  }
  stopifnot(
    "'weights' must be a finite numeric vector" =
      (isTRUE(is.numeric(weights) && all(is.finite(weights)))),
    "'weights' must be positive" = (all(weights > 0))
  )
  if (length(weights) != n) {
    stop(gettextf("length of 'weights' must be %d", n, domain = NA))
  }
  weights <- (n / sum(weights)) * weights
  weights
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

#' Validate `alpha`
#'
#' Validate `alpha` in [elt()].
#'
#' @param alpha A single numeric.
#' @return A single numeric.
#' @noRd
validate_alpha <- function(alpha) {
  stopifnot(
    "'alpha' must be a finite single numeric" =
      (isTRUE(is.numeric(alpha) && length(alpha) == 1L && is.finite(alpha))),
    "'alpha' must be between 0 and 1" = (isTRUE(alpha > 0 && alpha < 1))
  )
  alpha
}

#' Validate `level`
#'
#' Validate `level` in [confint()] and [confreg()].
#'
#' @param level A single numeric.
#' @return A single numeric.
#' @noRd
validate_level <- function(level) {
  stopifnot(
    "'level' must be a finite single numeric" =
      (isTRUE(is.numeric(level) && length(level) == 1L && is.finite(level))),
    "'level' must be between 0 and 1" = (isTRUE(level >= 0 && level <= 1))
  )
  level
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
    "'cv' must be a finite single numeric" =
      (isTRUE(is.numeric(cv) && length(cv) == 1L && is.finite(cv))),
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
    rhs <- validate_rhs(rhs, p)
  } else if (is.null(rhs)) {
    lhs <- validate_lhs(lhs, p)
    rhs <- rep(0, nrow(lhs))
  } else {
    lhs <- validate_lhs(lhs, p)
    rhs <- validate_rhs(rhs, nrow(lhs))
  }
  list(l = lhs, r = rhs)
}

#' Validate `rhs`
#'
#' Validate `rhs` in [elt()].
#'
#' @param rhs A numeric vector or a column matrix.
#' @param p A single integer.
#' @return A numeric vector.
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
#' @return A numeric vector.
#' @noRd
validate_rhs.numeric <- function(rhs, p) {
  stopifnot("'rhs' must be a finite numeric vector" = (all(is.finite(rhs))))
  if (length(rhs) != p) {
    stop(gettextf("length of 'rhs' must be %d", p, domain = NA))
  }
  rhs
}

#' Validate `rhs`
#'
#' Validate `rhs` in [elt()].
#'
#' @param rhs A numeric matrix.
#' @param p A single integer.
#' @return A numeric vector.
#' @noRd
validate_rhs.matrix <- function(rhs, p) {
  stopifnot(
    "'rhs' must be a finite numeric vector" =
      (ncol(rhs) == 1L && all(is.finite(rhs)))
  )
  if (nrow(rhs) != p) {
    stop(gettextf("length of 'rhs' must be %d", p, domain = NA))
  }
  attr(rhs, "dim") <- NULL
  rhs
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
    "'lhs' must have full row rank" = (getRank(lhs) == 1L)
  )
  if (length(lhs) != p) {
    stop(gettextf("length of 'lhs' must be %d", p, domain = NA))
  }
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
    "'lhs' must have full row rank" =
      (isTRUE(q >= 1L && q <= p && getRank(lhs) == q))
  )
  if (ncol(lhs) != p) {
    stop(gettextf("'lhs' must have %d columns", p, domain = NA))
  }
  lhs
}

#' Validate `rhs` and `lhs`
#'
#' Validate `rhs` and `lhs` in [elmt()].
#'
#' @param rhs A numeric vector (column matrix) or a list of numeric vectors.
#' @param lhs A numeric matrix or a list of numeric matrices.
#' @param p A single integer.
#' @return A list.
#' @noRd
validate_hypotheses <- function(rhs, lhs, p) {
  if (isTRUE(is.null(rhs) && is.null(lhs))) {
    stop("either 'rhs' or 'lhs' must be provided")
  } else if (is.null(lhs)) {
    rhs <- validate_rhses(rhs, p)
    lhs <- matrix(rep(diag(1, nrow = p, ncol = p), attr(rhs, "m")),
      ncol = p,
      byrow = TRUE
    )
    q <- attr(rhs, "q")
    m <- attr(rhs, "m")
  } else if (is.null(rhs)) {
    lhs <- validate_lhses(lhs, p)
    rhs <- rep(0, nrow(lhs))
    q <- attr(lhs, "q")
    m <- attr(lhs, "m")
  } else {
    rhs <- validate_rhses(rhs, p)
    lhs <- validate_lhses(lhs, p)
    q <- attr(lhs, "q")
    m <- attr(lhs, "m")
    stopifnot(
      "'rhs' and 'lhs' have incompatible dimensions" =
        ((isTRUE(all.equal(attr(rhs, "q"), q)) && attr(rhs, "m") == m))
    )
  }
  stopifnot(
    "'rhs' and 'lhs' have incompatible dimensions" =
      (length(rhs) == nrow(lhs))
  )
  list(r = rhs, l = lhs, q = q, m = m)
}

#' Validate `rhs`
#'
#' Validate `rhs` in [elmt()].
#'
#' @param rhs A numeric vector (column matrix) or a list of numeric vectors.
#' @param p A single integer.
#' @return A numeric vector.
#' @noRd
validate_rhses <- function(rhs, p) {
  UseMethod("validate_rhses", rhs)
}

#' Validate `rhs`
#'
#' Validate `rhs` in [elmt()].
#'
#' @param rhs A numeric vector.
#' @param p A single integer.
#' @return A numeric vector.
#' @noRd
validate_rhses.numeric <- function(rhs, p) {
  stopifnot(
    "'rhs' must be a finite numeric vector" = (all(is.finite(rhs)))
  )
  m <- length(rhs)
  attr(rhs, "q") <- c(0L, cumsum(rep(1L, m)))
  attr(rhs, "m") <- m
  rhs
}

#' Validate `rhs`
#'
#' Validate `rhs` in [elmt()].
#'
#' @param rhs A numeric matrix.
#' @param p A single integer.
#' @return A numeric vector.
#' @noRd
validate_rhses.matrix <- function(rhs, p) {
  stopifnot(
    "'rhs' must be a finite numeric vector" =
      (isTRUE((ncol(rhs) == 1L) && all(is.finite(rhs))))
  )
  attr(rhs, "dim") <- NULL
  m <- length(rhs)
  attr(rhs, "q") <- c(0L, cumsum(rep(1L, m)))
  attr(rhs, "m") <- m
  rhs
}

#' Validate `rhs`
#'
#' Validate `rhs` in [elmt()].
#'
#' @param rhs A list of numeric vectors.
#' @param p A single integer.
#' @return A numeric vector.
#' @noRd
validate_rhses.list <- function(rhs, p) {
  m <- length(rhs)
  stopifnot(
    "'rhs' must specify multiple hypotheses" = (m >= 2L),
    "'rhs' must be a list of finite numeric vectors" =
      (all(vapply(rhs, is.vector, TRUE))),
    "'rhs' must be a list of finite numeric vectors" =
      (all(vapply(rhs, \(x) {
        is.numeric(x) && all(is.finite(x))
      }, TRUE)))
  )
  out <- do.call(c, rhs)
  attr(out, "q") <- c(0L, cumsum(vapply(rhs, length, 1L)))
  attr(out, "m") <- m
  out
}

#' Validate `lhs`
#'
#' Validate `lhs` in [elmt()].
#'
#' @param lhs A numeric matrix or a list of numeric matrices.
#' @param p A single integer.
#' @return A numeric matrix.
#' @noRd
validate_lhses <- function(lhs, p) {
  UseMethod("validate_lhses", lhs)
}

#' Validate `lhs`
#'
#' Validate `lhs` in [elmt()].
#'
#' @param lhs A numeric matrix.
#' @param p A single integer.
#' @return A numeric matrix.
#' @noRd
validate_lhses.matrix <- function(lhs, p) {
  m <- nrow(lhs)
  stopifnot(
    "'lhs' must specify multiple hypotheses" = (m >= 2L),
    "'lhs' must be a finite numeric matrix" =
      (isTRUE(is.numeric(lhs) && all(is.finite(lhs)))),
    "every row of 'lhs' must be a nonzero vector" =
      (all(apply(lhs, 1L, getRank)))
  )
  if (ncol(lhs) != p) {
    stop(gettextf("'lhs' must have %d columns", p, domain = NA))
  }
  attr(lhs, "q") <- c(0L, cumsum(rep(1L, m)))
  attr(lhs, "m") <- m
  lhs
}

#' Validate `lhs`
#'
#' Validate `lhs` in [elmt()].
#'
#' @param lhs A list of numeric matrices.
#' @param p A single integer.
#' @return A numeric matrix.
#' @noRd
validate_lhses.list <- function(lhs, p) {
  m <- length(lhs)
  stopifnot(
    "'lhs' must specify multiple hypotheses" = (m >= 2L),
    "'lhs' must be a list of finite numeric matrices" =
      (all(vapply(lhs, is.matrix, TRUE))),
    "'lhs' must be a list of finite numeric matrices" =
      (all(vapply(lhs, \(x) {
        is.numeric(x) && all(is.finite(x))
      }, TRUE))),
    "every matrix in 'lhs' must have full row rank" =
      (all(vapply(lhs, \(x) {
        nrow(x) >= 1L && nrow(x) <= p && getRank(x) == nrow(x)
      }, TRUE)))
  )
  if (any(vapply(lhs, \(x) {
    ncol(x) != p
  }, FALSE))) {
    stop(gettextf("every matrix in 'lhs' must have %d columns", p, domain = NA))
  }
  out <- do.call(rbind, lhs)
  attr(out, "q") <- c(0L, cumsum(vapply(lhs, nrow, 1L)))
  attr(out, "m") <- m
  out
}
