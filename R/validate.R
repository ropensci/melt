#' Validate `maxit`
#'
#' Validate `maxit` in [el_control()].
#'
#' @param maxit A single integer.
#' @return A single integer.
#' @noRd
validate_maxit <- function(maxit) {
  maxit <- tryCatch(as.integer(maxit), warning = \(x) NA, error = \(x) NA)
  stopifnot(
    "`maxit` must be a finite single integer." = isTRUE(!is.na(maxit)),
    "`maxit` must be a positive single integer." = maxit > 0L
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
  maxit_l <- tryCatch(as.integer(maxit_l), warning = \(x) NA, error = \(x) NA)
  stopifnot(
    "`maxit_l` must be a finite single integer." = isTRUE(!is.na(maxit_l)),
    "`maxit_l` must be a positive single integer." = maxit_l > 0L
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
  tol <- tryCatch(as.numeric(tol), warning = \(x) NA, error = \(x) NA)
  stopifnot(
    "`tol` must be a finite single numeric." =
      isTRUE(!is.na(tol) && is.finite(tol)),
    "`tol` is too small." = tol >= .Machine$double.eps
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
  tol_l <- tryCatch(as.numeric(tol_l), warning = \(x) NA, error = \(x) NA)
  stopifnot(
    "`tol_l` must be a finite single numeric." =
      isTRUE(!is.na(tol_l) && is.finite(tol_l)),
    "`tol_l` is too small." = tol_l >= .Machine$double.eps
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
  step <- tryCatch(as.numeric(step), warning = \(x) NA, error = \(x) NA)
  stopifnot(
    "`step` must be a finite single numeric." =
      isTRUE(!is.na(step) && is.finite(step)),
    "`step` is too small." = step >= .Machine$double.eps
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
  th <- tryCatch(as.numeric(th), warning = \(x) NA, error = \(x) NA)
  stopifnot(
    "`th` must be a finite single numeric." =
      isTRUE(!is.na(th) && is.finite(th)),
    "`th` is too small." = th >= .Machine$double.eps
  )
  th
}

#' Validate `verbose`
#'
#' Validate `verbose` in [el_control()].
#'
#' @param verbose A single logical.
#' @return A single logical.
#' @noRd
validate_verbose <- function(verbose) {
  stopifnot(
    "`verbose` must be a single logical." =
      isTRUE(is.logical(verbose) && length(verbose) == 1L)
  )
  verbose
}

#' Validate `keep_data`
#'
#' Validate `keep_data` in [el_control()].
#'
#' @param keep_data A single logical.
#' @return A single logical.
#' @noRd
validate_keep_data <- function(keep_data) {
  stopifnot(
    "`keep_data` must be a single logical." =
      isTRUE(is.logical(keep_data) && length(keep_data) == 1L)
  )
  keep_data
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
  nthreads <- tryCatch(as.integer(nthreads), warning = \(x) NA, error = \(x) NA)
  stopifnot("`nthreads` must be a single integer." = isTRUE(!is.na(nthreads)))
  if (nthreads < 1) {
    warning("`nthreads` is set to 1.")
    nthreads <- 1L
  }
  if (nthreads > max_threads) {
    warning("`nthreads` is set to the maximum number of threads available.")
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
  seed <- tryCatch(as.integer(seed), warning = \(x) NA, error = \(x) NA)
  stopifnot("`seed` must be a finite single integer." = isTRUE(!is.na(seed)))
  seed
}

#' Validate `b`
#'
#' Validate `b` in [el_control()].
#'
#' @param b A single integer.
#' @return A single integer.
#' @noRd
validate_b <- function(b) {
  b <- tryCatch(as.integer(b), warning = \(x) NA, error = \(x) NA)
  stopifnot(
    "`b` must be a finite single integer." = isTRUE(!is.na(b)),
    "`b` must be a positive single integer." = b > 0L
  )
  b
}

#' Validate `m`
#'
#' Validate `m` in [el_control()].
#'
#' @param m A single integer.
#' @return A single integer.
#' @noRd
validate_m <- function(m) {
  m <- tryCatch(as.integer(m), warning = \(x) NA, error = \(x) NA)
  stopifnot(
    "`m` must be a finite single integer." = isTRUE(!is.na(m)),
    "`m` must be a positive single integer." = m > 0L
  )
  m
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
  x <- as.matrix(x, rownames.force = TRUE)
  stopifnot(
    "`x` must have at least two observations." = (nrow(x) >= 2L),
    "`x` must must have larger number of rows than columns." =
      nrow(x) > ncol(x),
    "`x` must be a finite numeric matrix." =
      isTRUE(is.numeric(x) && all(is.finite(x))),
    "`x` must have full column rank." = get_rank(x) == ncol(x)
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
    "`weights` must be a finite numeric vector." =
      isTRUE(is.numeric(weights) && all(is.finite(weights))),
    "`weights` must be all positive." = all(weights > 0)
  )
  if (length(weights) != n) {
    stop(gettextf("Length of `weights` must be %d.", n, domain = NA))
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
          "`el_glm()` does not support %s family with %s link.",
          sQuote(f), sQuote(l)
        ), domain = NA)
      }
    },
    "binomial" = {
      if (!any(l == c("logit", "probit", "log"))) {
        stop(gettextf(
          "`el_glm()` does not support %s family with %s link.",
          sQuote(f), sQuote(l)
        ), domain = NA)
      }
    },
    "poisson" = {
      if (!any(l == c("log", "identity", "sqrt"))) {
        stop(gettextf(
          "`el_glm()` does not support %s family with %s link.",
          sQuote(f), sQuote(l)
        ), domain = NA)
      }
    },
    "quasipoisson" = {
      if (!any(l == c("log", "identity"))) {
        stop(gettextf(
          "`el_glm()` does not support %s family with %s link.",
          sQuote(f), sQuote(l)
        ), domain = NA)
      }
    },
    stop(gettextf("`el_glm()` does not support %s family.", sQuote(f)),
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
    "`alpha` must be a finite single numeric." =
      isTRUE(is.numeric(alpha) && length(alpha) == 1L && is.finite(alpha)),
    "`alpha` must be between 0 and 1." = isTRUE(alpha > 0 && alpha < 1)
  )
  alpha
}

#' Validate `calibrate`
#'
#' Validate `calibrate` in [elt()].
#'
#' @param calibrate A single character.
#' @return A single character.
#' @noRd
validate_calibrate <- function(calibrate) {
  stopifnot(
    "`calibrate` must be a single character." =
      isTRUE(is.character(calibrate) && length(calibrate) == 1L)
  )
  table <- c("chisq", "boot", "f")
  calibrate <- table[pmatch(tolower(calibrate), table = table)]
  if (isTRUE(is.na(calibrate))) {
    stop(gettextf(
      "`calibrate` must be one of %s, %s, or %s.",
      dQuote("chisq"), dQuote("boot"), dQuote("f")
    ), domain = NA)
  }
  calibrate
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
    "`level` must be a finite single numeric." =
      isTRUE(is.numeric(level) && length(level) == 1L && is.finite(level)),
    "`level` must be between 0 and 1." = isTRUE(level >= 0 && level <= 1)
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
    "`cv` must be a finite single numeric." =
      isTRUE(is.numeric(cv) && length(cv) == 1L && is.finite(cv)),
    "`cv` is too small." = (cv >= .Machine$double.eps)
  )
  if (is.null(th)) {
    if (cv > 400) {
      stop("`cv` is too large compared to `th`.")
    }
  } else {
    if (cv > 2 * th) {
      stop("`cv` is too large compared to `th`.")
    }
  }
  cv
}

#' Validate `npoints`
#'
#' Validate `npoints` in [confreg()].
#'
#' @param npoints A single integer.
#' @return A single integer.
#' @noRd
validate_npoints <- function(npoints) {
  npoints <- tryCatch(as.integer(npoints), warning = \(x) NA, error = \(x) NA)
  stopifnot(
    "`npoints` must be a finite single integer." = isTRUE(!is.na(npoints)),
    "`npoints` must be a positive single integer." = npoints > 0L
  )
  npoints
}

#' Validate `rhs` and `lhs`
#'
#' Validate `rhs` and `lhs` in [elt()].
#'
#' @param rhs A numeric vector or a column matrix.
#' @param lhs A numeric matrix or a vector (treated as a row matrix).
#' @param p A single integer.
#' @param pnames An optional character vector.
#' @return A list.
#' @noRd
validate_hypothesis <- function(rhs, lhs, p, pnames) {
  if (is.null(rhs) && is.null(lhs)) {
    stop("either `rhs` or `lhs` must be provided.")
  } else if (is.null(lhs)) {
    lhs <- diag(1L, nrow = p, ncol = p)
    rhs <- validate_rhs(rhs, p)
  } else if (is.null(rhs)) {
    lhs <- validate_lhs(lhs, p, pnames)
    rhs <- rep(0, nrow(lhs))
  } else {
    lhs <- validate_lhs(lhs, p, pnames)
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
  stopifnot("`rhs` must be a finite numeric vector." = all(is.finite(rhs)))
  if (length(rhs) != p) {
    stop(gettextf("Length of `rhs` must be %d.", p, domain = NA))
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
    "`rhs` must be a finite numeric vector." =
      ncol(rhs) == 1L && all(is.finite(rhs))
  )
  if (nrow(rhs) != p) {
    stop(gettextf("Length of `rhs` must be %d.", p, domain = NA))
  }
  attr(rhs, "dim") <- NULL
  message("`rhs` is converted to a vector.")
  rhs
}

#' Validate `lhs`
#'
#' Validate `lhs` in [elt()].
#'
#' @param lhs A numeric matrix or a vector (treated as a row matrix).
#' @param p A single integer.
#' @param pnames An optional character vector.
#' @return A numeric matrix.
#' @noRd
validate_lhs <- function(lhs, p, pnames) {
  UseMethod("validate_lhs", lhs)
}

#' Validate `lhs`
#'
#' Validate `lhs` in [elt()].
#'
#' @param lhs A character vector.
#' @param p A single integer.
#' @param pnames An optional character vector.
#' @return A numeric matrix.
#' @noRd
validate_lhs.character <- function(lhs, p, pnames) {
  if (is.null(pnames)) {
    pnames <- if (p == 1L) "par" else paste0("par", seq_len(p))
  }
  q <- length(lhs)
  stopifnot(
    "Length of `lhs` must be positive." = isTRUE(q > 0L),
    "Use `rhs` for equality comparison." = isFALSE(any(grepl("=", lhs)))
  )
  out <- matrix(NA, nrow = q, ncol = p)
  for (i in seq_len(q)) {
    idx <- vapply(pnames,
      FUN = \(j) {
        grepl(j, x = lhs[i], fixed = TRUE)
      },
      FUN.VALUE = logical(1L)
    )
    sub0 <- gsub(paste(pnames, collapse = "|"), "(0)", x = lhs[i])
    eval0 <- tryCatch(eval(parse(text = sub0)),
      warning = \(x) NA, error = \(x) NA
    )
    stopifnot(
      "Invalid `lhs` specified." = isTRUE(is.finite(eval0)),
      "Constants are not allowed in `lhs`." =
        isTRUE(abs(eval0) < sqrt(.Machine$double.eps))
    )
    for (j in seq_len(p)) {
      if (idx[j]) {
        sub1 <- gsub(pnames[j], "(1)", x = lhs[i], fixed = TRUE)
        sub10 <- gsub(paste(pnames, collapse = "|"), "(0)", x = sub1)
        eval10 <- tryCatch(eval(parse(text = sub10)),
          warning = \(x) NA, error = \(x) NA
        )
        stopifnot("Invalid `lhs` specified." = isTRUE(is.finite(eval10)))
        out[i, j] <- eval10
      } else {
        out[i, j] <- 0
      }
    }
  }
  stopifnot(
    "`lhs` matrix must have full row rank." =
      isTRUE(q >= 1L && q <= p && get_rank(out) == q)
  )
  out
}

#' Validate `lhs`
#'
#' Validate `lhs` in [elt()].
#'
#' @param lhs A numeric vector.
#' @param p A single integer.
#' @param pnames An optional character vector.
#' @return A numeric matrix.
#' @noRd
validate_lhs.numeric <- function(lhs, p, pnames) {
  stopifnot(
    "`lhs` must be a finite numeric vector." = all(is.finite(lhs)),
    "`lhs` must have full row rank." = get_rank(lhs) == 1L
  )
  if (length(lhs) != p) {
    stop(gettextf("Length of `lhs` must be %d.", p, domain = NA))
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
validate_lhs.matrix <- function(lhs, p, pnames) {
  q <- nrow(lhs)
  stopifnot(
    "`lhs` must be a finite numeric matrix." = all(is.finite(lhs)),
    "`lhs` must have full row rank." =
      isTRUE(q >= 1L && q <= p && get_rank(lhs) == q)
  )
  if (ncol(lhs) != p) {
    stop(gettextf("`lhs` must have %d columns.", p, domain = NA))
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
#' @param pnames An optional character vector.
#' @return A list.
#' @noRd
validate_hypotheses <- function(rhs, lhs, p, pnames) {
  if (isTRUE(is.null(rhs) && is.null(lhs))) {
    stop("either `rhs` or `lhs` must be provided.")
  } else if (is.null(lhs)) {
    rhs <- validate_rhses(rhs, p)
    lhs <- matrix(rep(diag(1, nrow = p, ncol = p), attr(rhs, "m")),
      ncol = p,
      byrow = TRUE
    )
    q <- attr(rhs, "q")
    m <- attr(rhs, "m")
  } else if (is.null(rhs)) {
    lhs <- validate_lhses(lhs, p, pnames)
    rhs <- rep(0, nrow(lhs))
    q <- attr(lhs, "q")
    m <- attr(lhs, "m")
  } else {
    rhs <- validate_rhses(rhs, p)
    lhs <- validate_lhses(lhs, p, pnames)
    q <- attr(lhs, "q")
    m <- attr(lhs, "m")
    stopifnot(
      "`rhs` and `lhs` have incompatible dimensions." =
        isTRUE(all.equal(attr(rhs, "q"), q)) && attr(rhs, "m") == m
    )
  }
  stopifnot(
    "`rhs` and `lhs` have incompatible dimensions." = length(rhs) == nrow(lhs)
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
#' @param rhs A list of numeric vectors.
#' @param p A single integer.
#' @return A numeric vector.
#' @noRd
validate_rhses.list <- function(rhs, p) {
  m <- length(rhs)
  stopifnot(
    "`rhs` must specify multiple hypotheses." = m >= 2L,
    "`rhs` must be a list of finite numeric vectors." =
      all(vapply(rhs, FUN = is.vector, FUN.VALUE = TRUE)),
    "`rhs` must be a list of finite numeric vectors." =
      all(vapply(rhs, FUN = \(x) {
        is.numeric(x) && all(is.finite(x))
      }, FUN.VALUE = TRUE))
  )
  out <- do.call(c, rhs)
  attr(out, "q") <- c(0L, cumsum(vapply(rhs, FUN = length, FUN.VALUE = 1L)))
  attr(out, "m") <- m
  out
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
    "`rhs` must be a finite numeric vector." =
      isTRUE((ncol(rhs) == 1L) && all(is.finite(rhs)))
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
#' @param rhs A numeric vector.
#' @param p A single integer.
#' @return A numeric vector.
#' @noRd
validate_rhses.numeric <- function(rhs, p) {
  stopifnot(
    "`rhs` must be a finite numeric vector." = all(is.finite(rhs))
  )
  m <- length(rhs)
  attr(rhs, "q") <- c(0L, cumsum(rep(1L, m)))
  attr(rhs, "m") <- m
  rhs
}

#' Validate `lhs`
#'
#' Validate `lhs` in [elmt()].
#'
#' @param lhs A numeric matrix or a list of numeric matrices.
#' @param p A single integer.
#' @param pnames An optional character vector.
#' @return A numeric matrix.
#' @noRd
validate_lhses <- function(lhs, p, pnames) {
  UseMethod("validate_lhses", lhs)
}

#' Validate `lhs`
#'
#' Validate `lhs` in [elmt()].
#'
#' @param lhs A numeric matrix.
#' @param p A single integer.
#' @param pnames An optional character vector.
#' @return A numeric matrix.
#' @noRd
validate_lhses.matrix <- function(lhs, p, pnames) {
  m <- nrow(lhs)
  stopifnot(
    "`lhs` must specify multiple hypotheses." = m >= 2L,
    "`lhs` must be a finite numeric matrix." =
      isTRUE(is.numeric(lhs) && all(is.finite(lhs))),
    "Every row of `lhs` must be a nonzero vector." =
      all(apply(lhs, 1L, get_rank))
  )
  if (ncol(lhs) != p) {
    stop(gettextf("`lhs` must have %d columns.", p, domain = NA))
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
#' @param pnames An optional character vector.
#' @return A numeric matrix.
#' @noRd
validate_lhses.list <- function(lhs, p, pnames) {
  m <- length(lhs)
  stopifnot(
    "`lhs` must specify multiple hypotheses." = m >= 2L,
    "Invalid `lhs` specified." = all(vapply(lhs, FUN = \(x) {
      isTRUE(is.matrix(x) || is.character(x) || is.numeric(x))
    }, FUN.VALUE = TRUE))
  )
  lhs <- lapply(lhs, \(x) {
    validate_lhs(x, p, pnames)
  })
  out <- do.call(rbind, lhs)
  attr(out, "q") <- c(0L, cumsum(vapply(lhs, FUN = nrow, FUN.VALUE = 1L)))
  attr(out, "m") <- m
  out
}

#' Validate `optim`
#'
#' Validate `optim` in model objects.
#'
#' @param optim A list of optimization results.
#' @return A list.
#' @noRd
validate_optim <- function(optim) {
  stopifnot(
    "NaN/Inf occured during the computation." =
      isTRUE(is.numeric(optim$lambda) && all(is.finite(optim$lambda)))
  )
  optim
}
