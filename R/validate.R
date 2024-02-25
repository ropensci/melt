#' Validate `calibrate`
#'
#' Validate `calibrate` in [elt()].
#'
#' @param calibrate A single character.
#' @return A single character.
#' @noRd
validate_calibrate <- function(calibrate) {
  assert_string(calibrate)
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

#' Validate `cv`
#'
#' Validate `cv` in [confint()] and [confreg()].
#'
#' @param cv A single numeric.
#' @param th A single numeric.
#' @return A single numeric.
#' @noRd
validate_cv <- function(cv, th) {
  assert_number(cv, lower = .Machine$double.eps, finite = TRUE)
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
      if (isFALSE(any(l == c("identity", "log", "inverse")))) {
        stop(gettextf(
          "`el_glm()` does not support %s family with %s link.",
          sQuote(f), sQuote(l)
        ), domain = NA)
      }
    },
    "binomial" = {
      if (isFALSE(any(l == c("logit", "probit", "log")))) {
        stop(gettextf(
          "`el_glm()` does not support %s family with %s link.",
          sQuote(f), sQuote(l)
        ), domain = NA)
      }
    },
    "poisson" = {
      if (isFALSE(any(l == c("log", "identity", "sqrt")))) {
        stop(gettextf(
          "`el_glm()` does not support %s family with %s link.",
          sQuote(f), sQuote(l)
        ), domain = NA)
      }
    },
    "quasipoisson" = {
      if (isFALSE(any(l == c("log", "sqrt", "identity")))) {
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
#' @param lhs A numeric matrix.
#' @param p A single integer.
#' @return A numeric matrix.
#' @noRd
validate_lhs.matrix <- function(lhs, p, pnames) {
  assert_matrix(lhs,
    mode = "numeric", any.missing = FALSE, all.missing = FALSE, min.rows = 1L,
    ncols = p
  )
  assert_numeric(lhs, finite = TRUE)
  stopifnot(
    "`lhs` must have full row rank." =
      isTRUE(nrow(lhs) <= p && get_rank(lhs) == nrow(lhs))
  )
  lhs
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
  assert_numeric(lhs,
    finite = TRUE, any.missing = FALSE, all.missing = FALSE, len = p,
    typed.missing = TRUE
  )
  stopifnot(
    "`lhs` must have full row rank." = get_rank(lhs) == 1L
  )
  matrix(lhs, nrow = 1L)
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
  assert_matrix(lhs,
    mode = "numeric", any.missing = FALSE, all.missing = FALSE, min.rows = 2L,
    ncols = p
  )
  assert_numeric(lhs, finite = TRUE)
  m <- nrow(lhs)
  stopifnot(
    "Every row of `lhs` must be a nonzero vector." =
      all(apply(lhs, 1L, get_rank))
  )
  attr(lhs, "q") <- c(0L, cumsum(rep(1L, m)))
  attr(lhs, "m") <- m
  lhs
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
      test_numeric(optim$lambda,
        finite = TRUE, any.missing = FALSE, all.missing = FALSE,
        typed.missing = TRUE
      )
  )
  optim
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
#' @param rhs A numeric matrix.
#' @param p A single integer.
#' @return A numeric vector.
#' @noRd
validate_rhs.matrix <- function(rhs, p) {
  assert_matrix(rhs,
    mode = "numeric", any.missing = FALSE, all.missing = FALSE, nrows = p,
    ncols = 1L
  )
  assert_numeric(rhs, finite = TRUE)
  attr(rhs, "dim") <- NULL
  message("`rhs` is converted to a vector.")
  rhs
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
  assert_numeric(rhs,
    finite = TRUE, any.missing = FALSE, all.missing = FALSE, len = p,
    typed.missing = TRUE
  )
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
  assert_matrix(rhs,
    mode = "numeric", any.missing = FALSE, all.missing = FALSE, min.rows = 1L,
    ncols = 1L
  )
  assert_numeric(rhs, finite = TRUE)
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
  assert_numeric(rhs,
    finite = TRUE, any.missing = FALSE, all.missing = FALSE,
    typed.missing = TRUE
  )
  m <- length(rhs)
  attr(rhs, "q") <- c(0L, cumsum(rep(1L, m)))
  attr(rhs, "m") <- m
  rhs
}

#' Validate `weights`
#'
#' Validate `weights` in [el_eval()], [el_glm()], [el_lm()], [el_mean()], and
#' [el_sd()].
#'
#' @param weights An optional numeric vector.
#' @param nw A single integer.
#' @return A numeric vector.
#' @noRd
validate_weights <- function(weights, n) {
  if (is.null(weights)) {
    return(numeric(length = 0L))
  }
  assert_numeric(weights,
    lower = 0, finite = TRUE, any.missing = FALSE, all.missing = FALSE, len = n,
    typed.missing = TRUE
  )
  weights <- (n / sum(weights)) * weights
  weights
}
