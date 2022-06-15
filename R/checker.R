check_weights <- function(weights, nw) {
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

check_rhs_ <- function(rhs, p) {
  UseMethod("check_rhs_", rhs)
}

check_rhs_.numeric <- function(rhs, p) {
  if (isFALSE(all(is.finite(rhs)))) {
    stop("'rhs' is not a finite numeric vector")
  }
  if (length(rhs) != p) {
    stop(gettextf("length of 'rhs' should be %d", p, domain = NA))
  }
  rhs
}

check_rhs_.matrix <- function(rhs, p) {
  if (ncol(rhs) == 1L) {
    stop(gettextf("'rhs' is not a vector"))
  }
  if (isFALSE(all(is.finite(rhs)))) {
    stop("'rhs' is not a finite numeric vector")
  }
  if (nrow(rhs) != p) {
    stop(gettextf("length of 'rhs' should be %d", p, domain = NA))
  }
  rhs
}

check_lhs_ <- function(lhs, p) {
  UseMethod("check_lhs_", lhs)
}

check_lhs_.numeric <- function(lhs, p) {
  if (isFALSE(all(is.finite(lhs)))) {
    stop("'lhs' is not a finite numeric matrix")
  }
  dim(lhs) <- c(1L, length(lhs))
  if (ncol(lhs) != p) {
    stop("'object' and 'lhs' have incompatible dimensions")
  }
  q <- nrow(lhs)
  if (q == 0L || q > p) {
    stop("'object' and 'lhs' have incompatible dimensions")
  }
  if (get_rank_(lhs) != q) {
    stop("'lhs' must have full row rank")
  }
  lhs
}

check_lhs_.matrix <- function(lhs, p) {
  if (isFALSE(all(is.finite(lhs)))) {
    stop("'lhs' is not a finite numeric matrix")
  }
  if (ncol(lhs) != p) {
    stop("'object' and 'lhs' have incompatible dimensions")
  }
  q <- nrow(lhs)
  if (q == 0L || q > p) {
    stop("'object' and 'lhs' have incompatible dimensions")
  }
  if (get_rank_(lhs) != q) {
    stop("'lhs' must have full row rank")
  }
  lhs
}

check_hypothesis <- function(lhs, rhs, p) {
  if (is.null(rhs) && is.null(lhs)) {
    stop("either 'rhs' or 'lhs' must be provided")
  } else if (is.null(lhs)) {
    rhs <- check_rhs_(rhs, p)
  } else if (is.null(rhs)) {
    lhs <- check_lhs_(lhs, p)
    rhs <- rep(0, nrow(lhs))
  } else {
    lhs <- check_lhs_(lhs, p)
    rhs <- check_rhs_(rhs, nrow(lhs))
  }
  list(l = lhs, r = rhs)
}

check_cv <- function(cv, th) {
  cv <- tryCatch(as.numeric(cv),
    warning = function(w) NA,
    error = function(e) NA
  )
  if (any(length(cv) != 1L, is.na(cv), is.infinite(cv))) {
    stop("'cv' is not a number")
  }
  if (cv < .Machine$double.eps) {
    stop("'cv' is too small")
  }
  if (!is.null(th) && cv > 2 * th) {
    stop("'cv' is too large")
  }
  cv
}

check_family <- function(family) {
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



# check_mht_ <- function(lhs, rhs, p) {
#   if (is.null(rhs) && is.null(lhs)) {
#     stop("either 'rhs' or 'lhs' must be provided")
#   } else if (is.null(lhs)) {
#     rhs <- check_rhs2_(rhs, p)
#     lhs <- matrix(rep(diag(1, nrow = p, ncol = p), attr(rhs, "m")),
#       ncol = p,
#       byrow = TRUE
#     )
#   } else if (is.null(rhs)) {
#     lhs <- check_lhs2_(lhs, p)
#     rhs <- rep(0, nrow(lhs))
#   } else {
#     rhs <- check_rhs2_(rhs, p)
#     lhs <- check_lhs2_(lhs, p)
#   }
#   q <- attr(lhs, "q")
#   m <- attr(lhs, "m")
#   if (any(
#     length(rhs) != nrow(lhs), isFALSE(all.equal(attr(rhs, "q"), q)),
#     attr(rhs, "m") != m
#   )) {
#     stop("'rhs' and 'lhs' are incompatible")
#   }
#   list(r = rhs, l = lhs, q = q, m = m)
# }
#
# check_rhs2_ <- function(rhs, p) {
#   UseMethod("check_rhs2_", rhs)
# }
#
# check_rhs2_.numeric <- function(rhs, p) {
#   if (isFALSE(all(is.finite(rhs)))) {
#     stop("'rhs' is not a finite numeric vector")
#   }
#   m <- length(rhs)
#   attr(rhs, "q") <- c(0L, cumsum(rep(1L, m)))
#   attr(rhs, "m") <- m
#   rhs
# }
#
# check_rhs2_.list <- function(rhs, p) {
#   if (isFALSE(all(vapply(rhs, is.vector, TRUE)))) {
#     stop("'rhs' is not a list of vectors")
#   }
#   m <- length(rhs)
#   if (m < 2L) {
#     stop("'rhs' does not specify multiple hypotheses")
#   }
#   if (any(vapply(rhs, \(x) {!is.numeric(x) || !all(is.finite(x))}, FALSE))) {
#     stop("'rhs' is not a list of finite numeric vectors")
#   }
#   if (any(vapply(rhs, \(x) {length(x) != p}, FALSE))) {
#     stop(gettextf("length of every vector in 'rhs' should be %s", sQuote(p),
#       domain = NA
#     ))
#   }
#   out <- do.call(c, rhs)
#   attr(out, "q") <- c(0L, cumsum(vapply(rhs, length, 1L)))
#   attr(out, "m") <- m
#   out
# }
#
# check_lhs2_ <- function(lhs, p) {
#   UseMethod("check_lhs2_", lhs)
# }
#
# check_lhs2_.matrix <- function(lhs, p) {
#   if (!is.numeric(lhs) || !all(is.finite(lhs))) {
#     stop("'lhs' is not a finite numeric matrix")
#   }
#   if (ncol(lhs) != p) {
#     stop("'object' and 'lhs' have incompatible dimensions")
#   }
#   m <- nrow(lhs)
#   if (m < 2L) {
#     stop("'lhs' does not specify multiple hypotheses")
#   }
#   if (isFALSE(all(apply(lhs, 1L, get_rank_)))) {
#     stop("every row of 'lhs' must be a nonzero vector")
#   }
#   attr(lhs, "q") <- c(0L, cumsum(rep(1L, m)))
#   attr(lhs, "m") <- m
#   lhs
# }
#
# check_lhs2_.list <- function(lhs, p) {
#   if (isFALSE(all(vapply(lhs, is.matrix, TRUE)))) {
#     stop("'lhs' is not a list of matrices")
#   }
#   m <- length(lhs)
#   if (m < 2L) {
#     stop("'lhs' does not specify multiple hypotheses")
#   }
#   if (any(vapply(lhs, \(x) {!is.numeric(x) || !all(is.finite(x))}, FALSE))) {
#     stop("'lhs' is not a list of finite numeric matrices")
#   }
#   if (any(vapply(lhs, \(x) {ncol(x) != p}, FALSE))) {
#     stop("'object' and 'lhs' have incompatible dimensions")
#   }
#   if (any(vapply(lhs, \(x) {nrow(x) == 0L || nrow(x) > p}, FALSE))) {
#     stop("'object' and 'lhs' have incompatible dimensions")
#   }
#   if (any(vapply(lhs, \(x) {get_rank_(x) != nrow(x)}, FALSE))) {
#     stop("every matrix in 'lhs' must have full row rank")
#   }
#   out <- do.call(rbind, lhs)
#   attr(out, "q") <- c(0L, cumsum(vapply(lhs, nrow, 1L)))
#   attr(out, "m") <- m
#   out
# }
