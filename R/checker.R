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
  if (any(w < 0)) {
    stop("negative 'weights' not allowed")
  }
  if (length(w) != nw) {
    stop("length of 'weights' is incompatible with data")
  }
  w <- (nw / sum(w)) * w
  w
}

check_rhs <- function(rhs, p) {
  rhs <- as.vector(rhs, mode = "numeric")
  if (!is.numeric(rhs) || !all(is.finite(rhs))) {
    stop("'rhs' is not a finite numeric vector")
  }
  if (length(rhs) != p) {
    stop(gettextf("length of 'rhs' should be %d", p, domain = NA))
  }
  rhs
}

check_lhs <- function(lhs, p) {
  lhs <- as.matrix(lhs)
  if (!is.numeric(lhs) || !all(is.finite(lhs))) {
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
    stop("'lhs' does not have full row rank")
  }
  lhs
}

check_hypothesis <- function(lhs, rhs, p) {
  if (is.null(rhs) && is.null(lhs)) {
    stop("either 'rhs' or 'lhs' must be provided")
  } else if (is.null(lhs)) {
    rhs <- check_rhs(rhs, p)
  } else if (is.null(rhs)) {
    lhs <- check_lhs(lhs, p)
    rhs <- rep(0, nrow(lhs))
  } else {
    lhs <- check_lhs(lhs, p)
    rhs <- check_rhs(rhs, nrow(lhs))
  }
  list(l = lhs, r = rhs)
}

check_cv <- function(cv, th) {
  cv <- tryCatch(as.numeric(cv), warning = function(w) NA,
                 error = function(e) NA)
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
        stop(gettextf("%s family with %s link not supported by 'el_glm'",
                      sQuote(f), sQuote(l)), domain = NA)
      }
    },
    "binomial" = {
      if (!any(l == c("logit", "probit", "log"))) {
        stop(gettextf("%s family with %s link not supported by 'el_glm'",
                      sQuote(f), sQuote(l)), domain = NA)
      }
    },
    "poisson" = {
      if (!any(l == c("log", "identity", "sqrt"))) {
        stop(gettextf("%s family with %s link not supported by 'el_glm'",
                      sQuote(f), sQuote(l)), domain = NA)
      }
    },
    stop(gettextf("%s family not supported by 'el_glm'", sQuote(f)),
         domain = NA)
  )
  list(family = f, link = l)
}
