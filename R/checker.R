check_weights <- function(weights, nw) {
  if (is.null(weights)) {
    return(numeric(length = 0L))
  }
  if (!is.numeric(weights))
    stop(gettextf("%s is not a numeric vector", sQuote("weights")), domain = NA)
  w <- as.vector(weights, mode = "numeric")
  if (!all(is.finite(w)))
    stop(gettextf("%s is not a finite numeric vector", sQuote("weights")),
         domain = NA)
  if (any(w < 0))
    stop(gettextf("negative %s not allowed", sQuote("weights")), domain = NA)
  if (length(w) != nw)
    stop(gettextf("length of %s is incompatible with data", sQuote("weights")),
         domain = NA)
  w <- (nw / sum(w)) * w
  w
}

check_rhs <- function(rhs, p) {
  rhs <- as.vector(rhs, mode = "numeric")
  if (!is.numeric(rhs) || !all(is.finite(rhs)))
    stop(gettextf("%s is not a finite numeric vector", sQuote("rhs")),
         domain = NA)
  if (length(rhs) != p)
    stop(paste("length of 'rhs' should be "), p)
  rhs
}

check_lhs <- function(lhs, p) {
  lhs <- as.matrix(lhs)
  if (!is.numeric(lhs) || !all(is.finite(lhs)))
    stop(gettextf("%s is not a finite numeric matrix", sQuote("lhs")),
         domain = NA)
  if (ncol(lhs) != p)
    stop(gettextf("%s and %s have incompatible dimensions", sQuote("object"),
                  sQuote("lhs")), domain = NA)
  q <- nrow(lhs)
  if (q == 0L || q > p)
    stop(gettextf("%s and %s have incompatible dimensions", sQuote("object"),
                  sQuote("lhs")), domain = NA)
  if (get_rank_(lhs) != q) {
    stop(gettextf("%s does not have full row rank", sQuote("lhs")), domain = NA)
  }
  lhs
}

check_hypothesis <- function(lhs, rhs, p) {
  if (is.null(rhs) && is.null(lhs)) {
    stop(gettextf("either %s or %s must be provided", sQuote("rhs"),
                  sQuote("lhs")), domain = NA)
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

check_family <- function(family) {
  f <- family$family
  l <- family$link
  switch(f,
    "gaussian" = {
      if (!any(l == c("identity", "log", "inverse")))
        stop(gettextf("%s family with %s link not supported by el_glm",
                      sQuote(f), sQuote(l)), domain = NA)
      },
    "binomial" = {
      if (!any(l == c("logit", "probit")))
        stop(gettextf("%s family with %s link not supported by el_glm",
                      sQuote(f), sQuote(l)), domain = NA)
      },
    stop(gettextf("%s family not supported by el_glm", sQuote(f)),
         domain = NA)
  )
  list(family = f, link = l)
}



