#' @importFrom stats qchisq
#' @rdname confint
#' @srrstats {RE4.3} `confint()` method is exported.
setMethod("confint", "EL", function(object,
                                    parm,
                                    level = 0.95,
                                    cv = NULL,
                                    control = el_control()) {
  stopifnot("Invalid `control` specified." = (is(control, "ControlEL")))
  # No confidence interval for an empty model
  if (getNumPar(object) == 0L) {
    ci <- matrix(, nrow = 0L, ncol = 2L)
    colnames(ci) <- c("lower", "upper")
    return(ci)
  }
  est <- coef(object)
  # Index for tracking the parameters
  idx <- seq(length(est))
  # Row names of the confidence interval matrix
  pnames <- if (is.null(names(est))) idx else names(est)
  # If `parm` is supplied, modify idx and pnames accordingly
  if (!missing(parm)) {
    if (is.numeric(parm) && all(is.finite(parm))) {
      pnames <- pnames[parm]
      if (is.null(names(est))) {
        idx <- match(pnames, idx)
      } else {
        idx <- match(pnames, names(est))
      }
    } else if (is.character(parm)) {
      idx <- match(parm, pnames)
      pnames <- parm
    } else {
      stop("Invalid `parm` specified.")
    }
  }
  # Number of rows of the confidence interval matrix
  p <- length(idx)
  level <- validate_level(level)
  if (isTRUE(all.equal(level, 0))) {
    ci <- matrix(rep(est[idx], 2L), ncol = 2L)
    colnames(ci) <- c("lower", "upper")
    return(ci)
  } else if (isTRUE(all.equal(level, 1))) {
    ci <- matrix(NA_real_, nrow = p, ncol = 2L)
    ci[which(!is.na(idx)), ] <- c(-Inf, Inf)
    colnames(ci) <- c("lower", "upper")
    return(ci)
  }

  method <- getMethodEL(object)
  maxit <- control@maxit
  maxit_l <- control@maxit_l
  tol <- control@tol
  tol_l <- control@tol_l
  step <- control@step
  th <- control@th
  nthreads <- control@nthreads
  w <- getWeights(object)
  cv <- if (is.null(cv)) qchisq(level, 1L) else validate_cv(cv, th)
  if (all(is.na(idx))) {
    ci <- matrix(NA_real_, nrow = p, ncol = 2L)
  } else if (any(is.na(idx))) {
    idx_na <- which(is.na(idx))
    ci <- matrix(NA_real_, nrow = p, ncol = 2L)
    ci[-idx_na, ] <- compute_confidence_intervals(
      method, est, getDataMatrix(object), cv, idx[-idx_na], maxit, maxit_l,
      tol, tol_l, step, th, nthreads, w
    )
  } else {
    ci <- compute_confidence_intervals(
      method, est, getDataMatrix(object), cv, idx, maxit,
      maxit_l, tol, tol_l, step, th, nthreads, w
    )
  }
  dimnames(ci) <- list(pnames, c("lower", "upper"))
  ci
})

#' @importFrom stats qchisq
#' @rdname confint
setMethod("confint", "SD", function(object,
                                    parm,
                                    level = 0.95,
                                    cv = NULL,
                                    control = el_control()) {
  stopifnot("Invalid `control` specified." = (is(control, "ControlEL")))
  est <- coef(object)
  idx <- seq(length(est))
  pnames <- if (is.null(names(est))) idx else names(est)
  if (!missing(parm)) {
    if (is.numeric(parm) && all(is.finite(parm))) {
      pnames <- pnames[parm]
      if (is.null(names(est))) {
        idx <- match(pnames, idx)
      } else {
        idx <- match(pnames, names(est))
      }
    } else if (is.character(parm)) {
      idx <- match(parm, pnames)
      pnames <- parm
    } else {
      stop("Invalid `parm` specified.")
    }
  }
  p <- length(idx)
  level <- validate_level(level)
  if (isTRUE(all.equal(level, 0))) {
    ci <- matrix(rep(est[idx], 2L), ncol = 2L)
    colnames(ci) <- c("lower", "upper")
    return(ci)
  } else if (isTRUE(all.equal(level, 1))) {
    ci <- matrix(NA_real_, nrow = p, ncol = 2L)
    ci[which(!is.na(idx)), ] <- c(0, Inf)
    colnames(ci) <- c("lower", "upper")
    return(ci)
  }

  maxit <- control@maxit
  maxit_l <- control@maxit_l
  tol <- control@tol
  tol_l <- control@tol_l
  step <- control@step
  th <- control@th
  nthreads <- control@nthreads
  w <- getWeights(object)
  cv <- if (is.null(cv)) qchisq(level, 1L) else validate_cv(cv, th)
  if (all(is.na(idx))) {
    ci <- matrix(NA_real_, nrow = p, ncol = 2L)
  } else if (any(is.na(idx))) {
    idx_na <- which(is.na(idx))
    ci <- matrix(NA_real_, nrow = p, ncol = 2L)
    ci[-idx_na, ] <- compute_confidence_intervals(
      "sd", est, getDataMatrix(object), cv, idx[-idx_na], maxit, maxit_l,
      tol, tol_l, step, th, nthreads, w
    )
  } else {
    ci <- compute_confidence_intervals(
      "sd", est, getDataMatrix(object), cv, idx, maxit,
      maxit_l, tol, tol_l, step, th, nthreads, w
    )
  }
  dimnames(ci) <- list(pnames, c("lower", "upper"))
  ci
})
