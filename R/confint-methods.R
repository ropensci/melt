#' @rdname confint
setMethod("confint", "EL", function(object,
                                    parm,
                                    level = 0.95,
                                    cv = NULL,
                                    control = NULL) {
  # No confidence interval for an empty model
  if (getDF(object) == 0L) {
    ci <- matrix(numeric(0), nrow = 0L, ncol = 2L)
    colnames(ci) <- c("lower", "upper")
    return(ci)
  }
  stopifnot(
    "`object` has no `data`. Fit the model with `keep_data == TRUE`." =
      (isFALSE(is.null(getData(object))))
  )
  if (is.null(control)) {
    control <- getControlEL(object)
  } else {
    stopifnot("Invalid `control` specified." = is(control, "ControlEL"))
  }
  est <- coef(object)
  # Index for tracking the parameters
  idx <- seq(length(est))
  # Row names of the confidence interval matrix
  nm <- names(est)
  pnames <- if (is.null(nm)) idx else nm
  # If `parm` is supplied, modify idx and pnames accordingly
  if (!missing(parm)) {
    if (is.numeric(parm) && all(is.finite(parm))) {
      pnames <- pnames[parm]
      idx <- if (is.null(nm)) match(pnames, idx) else match(pnames, nm)
    } else if (is.character(parm)) {
      idx <- match(parm, pnames)
      pnames <- parm
    } else {
      stop("Invalid `parm` specified.")
    }
  }
  # Number of rows of the confidence interval matrix
  p <- length(idx)
  level <- assert_number(level, lower = 0, upper = 1, finite = TRUE)
  if (isTRUE(all.equal(level, 0))) {
    ci <- matrix(rep(est[idx], 2L), ncol = 2L)
    colnames(ci) <- c("lower", "upper")
    return(ci)
  } else if (isTRUE(all.equal(level, 1))) {
    ci <- matrix(NA_real_, nrow = p, ncol = 2L)
    ci[which(!is.na(idx)), ] <- c(if (is(object, "SD")) 0 else -Inf, Inf)
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
      method, getEstimates(object), getData(object), cv, idx[-idx_na], maxit,
      maxit_l, tol, tol_l, step, th, nthreads, w
    )
  } else {
    ci <- compute_confidence_intervals(
      method, getEstimates(object), getData(object), cv, idx, maxit, maxit_l,
      tol, tol_l, step, th, nthreads, w
    )
  }
  dimnames(ci) <- list(pnames, c("lower", "upper"))
  ci
})

#' @rdname confint
setMethod("confint", "ELMT", function(object,
                                      cv = NULL,
                                      control = NULL) {
  stopifnot(
    "Each hypothesis must correspond to a linear combination of parameters." =
      (isTRUE(all(getDF(object) == 1L))),
    "`object` has no `data`. Fit the model with `keep_data == TRUE`." =
      (isFALSE(is.null(getData(object))))
  )
  if (is.null(control)) {
    control <- getControlEL(object)
  } else {
    stopifnot("Invalid `control` specified." = is(control, "ControlEL"))
  }
  method <- getMethodEL(object)
  maxit <- control@maxit
  maxit_l <- control@maxit_l
  tol <- control@tol
  tol_l <- control@tol_l
  step <- control@step
  th <- control@th
  cv <- if (is.null(cv)) object@cv else validate_cv(cv, th)
  nthreads <- control@nthreads
  w <- getWeights(object)
  estimates <- unlist(getEstimates(object))
  ci <- compute_confidence_intervals_EMLT(
    method, getData(object), coef(object), object@lhs, estimates, cv, maxit,
    maxit_l, tol, tol_l, step, th, nthreads, w
  )
  dimnames(ci) <-
    list(
      describe_hypothesis(object@rhs, object@lhs, colnames(object@lhs)),
      c("lower", "upper")
    )
  ci
})
