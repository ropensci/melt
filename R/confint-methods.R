#' @rdname confint
setMethod(
  "confint", "EL",
  function(object, parm, level = 0.95, ..., cv = qchisq(level, 1L),
           control = el_control()) {
    est <- object@coefficients
    # no confidence interval for empty model
    if (length(est) == 0L) {
      ci <- matrix(, nrow = 0L, ncol = 2L)
      colnames(ci) <- c("lower", "upper")
      return(ci)
    }
    # index for tracking the parameters
    idx <- seq(length(est))
    # rownames of the confidence interval matrix
    pnames <- if (is.null(names(est))) idx else names(est)
    # if parm is supplied, modify idx and pnames accordingly
    if (!missing(parm)) {
      if (is.numeric(parm) && all(is.finite(parm))) {
        pnames <- pnames[parm]
        # idx <- match(pnames, names(est))
        idx <- if (is.null(names(est))) {
          match(pnames, idx)
        } else {
          match(pnames, names(est))
        }
      } else if (is.character(parm)) {
        idx <- match(parm, pnames)
        pnames <- parm
      } else {
        stop("invalid 'parm' specified")
      }
    }
    # number of rows of the confidence interval matrix
    p <- length(idx)
    # check level
    if (!missing(level) &&
        (length(level) != 1L || !is.finite(level) || level < 0 || level > 1)) {
      stop("'level' must be a number between 0 and 1")
    }
    if (isTRUE(all.equal(level, 0))) {
      ci <- matrix(rep(est[idx], 2L), ncol = 2L)
      colnames(ci) <- c("lower", "upper")
      return(ci)
    } else if (isTRUE(all.equal(level, 1))) {
      ci <- matrix(NA, nrow = p, ncol = 2L)
      ci[which(!is.na(idx)), ] <- c(-Inf, Inf)
      colnames(ci) <- c("lower", "upper")
      return(ci)
    }

    if (!is(control, "ControlEL")) {
      stop("invalid 'control' specified")
    }
    method <- getMethodEL(object)
    maxit <- control@maxit
    maxit_l <- control@maxit_l
    tol <- control@tol
    tol_l <- control@tol_l
    step <- control@step
    th <- control@th
    nthreads <- control@nthreads
    w <- if (is.null(object@weights)) numeric(length = 0L) else object@weights
    cv <- check_cv(cv, th)
    # compute the confidence interval matrix
    if (isTRUE(all.equal(level, 0))) {
      ci <- matrix(rep(est[idx], 2L), ncol = 2L)
    } else if (isTRUE(all.equal(level, 1))) {
      ci <- matrix(NA, nrow = p, ncol = 2L)
      ci[which(!is.na(idx)), ] <- c(-Inf, Inf)
    } else if (all(is.na(idx))) {
      ci <- matrix(NA, nrow = p, ncol = 2L)
    } else if (any(is.na(idx))) {
      idx_na <- which(is.na(idx))
      ci <- matrix(NA, nrow = p, ncol = 2L)
      ci[-idx_na, ] <- confint_(method, est, object@data, cv, idx[-idx_na],
                                maxit, maxit_l, tol, tol_l, step, th, nthreads,
                                w)
    } else {
      ci <- confint_(method, est, object@data, cv, idx, maxit, maxit_l,
                     tol, tol_l, step, th, nthreads, w)
    }
    dimnames(ci) <- list(pnames, c("lower", "upper"))
    ci
  }
)




