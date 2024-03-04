#' @rdname elt
setMethod("elt", "EL", function(object,
                                rhs = NULL,
                                lhs = NULL,
                                alpha = 0.05,
                                calibrate = "chisq",
                                control = NULL) {
  stopifnot(
    "`elt()` is not applicable to an empty model." = getDF(object) >= 1L,
    "`object` has no `data`. Fit the model with `keep_data == TRUE`." =
      isFALSE(is.null(getData(object)))
  )
  if (is.null(control)) {
    control <- getControlEL(object)
  } else {
    assert_class(control, "ControlEL")
  }
  onames <- names(logProb(object))
  pnames <- names(getOptim(object)$par)
  h <- validate_hypothesis(rhs, lhs, getNumPar(object), pnames)
  alpha <- assert_number(alpha, lower = 0, upper = 1, finite = TRUE)
  calibrate <- validate_calibrate(calibrate)
  method <- getMethodEL(object)
  maxit <- control@maxit
  maxit_l <- control@maxit_l
  tol <- control@tol
  tol_l <- control@tol_l
  step <- control@step
  th <- control@th
  w <- getWeights(object)
  if (is.null(lhs)) {
    stopifnot(
      "F calibration is applicable only to the mean." =
        (isTRUE(calibrate != "f" || method == "mean"))
    )
    par <- h$r
    out <- compute_EL(method, par, getData(object), maxit_l, tol_l, th, w)
    optim <- validate_optim(out$optim)
    names(optim$par) <- pnames
    optim$cstr <- FALSE
    cal <- calibrate(
      calibrate, alpha, out$statistic, length(par), par, object, control
    )
    return(new("ELT",
      optim = optim, logp = setNames(out$logp, onames), logl = out$logl,
      loglr = out$loglr, statistic = out$statistic, df = length(par),
      pval = unname(cal["pval"]), cv = unname(cal["cv"]), rhs = par, lhs = h$l,
      alpha = alpha, calibrate = calibrate, control = control
    ))
  }
  # Proceed with chi-square calibration for non-NULL `lhs`
  stopifnot(
    "Bootstrap calibration is applicable only when `lhs` is NULL." =
      (calibrate != "boot"),
    "F calibration is applicable only when `lhs` is NULL." =
      (calibrate != "f")
  )
  out <- test_hypothesis(
    method, coef(object), getData(object), h$l, h$r, maxit, maxit_l, tol, tol_l,
    step, th, w
  )
  optim <- validate_optim(out$optim)
  names(optim$par) <- pnames
  optim$cstr <- TRUE
  q <- nrow(h$l)
  new("ELT",
    optim = optim, logp = setNames(out$logp, onames), logl = out$logl,
    loglr = out$loglr, statistic = out$statistic, df = q,
    pval = pchisq(out$statistic, df = q, lower.tail = FALSE),
    cv = qchisq(1 - alpha, df = q), rhs = h$r, lhs = h$l, alpha = alpha,
    calibrate = calibrate, control = control
  )
})

#' @rdname elt
#' @usage NULL
setMethod("elt", "QGLM", function(object,
                                  rhs = NULL,
                                  lhs = NULL,
                                  alpha = 0.05,
                                  calibrate = "chisq",
                                  control = NULL) {
  stopifnot(
    "`elt()` is not applicable to an empty model." = getDF(object) >= 1L,
    "`object` has no `data`. Fit the model with `keep_data == TRUE`." =
      isFALSE(is.null(getData(object)))
  )
  if (is.null(control)) {
    control <- getControlEL(object)
  } else {
    assert_class(control, "ControlEL")
  }
  # `npar` includes the dispersion parameter
  onames <- names(logProb(object))
  nm <- names(getOptim(object)$par)
  pnames <- nm[-getNumPar(object)]
  h <- validate_hypothesis(rhs, lhs, getNumPar(object) - 1L, pnames)
  alpha <- assert_number(alpha, lower = 0, upper = 1, finite = TRUE)
  calibrate <- validate_calibrate(calibrate)
  method <- getMethodEL(object)
  maxit <- control@maxit
  maxit_l <- control@maxit_l
  tol <- control@tol
  tol_l <- control@tol_l
  step <- control@step
  th <- control@th
  w <- getWeights(object)
  if (is.null(lhs)) {
    stopifnot(
      "F calibration is applicable only to the mean." = (calibrate != "f")
    )
    if (calibrate == "boot") {
      par <- c(h$r, compute_dispersion(object, h$r))
    } else {
      par <- c(h$r, object@dispersion)
    }
    out <- compute_EL(method, par, getData(object), maxit_l, tol_l, th, w)
    optim <- validate_optim(out$optim)
    names(optim$par) <- nm
    optim$cstr <- TRUE
    cal <- calibrate(
      calibrate, alpha, out$statistic, length(par), par, object, control
    )
    return(new("ELT",
      optim = optim, logp = setNames(out$logp, onames), logl = out$logl,
      loglr = out$loglr, statistic = out$statistic, df = length(par),
      pval = unname(cal["pval"]), cv = unname(cal["cv"]), rhs = h$r, lhs = h$l,
      alpha = alpha, calibrate = calibrate, control = control
    ))
  }
  stopifnot(
    "Bootstrap calibration is applicable only when `lhs` is NULL." =
      (calibrate != "boot"),
    "F calibration is applicable only to the mean." = (calibrate != "f")
  )
  l <- cbind(h$l, 0)
  out <- test_hypothesis(
    method, c(coef(object), object@dispersion), getData(object), l, h$r, maxit,
    maxit_l, tol, tol_l, step, th, w
  )
  optim <- validate_optim(out$optim)
  names(optim$par) <- nm
  optim$cstr <- TRUE
  q <- nrow(l)
  new("ELT",
    optim = optim, logp = setNames(out$logp, onames), logl = out$logl,
    loglr = out$loglr, statistic = out$statistic, df = q,
    pval = pchisq(out$statistic, df = q, lower.tail = FALSE),
    cv = qchisq(1 - alpha, df = q), rhs = h$r, lhs = h$l, alpha = alpha,
    calibrate = calibrate, control = control
  )
})

#' @rdname elt
#' @usage NULL
setMethod("elt", "SD", function(object,
                                rhs = NULL,
                                lhs = NULL,
                                alpha = 0.05,
                                calibrate = "chisq",
                                control = NULL) {
  stopifnot(
    "`elt()` is not applicable to an empty model." = getDF(object) >= 1L,
    "`object` has no `data`. Fit the model with `keep_data == TRUE`." =
      isFALSE(is.null(getData(object)))
  )
  if (is.null(control)) {
    control <- getControlEL(object)
  } else {
    assert_class(control, "ControlEL")
  }
  onames <- names(logProb(object))
  pnames <- names(getOptim(object)$par)
  h <- validate_hypothesis(rhs, lhs, getNumPar(object), pnames)
  alpha <- assert_number(alpha, lower = 0, upper = 1, finite = TRUE)
  calibrate <- validate_calibrate(calibrate)
  maxit_l <- control@maxit_l
  tol_l <- control@tol_l
  th <- control@th
  w <- getWeights(object)
  if (is.null(lhs)) {
    par <- h$r
    stopifnot(
      "F calibration is applicable only to the mean." = (calibrate != "f"),
      "Parametr value must be a finite single numeric." =
        (isTRUE(is.numeric(par) && length(par) == 1L && is.finite(par))),
      "Parametr value must be a positive single numeric." =
        (par >= .Machine$double.eps)
    )
    out <- compute_EL("sd", par, getData(object), maxit_l, tol_l, th, w)
    optim <- validate_optim(out$optim)
    names(optim$par) <- pnames
    optim$cstr <- FALSE
    cal <- calibrate(calibrate, alpha, out$statistic, 1L, par, object, control)
    return(new("ELT",
      optim = optim, logp = setNames(out$logp, onames), logl = out$logl,
      loglr = out$loglr, statistic = out$statistic, df = 1L,
      pval = unname(cal["pval"]), cv = unname(cal["cv"]), rhs = h$r, lhs = h$l,
      alpha = alpha, calibrate = calibrate, control = control
    ))
  }
  stopifnot(
    "Bootstrap calibration is applicable only when `lhs` is `NULL`." =
      (calibrate != "boot"),
    "F calibration is applicable only when `lhs` is `NULL`." =
      (calibrate != "f")
  )
  par <- solve(h$l, h$r)
  stopifnot(
    "Parametr value must be a finite single numeric." =
      (isTRUE(is.numeric(par) && length(par) == 1L && is.finite(par))),
    "Parametr value must be a positive single numeric." =
      (par >= .Machine$double.eps)
  )
  out <- compute_EL("sd", par, getData(object), maxit_l, tol_l, th, w)
  optim <- validate_optim(out$optim)
  names(optim$par) <- pnames
  optim$cstr <- FALSE
  new("ELT",
    optim = optim, logp = setNames(out$logp, onames), logl = out$logl,
    loglr = out$loglr, statistic = out$statistic, df = 1L,
    pval = pchisq(out$statistic, df = 1L, lower.tail = FALSE),
    cv = qchisq(1 - alpha, df = 1L), rhs = h$r, lhs = h$l, alpha = alpha,
    calibrate = calibrate, control = control
  )
})

#' @rdname elt
#' @usage NULL
setMethod("elt", "missing", function(object,
                                     rhs = NULL,
                                     lhs = NULL,
                                     alpha = 0.05,
                                     calibrate = "chisq",
                                     control = NULL) {
  alpha <- assert_number(alpha, lower = 0, upper = 1, finite = TRUE)
  calibrate <- validate_calibrate(calibrate)
  if (is.null(control)) {
    control <- el_control()
  } else {
    assert_class(control, "ControlEL")
  }
  NULL
})
