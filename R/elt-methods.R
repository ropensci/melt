#' @rdname elt
setMethod("elt", "EL", function(object,
                                rhs = NULL,
                                lhs = NULL,
                                alpha = 0.05,
                                calibrate = "chisq",
                                control = el_control()) {
  npar <- getNumPar(object)
  stopifnot(
    "`elt()` method is not applicable to an empty model." = (npar != 0L),
    "`object` has no `data`. Fit the model with `keep_data == TRUE`." =
      (isFALSE(is.null(getData(object)))),
    "Invalid `control` specified." = (is(control, "ControlEL"))
  )
  h <- validate_hypothesis(rhs, lhs, npar)
  alpha <- validate_alpha(alpha)
  calibrate <- validate_calibrate(calibrate)
  nm <- names(getOptim(object)$par)
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
    names(optim$par) <- nm
    p <- length(par)
    cal <- calibrate(calibrate, alpha, out$statistic, p, par, object, control)
    return(new("ELT",
      optim = optim, alpha = alpha, logl = out$logl, loglr = out$loglr,
      statistic = out$statistic, cv = unname(cal["cv"]),
      pval = unname(cal["pval"]), calibrate = calibrate
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
  names(optim$par) <- nm
  new("ELT",
    optim = optim, alpha = alpha, logl = out$logl, loglr = out$loglr,
    statistic = out$statistic, cv = qchisq(1 - alpha, df = nrow(h$l)),
    pval = pchisq(out$statistic, df = nrow(h$l), lower.tail = FALSE),
    calibrate = calibrate
  )
})

#' @rdname elt
#' @usage NULL
setMethod("elt", "QGLM", function(object,
                                  rhs = NULL,
                                  lhs = NULL,
                                  alpha = 0.05,
                                  calibrate = "chisq",
                                  control = el_control()) {
  # `npar` includes the dispersion parameter
  p <- getNumPar(object) - 1L
  stopifnot(
    "`elt()` method is not applicable to an empty model." = (p != 0L),
    "`object` has no `data`. Fit the model with `keep_data == TRUE`." =
      (isFALSE(is.null(getData(object)))),
    "Invalid `control` specified." = (is(control, "ControlEL"))
  )
  h <- validate_hypothesis(rhs, lhs, p)
  alpha <- validate_alpha(alpha)
  calibrate <- validate_calibrate(calibrate)
  nm <- names(getOptim(object)$par)
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
    par <- c(h$r, object@dispersion)
    out <- compute_EL(method, par, getData(object), maxit_l, tol_l, th, w)
    optim <- validate_optim(out$optim)
    names(optim$par) <- nm
    npar <- length(par)
    cal <- calibrate(
      calibrate, alpha, out$statistic, npar, par, object, control
    )
    return(new("ELT",
      optim = optim, alpha = alpha, logl = out$logl, loglr = out$loglr,
      statistic = out$statistic, cv = unname(cal["cv"]),
      pval = unname(cal["pval"]), calibrate = calibrate
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
  new("ELT",
    optim = optim, alpha = alpha, logl = out$logl, loglr = out$loglr,
    statistic = out$statistic, cv = qchisq(1 - alpha, df = nrow(l)),
    pval = pchisq(out$statistic, df = nrow(l), lower.tail = FALSE),
    calibrate = calibrate
  )
})

#' @rdname elt
#' @usage NULL
setMethod("elt", "SD", function(object,
                                rhs = NULL,
                                lhs = NULL,
                                alpha = 0.05,
                                calibrate = "chisq",
                                control = el_control()) {
  npar <- getNumPar(object)
  stopifnot(
    "`elt()` method is not applicable to an empty model." = (npar != 0L),
    "`object` has no `data`. Fit the model with `keep_data == TRUE`." =
      (isFALSE(is.null(getData(object)))),
    "Invalid `control` specified." = (is(control, "ControlEL"))
  )
  h <- validate_hypothesis(rhs, lhs, npar)
  alpha <- validate_alpha(alpha)
  calibrate <- validate_calibrate(calibrate)
  nm <- names(getOptim(object)$par)
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
    names(optim$par) <- nm
    cal <- calibrate(calibrate, alpha, out$statistic, 1L, par, object, control)
    return(new("ELT",
      optim = optim, alpha = alpha, logl = out$logl, loglr = out$loglr,
      statistic = out$statistic, cv = unname(cal["cv"]),
      pval = unname(cal["pval"]), calibrate = calibrate
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
  names(optim$par) <- nm
  new("ELT",
    optim = optim, alpha = alpha, logl = out$logl, loglr = out$loglr,
    statistic = out$statistic, cv = qchisq(1 - alpha, df = 1L),
    pval = pchisq(out$statistic, df = 1L, lower.tail = FALSE),
    calibrate = calibrate
  )
})

#' @rdname elt
#' @usage NULL
setMethod("elt", "missing", function(object,
                                     rhs = NULL,
                                     lhs = NULL,
                                     alpha = 0.05,
                                     calibrate = "chisq",
                                     control = el_control()) {
  alpha <- validate_alpha(alpha)
  calibrate <- validate_calibrate(calibrate)
  stopifnot("Invalid `control` specified." = (is(control, "ControlEL")))
  alpha <- validate_alpha(alpha)
  NULL
})
