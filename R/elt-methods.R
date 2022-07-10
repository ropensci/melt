#' @importFrom methods is
#' @importFrom stats pchisq
#' @rdname elt
setMethod("elt", "EL", function(object,
                                rhs = NULL,
                                lhs = NULL,
                                alpha = 0.05,
                                calibrate = "chisq",
                                control = el_control()) {
  npar <- getNumPar(object)
  stopifnot(
    "`elt()` method is not applicable to to an empty model." = (npar != 0L),
    "`object` has no `data`. Fit the model with `keep_data == TRUE`." =
      (!is.null(getDataMatrix(object))),
    "Invalid `control` specified." = (is(control, "ControlEL"))
  )
  h <- validate_hypothesis(rhs, lhs, npar)
  alpha <- validate_alpha(alpha)
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
    if (calibrate == "f" && method != "mean") {
      stop("F calibration is applicable only to the mean")
    }
    par <- h$r
    out <- compute_EL(method, par, getDataMatrix(object), maxit_l, tol_l, th, w)
    optim <- validate_optim(out$optim)
    names(optim$par) <- names(coef(object))
    p <- length(par)
    cal <- calibrate(calibrate, alpha, out$statistic, p, par, object, control)
    return(new("ELT",
      optim = optim, alpha = alpha, logl = out$logl, loglr = out$loglr,
      statistic = out$statistic, cv = cal["cv"], pval = cal["pval"],
      calibrate = calibrate
    ))
  }
  # Proceed with chi-square calibration for non-NULL `lhs`
  stopifnot(
    "Bootstrap calibration is applicable only when `lhs` is `NULL`." =
      (calibrate != "boot"),
    "F calibration is applicable only when `lhs` is `NULL`." = (calibrate != "f")
  )
  out <- test_hypothesis(
    method, coef(object), getDataMatrix(object), h$l, h$r,
    maxit, maxit_l, tol, tol_l, step, th, w
  )
  optim <- validate_optim(out$optim)
  names(optim$par) <- names(coef(object))
  new("ELT",
    optim = optim, alpha = alpha, logl = out$logl, loglr = out$loglr,
    statistic = out$statistic, cv = qchisq(p = 1 - alpha, df = nrow(h$l)),
    pval = pchisq(out$statistic, df = nrow(h$l), lower.tail = FALSE),
    calibrate = calibrate
  )
})

#' @importFrom methods is
#' @importFrom stats pchisq
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
    "`elt()` method is not applicable to to an empty model." = (npar != 0L),
    "`object` has no `data`. Fit the model with `keep_data == TRUE`." =
      (!is.null(getDataMatrix(object))),
    "Invalid `control` specified." = (is(control, "ControlEL"))
  )
  h <- validate_hypothesis(rhs, lhs, npar)
  alpha <- validate_alpha(alpha)
  calibrate <- validate_calibrate(calibrate)
  nm <- names(coef(object))
  maxit_l <- control@maxit_l
  tol_l <- control@tol_l
  th <- control@th
  w <- getWeights(object)
  if (is.null(lhs)) {
    if (isTRUE(calibrate == "f")) {
      stop("F calibration is applicable only to the mean")
    }
    par <- h$r
    stopifnot(
      "Parametr value must be a finite single numeric." =
        (isTRUE(is.numeric(par) && length(par) == 1L && is.finite(par))),
      "Parametr value must be a positive single numeric." =
        (par >= .Machine$double.eps)
    )
    out <- compute_EL("sd", par, getDataMatrix(object), maxit_l, tol_l, th, w)
    optim <- validate_optim(out$optim)
    names(optim$par) <- nm
    cal <- calibrate(calibrate, alpha, out$statistic, 1L, par, object, control)
    return(new("ELT",
      optim = optim, alpha = alpha, logl = out$logl, loglr = out$loglr,
      statistic = out$statistic, cv = cal["cv"], pval = cal["pval"],
      calibrate = calibrate
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
  out <- compute_EL("sd", par, getDataMatrix(object), maxit_l, tol_l, th, w)
  optim <- validate_optim(out$optim)
  names(optim$par) <- nm
  new("ELT",
    optim = optim, alpha = alpha, logl = out$logl, loglr = out$loglr,
    statistic = out$statistic, cv = qchisq(p = 1 - alpha, df = 1L),
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
  stopifnot(
    "Invalid `control` specified." = (is(control, "ControlEL"))
  )
  alpha <- validate_alpha(alpha)
  NULL
})
