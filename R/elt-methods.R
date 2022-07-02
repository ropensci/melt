#' @importFrom methods is
#' @importFrom stats pchisq
#' @rdname elt
#' @srrstats {G2.3a} `match.arg()` is used to the `calibrate` argument.
setMethod("elt", "EL", function(object,
                                rhs = NULL,
                                lhs = NULL,
                                alpha = 0.05,
                                calibrate = c("chisq", "boot", "f"),
                                control = el_control()) {
  if (is.character(calibrate)) {
    calibrate <- tolower(calibrate)
  }

  stopifnot(
    "'object' has no 'data'; fit the model with 'model' = TRUE" =
      (length(getDataMatrix(object)) > 1L),
    "invalid 'control' specified" = (is(control, "ControlEL"))
  )
  h <- validate_hypothesis(rhs, lhs, object@npar)
  alpha <- validate_alpha(alpha)
  calibrate <- match.arg(calibrate)
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
    out <- computeEL(method, par, getDataMatrix(object), maxit_l, tol_l, th, w)
    p <- length(par)
    cal <- calibrate(calibrate, alpha, out$statistic, p, par, object, control)
    return(new("ELT",
      optim = out$optim, alpha = alpha, logl = out$logl,
      statistic = out$statistic, cv = cal["cv"], pval = cal["pval"],
      calibrate = calibrate
    ))
  }
  # proceed with chi-square calibration for non-NULL 'lhs'
  stopifnot(
    "bootstrap calibration is applicable only when 'lhs' is NULL" =
      (calibrate != "boot"),
    "F calibration is applicable only when 'lhs' is NULL" = (calibrate != "f")
  )
  out <- testHypothesis(
    method, coef(object), getDataMatrix(object), h$l, h$r,
    maxit, maxit_l, tol, tol_l, step, th, w
  )
  new("ELT",
    optim = out$optim, alpha = alpha, logl = out$logl,
    statistic = out$statistic, cv = qchisq(p = 1 - alpha, df = nrow(h$l)),
    pval = pchisq(out$statistic, df = nrow(h$l), lower.tail = FALSE),
    calibrate = calibrate
  )
})

#' @rdname elt
setMethod("elt", "missing", function(object,
                                     rhs = NULL,
                                     lhs = NULL,
                                     alpha,
                                     calibrate = c("chisq", "boot", "f"),
                                     control = el_control()) {
  stopifnot(
    # "'object' has no 'data'; fit the model with 'model' = TRUE" =
    #   (length(getDataMatrix(object)) > 1L),
    "invalid 'control' specified" = (is(control, "ControlEL"))
  )
  # h <- validate_hypothesis(rhs, lhs, object@npar)
  alpha <- validate_alpha(alpha)
  calibrate <- match.arg(calibrate)

  if (length(rhs) != 1L) {
    stop("sefse")
  }
  NULL
})
