#' @importFrom methods is
#' @importFrom stats pchisq
#' @rdname elt
setMethod("elt", "EL", function(object,
                                rhs = NULL,
                                lhs = NULL,
                                alpha = 0.05,
                                calibrate = c("chisq", "boot", "f"),
                                control = el_control()) {
  stopifnot(
    "invalid 'object' supplied" = (is(object, "EL")),
    "'object' has no 'data'; fit the model with 'model' = TRUE" =
      (length(getDataMatrix(object)) > 1L),
    "invalid 'control' specified" = (is(control, "ControlEL"))
  )
  h <- validate_hypothesis(lhs, rhs, object@npar)
  validate_alpha(alpha)
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
    el <- eval_(method, par, getDataMatrix(object), maxit_l, tol_l, th, w)
    p <- length(par)
    cal <- calibrate(alpha, el$statistic, calibrate, p, par, object, control)
    return(new("ELT",
               optim = el$optim, alpha = alpha, logl = el$logl, statistic = el$statistic,
               cv = cal["cv"], pval = cal["pval"], calibrate = calibrate
    ))
  }
  # proceed with chi-square calibration for non-NULL 'lhs'
  stopifnot(
    "bootstrap calibration is applicable only when 'lhs' is NULL" =
      (calibrate != "boot"),
    "F calibration is applicable only when 'lhs' is NULL" = (calibrate != "f")
  )
  el <- elt_(
    method, coef(object), getDataMatrix(object), h$l, h$r,
    maxit, maxit_l, tol, tol_l, step, th, w
  )
  new("ELT",
      optim = el$optim, alpha = alpha, logl = el$logl, statistic = el$statistic,
      cv = qchisq(p = 1 - alpha, df = nrow(h$l)),
      pval = pchisq(el$statistic, df = nrow(h$l), lower.tail = FALSE),
      calibrate = calibrate
  )
})
