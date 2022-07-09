#' @importFrom methods is
#' @rdname elmt
setMethod("elmt", "EL", function(object,
                                 rhs = NULL,
                                 lhs = NULL,
                                 alpha = 0.05,
                                 control = el_control()) {
  npar <- getNumPar(object)
  stopifnot(
    "`elmt()` method is not applicable to to an empty model." = (npar != 0L),
    "`elmt()` method is not applicable to a model with one parameter." =
      (npar != 1L),
    "`object` has no `data`. Fit the model with `keep_data == TRUE`." =
      (!is.null(getDataMatrix(object))),
    "Invalid `control` specified." = (is(control, "ControlEL"))
  )
  h <- validate_hypotheses(rhs, lhs, npar)
  alpha <- validate_alpha(alpha)
  method <- getMethodEL(object)
  maxit <- control@maxit
  maxit_l <- control@maxit_l
  tol <- control@tol
  tol_l <- control@tol_l
  step <- control@step
  th <- control@th
  m <- control@m
  out <- test_multiple_hypotheses(
    alpha, h$q, h$m, m, method, coef(object), getDataMatrix(object), h$r,
    h$l, maxit, maxit_l, tol, tol_l, step, th, getWeights(object)
  )
  new("ELMT",
    alpha = alpha, statistic = out$statistic, cv = out$cv, pval = out$pval,
    calibrate = "mvchisq"
  )
})
