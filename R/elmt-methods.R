#' @importFrom methods is
#' @importFrom stats pchisq
#' @rdname elmt
setMethod("elmt", "EL", function(object,
                                 rhs = NULL,
                                 lhs = NULL,
                                 alpha = 0.05,
                                 control = el_control()) {
  stopifnot(
    "'object' has no 'data'; fit the model with 'model' = TRUE" =
      (length(getDataMatrix(object)) > 1L),
    "invalid 'control' specified" = (is(control, "ControlEL"))
  )
  h <- validate_hypotheses(rhs, lhs, object@npar)
  # return(h)

  alpha <- validate_alpha(alpha)
  method <- getMethodEL(object)
  maxit <- control@maxit
  maxit_l <- control@maxit_l
  tol <- control@tol
  tol_l <- control@tol_l
  step <- control@step
  th <- control@th


  out <- testMultipleHypotheses(
    h$q, h$m, 10000, method, coef(object), getDataMatrix(object), h$r, h$l,
    maxit, maxit_l, tol, tol_l, step, th, getWeights(object)
  )
  new("MELT",
      statistic = out$tmp,
      cutoff = out$tmp2
  )
})
