# el_var <- function(par, x, mean = NULL, weights = NULL, control = el_control(),
#                    model = TRUE) {
#   if (isTRUE(!is.numeric(x) || !all(is.finite(x)) || length(x) != 1L)) {
#     stop("'par' must be a finite number")
#   }
#   if (isFALSE(x >= .Machine$double.eps)) {
#     stop("'par' must be a positive number")
#   }
#   if (isTRUE(!is.numeric(x) || !all(is.finite(x)))) {
#     stop("'x' must be a finite numeric vector")
#   }
#   mm <- as.matrix(x)
#   n <- nrow(mm)
#   p <- ncol(mm)
#   # observations in x
#   if (n < 2L) {
#     stop("not enough 'x' observations")
#   }
#   if (p != 1L) {
#     stop("'x' must be a vector")
#   }
#
#   w <- check_weights(weights, n)
#   if (!is.null(weights)) {
#     #
#     est <- colSums(mm * w) / n
#   } else {
#     #
#     est <- colMeans(mm)
#   }
#   if (!is(control, "ControlEL")) {
#     stop("invalid 'control' specified")
#   }
#
#
#   # always one parameter regardless of whether we know mean or not
#
#   if (is.null(mean)) {
#     stop("sfesf")
#   } else {
#     # known mean provided
#     # check mean
#     if (isTRUE(
#       !is.numeric(mean) || !all(is.finite(mean)) || length(mean) != 1L
#     )) {
#       stop("'mean' must be a finite number")
#     }
#
#     # g
#     (mm - mean)^2L - par
#     # apply eval
#     par2 <- c(par, mean)
#     el <- eval_("var?", par2, mm, control@maxit_l, control@tol_l, control@th, w)
#
#   }
#
#
#   el <- eval_("mean", par, mm, control@maxit_l, control@tol_l, control@th, w)
#   new("EL",
#     optim = el$optim, logp = el$logp, logl = el$logl, loglr = el$loglr,
#     statistic = el$statistic, df = p,
#     pval = pchisq(el$statistic, df = p, lower.tail = FALSE), npar = p,
#     weights = w,
#     data = if (model) mm else matrix(NA_real_, nrow = 0L, ncol = 0L),
#     coefficients = est, method = "variance"
#   )
# }
