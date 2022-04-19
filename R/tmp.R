#' #' @export
#' el_lm2 <- function(formula, data, weights = NULL, na.action,
#'                    control = control_el(), model = TRUE, ...) {
#'   cl <- match.call()
#'   if (missing(data)) {
#'     data <- environment(formula)
#'   }
#'   mf <- match.call(expand.dots = FALSE)
#'   m <- match(c("formula", "data", "weights", "na.action"), names(mf), 0L)
#'   mf <- mf[c(1L, m)]
#'   mf$drop.unused.levels <- TRUE
#'   mf[[1L]] <- quote(stats::model.frame)
#'   mf <- eval(mf, parent.frame())
#'   mt <- attr(mf, "terms")
#'   y <- model.response(mf, "numeric")
#'   w <- as.vector(model.weights(mf))
#'   if (!is.null(w) && !is.numeric(w)) {
#'     stop("'weights' must be a numeric vector")
#'   }
#'   if (!is.null(w) && any(w < 0)) {
#'     stop("negative weights not allowed")
#'   }
#'   if (is.matrix(y)) {
#'     stop("'el_lm' does not support multiple responses")
#'   }
#'   if (is.empty.model(mt)) {
#'     x <- NULL
#'     out <- list(
#'       optim = list(), log.prob = numeric(), loglik = numeric(),
#'       statistic = numeric(), df = 0L, p.value = numeric(), npar = 0L,
#'       coefficients = numeric(), na.action = attr(mf, "na.action"),
#'       xlevels = .getXlevels(mt, mf), call = cl, terms = mt
#'     )
#'     class(out) <- c("el_lm", "el")
#'     return(out)
#'   } else {
#'     x <- model.matrix(mt, mf, NULL)
#'     z <- if (is.null(w)) {
#'       lm.fit(x, y, offset = NULL, singular.ok = FALSE, ...)
#'     } else {
#'       lm.wfit(x, y, w, offset = NULL, singular.ok = FALSE, ...)
#'     }
#'   }
#'   if (!inherits(control, "control_el") || !is.list(control)) {
#'     stop("invalid 'control' supplied")
#'   }
#'   intercept <- attr(mt, "intercept")
#'   mm <- cbind(y, x)
#'   p <- ncol(x)
#'   w <- check_weights(w, nrow(mm))
#'   out <- tmp_(
#'     mm, z$coefficients, intercept, control$maxit, control$maxit_l,
#'     control$tol, control$tol_l, control$th, control$nthreads, w
#'   )
#'   out
#' }
