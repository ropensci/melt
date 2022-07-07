#' #' aaa
#' #'
#' #' aaa
#' #'
#' #' @param object An object
#' #' @export
#' aa <- function(object) {
#'   # take GLM object and run for loop to extract loo estimates
#'   # object is a valid GLM, no need to check
#'   # in order to run the loop, we need data for glm.fit
#'   #
#'   mm <- getDataMatrix(object)
#'   fit <- glm.fit(
#'     x = mm[, -1L],
#'     y = mm[, 1L],
#'     weights = weights(object),
#'     family = object@misc$family,
#'     control = object@misc$control,
#'     intercept = object@misc$intercept,
#'     singular.ok = FALSE
#'   )
#'   n <- 10000
#'   x <- vapply(seq_len(n), function(i) {
#'     ## ...
#'     }, integer(1))
#'     #'   return(fit$coefficients)
#' }
