el_gbd <- function(data) {
  # take data.frame with 3 variables
  if (!is.data.frame(data) || length(data) != 3) {
    stop("data must be a data.frame with three variables")
  }

  # 2 variables for treatments and blocks & 1 variable for observations
  type <- sapply(data, class)
  if (!setequal(type, c("factor", "factor", "numeric"))) {
    stop("error")
  }

  # infer variable types and name them as "block", "trt, and "obs"
  obs.loc <- unname(which(type == "numeric"))
  levels_lengths <-
    sapply(data[, which(type == "factor")], function(x) length(levels(x)))
  block_loc <- setdiff(c(1, 2, 3), obs.loc)[which.max(levels_lengths)]
  # block_loc <- which.max(levels_lengths)
  trt_loc <- setdiff(c(1, 2, 3), obs.loc)[which.min(levels_lengths)]
  if (block_loc == trt_loc) {
    stop("error")
  }
  trt <- levels(data[, trt_loc])

  # incidence matrix
  c <- unclass(table(data[[block_loc]], data[[trt_loc]]))
  # data matrix
  x <- reshape(as.data.frame(data)[order(data[[trt_loc]]), ],
               idvar = names(data)[block_loc],
               timevar = names(data)[trt_loc],
               v.names = names(data)[obs.loc],
               direction = "wide")
  # remove block variable and convert to matrix
  x <- x[order(x[[block_loc]]), ]
  x[, names(data)[block_loc]] <- NULL
  x <- as.matrix(x)
  # replace NA with 0
  x[is.na(x)] <- 0
  dimnames(x) <- list(levels(data[, block_loc]), trt)

  gbd <- list("model_matrix" = x, "incidence_matrix" = c, "trt" = trt)
  class(gbd) <- c("gbd", "elmulttest")
  gbd
}
