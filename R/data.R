#' Clothianidin concentration in maize plants
#'
#' A dataset summarizing field experiments result of seed treatments on
#'   clothianidin concentration.
#'
#' @format A data frame with 102 rows and 3 variables:
#' \describe{
#'   \item{blk}{New blocks constructed from original data. The format is 'days
#'     post planting_original block_year'.}
#'   \item{trt}{Seed treatment.}
#'   \item{clo}{Log transformed clothianidin concentration (µg).}
#' }
#' @usage data("clothianidin")
#'
#' @details The original data is provided by Alford and Krupke (2017). Only some
#'   of the shoot region observations are taken from the original data and
#'   processed for illustration.
#'
#' @source Alford, Adam, and Christian H. Krupke. 2017. “Translocation of the
#'   Neonicotinoid Seed Treatment Clothianidin in Maize.”
#'   Edited by Michael J. Stout. PLOS ONE 12 (3): e0173836.
#'   \doi{10.1371/journal.pone.0173836}.
"clothianidin"
