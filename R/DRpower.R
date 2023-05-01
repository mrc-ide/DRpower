#------------------------------------------------
#' @title DRpower
#'
#' @description TODO
#'
#' @docType package
#' @name DRpower
NULL

#------------------------------------------------
# link to Rcpp
#' @useDynLib DRpower, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom magrittr %>%
NULL

#------------------------------------------------
# unload dll when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("DRpower", libpath)
}
