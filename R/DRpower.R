#------------------------------------------------
#' @title DRpower
#'
#' @description This package can be used in the design and/or analysis stages of
#'   Plasmodium falciparum pfhrp2/3 deletion prevalence studies. We assume that
#'   the study takes the form of a clustered prevalence survey, meaning the data
#'   consists of a numerator (number tested) and denominator (number of
#'   deletions found) over multiple clusters. We are interested in estimating
#'   the study-level prevalence, i.e. over all clusters, while accounting for
#'   the possibility of high intra-cluster correlation. The analysis approach
#'   uses a Bayesian random effects model to estimate prevalence and
#'   intra-cluster correlation. The approach to power analysis is
#'   simulation-based, running the analysis many times on simulated data and
#'   estimating empirical power. This method can be used to establish a minimum
#'   sample size required to achieve a given target power.
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
