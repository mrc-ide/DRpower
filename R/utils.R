
#------------------------------------------------
#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#------------------------------------------------
#' @title Get pre-computed sample size tables
#'
#' @description Produce a sample size table giving the minimum sample size per
#'   cluster for given values of the ICC and the prevalence threshold against
#'   which we are comparing.
#' 
#' @details The function \code{get_power_threshold()} was run over a large range
#'   of parameter combinations and results were stored within the \code{df_sim}
#'   object (see \code{?df_sim}). These simulations were then used to produce
#'   minimum sample size estimates by linear interpolation that were stored
#'   within the \code{df_ss} object (see \code{?df_ss}). This function provides
#'   a simple way of querying the \code{df_ss} object for given parameter
#'   values.
#' 
#' @param prevalence the assumed true prevalence of pfhrp2/3 deletions in the
#'   domain. Allowed values are anything in \code{seq(0, 0.2, 0.01)}, including
#'   vectors of values.
#' @param ICC the assumed intra-cluster correlation. Allowed values are" \{0,
#'   0.01, 0.02, 0.05, 0.1, 0.2\}.
#' @param prev_thresh the prevalence threshold against which we are comparing.
#'   Allowed values are: \{0.05, 0.08, 0.1\}.
#'
#' @examples
#' get_sample_size_table()
#' 
#' @importFrom tidyr pivot_wider
#' @importFrom utils data
#' @export

get_sample_size_table <- function(prevalence = seq(0, 0.2, 0.01),
                                  ICC = 0.05,
                                  prev_thresh = 0.05) {
  
  # avoid "no visible binding" note
  df_ss <- n_clust <- N_opt <- NULL
  
  # check inputs
  assert_vector_bounded(prevalence)
  assert_in(prevalence, seq(0, 0.2, 0.01))
  assert_single_bounded(ICC)
  assert_in(ICC, c(0, 0.01, 0.02, 0.05, 0.1, 0.2))
  assert_single_bounded(prev_thresh)
  assert_in(prev_thresh, c(0.05, 0.08, 0.1))
  
  # load data
  data("df_ss", envir = environment())
  
  # rename variables for comparison
  p <- prevalence
  r <- ICC
  p_thresh <- prev_thresh
  
  # filter
  df_ss %>%
    filter(mapply(function(x) any(near(x, p)), prevalence)) %>%
    filter(mapply(function(x) any(near(x, r)), ICC)) %>%
    filter(mapply(function(x) any(near(x, p_thresh)), prev_thresh)) %>%
    select(prevalence, n_clust, N_opt) %>%
    pivot_wider(names_from = prevalence, values_from = N_opt)
}

#------------------------------------------------
# produce Clopper-Pearson upper and lower intervals
#' @noRd

ClopperPearson <- function(n_success, n_total, alpha = 0.05) {
  p_lower <- qbeta(p = alpha / 2, shape1 = n_success, shape2 = n_total - n_success + 1)
  p_upper <- qbeta(p = 1 - alpha / 2, shape1 = n_success + 1, shape2 = n_total - n_success)
  ret <- data.frame(lower = p_lower, upper = p_upper)
  return(ret)
}
