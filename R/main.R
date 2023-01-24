
#------------------------------------------------
# produce Clopper-Pearson upper and lower intervals
#' @noRd

ClopperPearson <- function(n_success, n_total, alpha = 0.05) {
  p_lower <- qbeta(p = alpha / 2, shape1 = n_success, shape2 = n_total - n_success + 1)
  p_upper <- qbeta(p = 1 - alpha / 2, shape1 = n_success + 1, shape2 = n_total - n_success)
  ret <- c(lower = p_lower, upper = p_upper)
  return(ret)
}

#------------------------------------------------
#' @title Estimate power via simulation
#'
#' @description Estimates power empirically via repeated simulation. A range of
#'   input parameters can be specified (e.g. a range of sample sizes per
#'   cluster), power is then calculated on all combinations of input parameters.
#' 
#' @details
#' Estimates power using the following approach:
#' \enumerate{
#'  \item Simulate repeatedly from the function \code{rbbinom_reparam()} using
#'  known values (e.g. a known "true" prevalence and intra-cluster correlation).
#'  \item Construct credible intervals (CrIs) on the prevalence from simulated
#'  data using \code{get_credible_prevalence()}.
#'  \item Make a decision as to whether prevalence is above the threshold,
#'  below, or inconclusive based on CrIs.
#'  \item Count the number of simulations for which the correct conclusion is
#'  reached. This gives an estimate of empirical power, and upper and lower 95%
#'  binomial CIs on the power are produced using the method of Clopper and
#'  Pearson (1934).
#' }
#' 
#' @param cluster_size vector giving the number of samples obtained from each
#'   cluster.
#' @param prevalence assumed true prevalence of pfhrp2 deletions. Input as
#'   proportion between 0 and 1.
#' @param ICC assumed true intra-cluster correlation (ICC), between 0 and 1.
#' @param prevalence_threshold threshold used in decision-making. Input as
#'   proportion between 0 and 1.
#' @param alpha the significance level of the credible interval - for example,
#'   use \code{alpha = 0.05} for a 95\% interval.
#' @param reps number of times to repeat simulation per parameter combination.
#'
#' @references
#' Clopper, C.J. and Pearson, E.S., 1934. The use of confidence or fiducial
#' limits illustrated in the case of the binomial. Biometrika, 26, 404–413. doi:
#' 10.2307/2331986.
#'
#' @export

get_power <- function(cluster_size, prevalence, ICC,
                      prevalence_threshold = 0.05, alpha = 0.05, reps = 1e2) {
  
  # check inputs
  assert_vector_pos_int(cluster_size)
  assert_vector_bounded(prevalence)
  assert_vector_bounded(ICC)
  assert_single_bounded(prevalence_threshold)
  assert_single_bounded(alpha)
  assert_single_pos_int(reps)
  
  # simulate
  sim_correct <- rep(NA, reps)
  for (i in 1:reps) {
    n <- rbbinom_reparam(n = length(cluster_size), m = cluster_size,
                         p = prevalence, rho = ICC)
    p_est <- get_credible_prevalence(n = n, N = cluster_size, alpha = alpha)
    if (prevalence > prevalence_threshold) {
      sim_correct[i] <-  (p_est[1] > prevalence_threshold)
    } else {
      sim_correct[i] <-  (p_est[2] < prevalence_threshold)
    }
  }
  
  # get 95% CIs on power
  power_CI <- ClopperPearson(n_success = sum(sim_correct), n_total = reps, alpha = 0.05)
  ret <- c(power = mean(sim_correct), power_CI)
  
  return(ret)
}

