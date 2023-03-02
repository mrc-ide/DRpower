
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
#' @description Estimates power empirically via repeated simulation.
#' 
#' @details
#' Estimates power using the following approach:
#' \enumerate{
#'  \item Simulate repeatedly from the function \code{rbbinom_reparam()} using
#'  known values (e.g. a known "true" prevalence and intra-cluster correlation).
#'  \item Analyse data using \code{get_prevalence()} to determine the
#'  probability of being above \code{prev_thresh}.
#'  \item If this probability is above \code{rejection_threshold} then
#'  reject the null hypothesis, and encode this as a single correct conclusion.
#'  \item Count the number of simulations for which the correct conclusion is
#'  reached. This gives an estimate of empirical power, along with upper and
#'  lower 95% binomial CIs on the power via the method of Clopper and Pearson
#'  (1934).
#' }
#' 
#' @param N vector giving the number of samples obtained from each cluster.
#' @param prevalence assumed true prevalence of pfhrp2 deletions. Input as
#'   proportion between 0 and 1.
#' @param ICC assumed true intra-cluster correlation (ICC), between 0 and 1.
#' @param prev_thresh the threshold prevalence that we are testing against.
#' @param rejection_threshold the posterior probability of being above the
#'   prevalence threshold needs to be greater than \code{rejection_threshold} in
#'   order to reject the null hypothesis.
#' @param prior_prev_shape1,prior_prev_shape2,prior_ICC_shape1,prior_ICC_shape2
#'   parameters that dictate the shape of the Beta priors on prevalence and the
#'   ICC. Increasing the first shape parameter (e.g. \code{prior_prev_shape1})
#'   pushes the distribution towards 1, increasing the second shape parameter
#'   (e.g. \code{prior_prev_shape2}) pushes the distribution towards 0.
#'   Increasing both shape parameters squeezes the distribution towards the
#'   centre and therefore makes it narrower.
#' @param n_intervals the number of intervals used in the adaptive quadrature
#'   method. Increasing this value gives a more accurate representation of the
#'   true posterior, but comes at the cost of reduced speed.
#' @param reps number of times to repeat simulation per parameter combination.
#'
#' @references
#' Clopper, C.J. and Pearson, E.S., 1934. The use of confidence or fiducial
#' limits illustrated in the case of the binomial. Biometrika, 26, 404–413. doi:
#' 10.2307/2331986.
#'
#' @export

get_power <- function(N, prevalence = 0.10, ICC = 0.25,
                      prev_thresh = 0.05,
                      rejection_threshold = 0.95,
                      prior_prev_shape1 = 1, prior_prev_shape2 = 1,
                      prior_ICC_shape1 = 1, prior_ICC_shape2 = 3,
                      n_intervals = 20, reps = 1e2) {
  
  # check inputs
  assert_vector_pos_int(N)
  assert_single_bounded(prevalence)
  assert_single_bounded(ICC)
  assert_single_bounded(prev_thresh)
  assert_greq(prevalence, prev_thresh)
  assert_single_bounded(rejection_threshold)
  assert_single_bounded(prior_prev_shape1, left = 1, right = 1e3)
  assert_single_bounded(prior_prev_shape2, left = 1, right = 1e3)
  assert_single_bounded(prior_ICC_shape1, left = 1, right = 1e3)
  assert_single_bounded(prior_ICC_shape2, left = 1, right = 1e3)
  assert_single_pos_int(n_intervals)
  assert_greq(n_intervals, 5)
  assert_single_pos_int(reps)
  
  # simulate
  sim_correct <- rep(NA, reps)
  for (i in 1:reps) {
    n <- rbbinom_reparam(n_clust = length(N), N = N,
                         p = prevalence, rho = ICC)
    
    p_est <- get_prevalence(n = n, N = N,
                            prior_prev_shape1 = prior_prev_shape1,
                            prior_prev_shape2 = prior_prev_shape2,
                            prior_ICC_shape1 = prior_ICC_shape1,
                            prior_ICC_shape2 = prior_ICC_shape2,
                            prev_thresh = prev_thresh,
                            return_type = list(mean_on = FALSE,
                                               median_on = FALSE,
                                               CrI_on = FALSE,
                                               thresh_on = TRUE,
                                               full_on = FALSE),
                            n_intervals = n_intervals)
    
    sim_correct[i] <- (p_est$prob_above_threshold > rejection_threshold)
  }
  
  # get 95% CIs on power
  power_CI <- ClopperPearson(n_success = sum(sim_correct), n_total = reps, alpha = 0.05)
  ret <- data.frame(power = mean(sim_correct),
                    lower = power_CI["lower"],
                    upper = power_CI["upper"])
  rownames(ret) <- NULL
  return(ret)
}

