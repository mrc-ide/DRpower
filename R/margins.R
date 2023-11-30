#------------------------------------------------
##' @name get_margins
##' @rdname get_margins
##'
#' @title Margin of error calculations when estimating prevalence from a
#'   clustered survey
#'
#' @description Calculate the expected margin of error when estimating
#'   prevalence from a clustered survey. Alternatively, calculate the sample
#'   size required to achieve a given target margin of error.
#'
#' @details A very common approach when constructing confidence intervals (CIs)
#'   from prevalence data is to use the Wald interval:
#'   
#'   \deqn{\hat{p} \pm z\sqrt{\frac{\hat{p}(1 - \hat{p})}{N}}}
#'   
#'   where \eqn{\hat{p}} is our estimate of the prevalence, \eqn{z} is the
#'   critical value of the normal distribution (\eqn{z=1.96} for a 95\%
#'   interval) and \eqn{N} is the sample size. When estimating prevalence from a
#'   clustered survey we need to modify this formula as follows:
#'   
#'   \deqn{\hat{p} \pm z\sqrt{\frac{\hat{p}(1 - \hat{p})}{Nc}(1 + (n - 1) r)}}
#'   
#'   where \eqn{\hat{p}} is the \emph{mean} prevalence over clusters, \eqn{c} is
#'   the number of clusters, and \eqn{r} is the intra-cluster correlation (ICC,
#'   a value between 0 and 1). The term to the right of the \eqn{\pm} symbol is
#'   called the \emph{margin of error} (MOE). The function \code{get_margin()}
#'   returns this value.
#'   
#'   We can also rearrange this formula to get the sample size (\eqn{N})
#'   required to reach any given MOE:
#'   
#'   \deqn{ N = \frac{ z^2p(1-p)(1-r) }{ cd^2 - z^2p(1-p)r } }
#'   
#'   where \eqn{d} is the desired MOE. The function
#'   \code{get_sample_size_margin()} returns this value. Note that in some cases
#'   it might not be possible to achieve the specified MOE for any finite sample
#'   size due to the ICC introducing too much variation, in which case this
#'   formula will return a negative value and the function will return an error.
#'   
#'   Although this is a very common approach, it has a number of weaknesses.
#'   First, notice that we sneakily replaced \eqn{\hat{p}} with \eqn{p} when
#'   moving to the sample size formula above. This implies that there is no
#'   uncertainty in our prevalence estimation, which is not true. Also note that
#'   the Wald interval assumes that the sampling distribution of our estimator
#'   is Gaussian, which is also not true. The difference between the Gaussian
#'   and the true distribution is particularly pronounced when prevalence is at
#'   the extremes of the range (near 0\% or 100\%). In these cases, the Wald
#'   interval can actually include values less than 0 or greater than 1, which
#'   is nonsensical.
#'   
#'   An arguably better approach is to construct CIs using the method of Clopper
#'   and Pearson (1934). This confidence interval guarantees that the false
#'   positive rate is \emph{at least} \code{alpha}, and in this sense is
#'   conservative. It is asymmetric and does not suffer from the problem of
#'   allowing values outside the [0,1] range. To make the Clopper-Pearson
#'   interval apply to a multi-cluster survey, we can use the idea of effective
#'   sample size, \eqn{N_e}:
#'   
#'   \deqn{ D_{eff} = 1 + (N - 1)r }
#'   \deqn{ N_e = \frac{N}{D_{eff}} }
#'   
#'   We then calculate the CI as normal but using \eqn{N_e} in place of \eqn{N}.
#'   The function \code{get_margin_CP()} returns the lower and upper MOE using
#'   the Clopper-Pearson interval, and the function
#'   \code{get_sample_size_margin_CP()} returns the corresponding sample size
#'   needed to achieve a certain MOE.
#'
#'   A third option is to use the DRpower model to estimate the credible
#'   interval of prevalence. See \code{?get_margin_Bayesian()} for how to
#'   estimate the margin of error under this method.
#' 
#' @returns the functions \code{get_margin()} and \code{get_margin_CP()} return
#'   the margin of error of the prevalence as a value between 0 and 1 (i.e. not
#'   as a \%). In the case of \code{get_margin_CP()} a different MOE is returned
#'   for the upper and lower limits.
#' 
#' @references
#' Clopper, C.J. and Pearson, E.S., 1934. The use of confidence or fiducial
#' limits illustrated in the case of the binomial. Biometrika, 26, 404–413. doi:
#' 10.2307/2331986.
NULL

##' @rdname get_margins
##' 
#' @param N the number of samples obtained from each cluster, assumed the same
#'   over all clusters.
#' @param n_clust the number of clusters.
#' @param prevalence the true prevalence of the marker in the population as a
#'   proportion between 0 and 1.
#' @param ICC assumed true intra-cluster correlation (ICC) between 0 and 1.
#' @param alpha the significance level of the CI.
#'
#' @importFrom stats qnorm
#'
#' @examples
#' get_margin(N = 60, n_clust = 3, prevalence = 0.2)
#' 
#' @export

get_margin <- function(N, n_clust, prevalence = 0.2, ICC = 0.05, alpha = 0.05) {
  
  # check inputs
  assert_single_pos_int(N, zero_allowed = FALSE)
  assert_single_pos_int(n_clust, zero_allowed = FALSE)
  assert_single_bounded(prevalence)
  assert_single_bounded(ICC)
  
  p <- prevalence
  ret <- qnorm(1 - alpha/2) * sqrt( p*(1 - p)/(N*n_clust)*(1 + (N - 1)*ICC) )
  
  return(ret)
}

##' @rdname get_margins
##' 
#' @param MOE the target margin of error.
#' @param n_clust the number of clusters.
#' @param prevalence the true prevalence of the marker in the population as a
#'   proportion between 0 and 1.
#' @param ICC assumed true intra-cluster correlation (ICC) between 0 and 1.
#' @param alpha the significance level of the CI.
#'
#' @importFrom stats qnorm
#'
#' @examples
#' get_sample_size_margin(MOE = 0.07, n_clust = 3, prevalence = 0.2, ICC = 0.01)
#' 
#' @export

get_sample_size_margin <- function(MOE, n_clust, prevalence = 0.2, ICC = 0.05, alpha = 0.05) {
  
  # check inputs
  assert_bounded(MOE, inclusive_left = FALSE, inclusive_right = FALSE)
  assert_single_pos_int(n_clust, zero_allowed = FALSE)
  assert_single_bounded(prevalence)
  assert_single_bounded(ICC)
  
  p <- prevalence
  z <- qnorm(1 - alpha/2)
  ret <- z^2*p*(1 - p)*(1 - ICC) / (n_clust*MOE^2 - z^2*p*(1 - p)*ICC)
  
  if (ret < 0) {
    stop("No finite sample size can achieve the desired margin of error. Consider decreasing the target MOE or changing other assumptions (e.g. number of clusters, ICC etc.)")
  }
  
  return(ceiling(ret))
}

##' @rdname get_margins
##' 
#' @param N the number of samples obtained from each cluster, assumed the same
#'   over all clusters.
#' @param n_clust the number of clusters.
#' @param prevalence the true prevalence of the marker in the population as a
#'   proportion between 0 and 1.
#' @param ICC assumed true intra-cluster correlation (ICC) between 0 and 1.
#' @param alpha the significance level of the CI.
#'
#'
#' @examples
#' get_margin_CP(N = 60, n_clust = 3, prevalence = 0.2)
#' 
#' @export

get_margin_CP <- function(N, n_clust, prevalence = 0.2, ICC = 0.05, alpha = 0.05) {
  
  # check inputs
  assert_single_pos_int(N, zero_allowed = FALSE)
  assert_single_pos_int(n_clust, zero_allowed = FALSE)
  assert_single_bounded(prevalence)
  assert_single_bounded(ICC)
  
  # get effective sample size
  Deff <- 1 + (N - 1)*ICC
  Ne <- N / Deff
  
  p <- prevalence
  moe_lower <- p - qbeta(p = alpha / 2, shape1 = Ne*p, shape2 = Ne*(1 - p) + 1)
  moe_upper <- qbeta(p = 1 - alpha / 2, shape1 = Ne*p + 1, shape2 = Ne*(1 - p)) - p
  ret <- c(lower = moe_lower, upper = moe_upper)
  
  return(ret)
}

##' @rdname get_margins
##' 
#' @param MOE the target margin of error.
#' @param n_clust the number of clusters.
#' @param prevalence the true prevalence of the marker in the population as a
#'   proportion between 0 and 1.
#' @param ICC assumed true intra-cluster correlation (ICC) between 0 and 1.
#' @param alpha the significance level of the CI.
#' @param N_max the largest value of \eqn{N} to consider.
#'
#'
#' @examples
#' get_sample_size_margin_CP(MOE = 0.14, n_clust = 3, prevalence = 0.2, ICC = 0.01)
#' 
#' @export

get_sample_size_margin_CP <- function(MOE, n_clust, prevalence = 0.2, ICC = 0.05,
                                      alpha = 0.05, N_max = 2e3) {
  
  # avoid "no visible binding" note
  lower <- upper <- NULL
  
  # check inputs
  assert_bounded(MOE, inclusive_left = FALSE, inclusive_right = FALSE)
  assert_single_pos_int(n_clust, zero_allowed = FALSE)
  assert_single_bounded(prevalence)
  assert_single_bounded(ICC)
  assert_single_pos_int(N_max)
  assert_greq(N_max, 10)
  
  # get MOE using CP method for all values of N up to N_max
  MOE_CP <- mapply(function(N) {
    get_margin_CP(N = N,
                  n_clust = n_clust,
                  prevalence = prevalence,
                  ICC = ICC,
                  alpha = alpha)
  }, 1:N_max) %>%
    t() %>%
    as.data.frame() %>%
    mutate(max = ifelse(lower > upper, lower, upper))
  
  # exit if no N achieves target
  if (!any(MOE_CP$max <= MOE)) {
    stop("No sample size up to N_max achieves the desired margin of error. Consider decreasing the target MOE, increasing N_max, or changing other assumptions (e.g. number of clusters, ICC etc.)")
  }
  
  # get smallest sample size to achieve target MOE
  ret <- which(MOE_CP$max <= MOE)[1]
  
  return(ret)
}

#------------------------------------------------
#' @title Margin of error calculations using the Bayesian DRpower model when
#'   estimating prevalence from a clustered survey
#'
#' @description As well as comparing against a threshold, the function
#'   \code{get_prevalence()} can be used to estimate a Bayesian credible
#'   interval (CrI) on the prevalence. We can estimate the margin of error (MOE)
#'   of this approach by simulation; repeatedly simulating datasets from known
#'   parameters, running \code{get_prevalence()} and measuring the MOE. Advantages of
#'   this method are that it accounts for uncertainty in the ICC, accounts for
#'   uncertainty in prevalence estimation by using simulation, and allows for
#'   incorporation of prior information.
#' 
#' @details
#' Estimates MOE using the following approach:
#' \enumerate{
#'  \item Simulate data via the function \code{rbbinom_reparam()} using known
#'  values (e.g. a known "true" prevalence and intra-cluster correlation).
#'  \item Analyse data using \code{get_prevalence()} to determine the
#'  upper and lower limits of the credible interval.
#'  \item Repeat steps 1-2 many times to obtain the distribution of upper and
#'  lower limits Return either the full distribution, or the mean of the
#'  distribution along with upper and lower 95\% CIs.
#' }
#' Note that we have not implemented a function that returns the sample size
#' needed to achieve a given MOE under the Bayesian model because this would
#' require repeated simulation over different values of \code{N}, which is
#' computationally costly. The appropriate value can be established manually if
#' needed by running \code{get_margin_Bayesian()} for different sample sizes.
#' 
#' @inheritParams get_power
#' @param prevalence assumed true prevalence of pfhrp2/3 deletions as a
#'   proportion between 0 and 1.
#' @param alpha the significance level of the credible interval - for example,
#'   use \code{alpha = 0.05} for a 95\% interval. See also \code{CrI_type}
#'   argument for how this is calculated.
#' @param CrI_type which method to use when computing credible intervals.
#'   Options are "ETI" (equal-tailed interval) or "HDI" (high-density interval).
#'   The ETI searches a distance \code{alpha/2} from either side of the [0,1]
#'   interval. The HDI method returns the narrowest interval that subtends a
#'   proportion \code{1-alpha} of the distribution. The HDI method is used by
#'   default as it guarantees that the MAP estimate is within the credible
#'   interval, which is not always the case for the ETI.
#' @param return_full if \code{TRUE} then return the complete distribution of
#'   lower and upper CrI limits in a data.frame. If \code{FALSE} (the default)
#'   return a summary including the mean and 95\% CI of these limits.
#'
#' @returns If \code{return_full = FALSE} (the default) returns a data.frame
#'   where the lower MOE is in the first row, and the upper MOE is in the second
#'   row. The first column gives the point estimate, and the subsequent columns
#'   give the 95\% CI on this estimate. If \code{return_full = TRUE} then
#'   returns a complete data.frame of all lower and upper MOE realisations over
#'   simulations.
#'
#' @importFrom stats var
#'
#' @examples
#' get_margin_Bayesian(N = c(120, 90, 150), prevalence = 0.15, ICC = 0.01 , reps = 1e2)
#' 
#' @export

get_margin_Bayesian <- function(N, prevalence = 0.2, ICC = 0.05, alpha = 0.05,
                                prior_prev_shape1 = 1,
                                prior_prev_shape2 = 1,
                                prior_ICC_shape1 = 1,
                                prior_ICC_shape2 = 9,
                                CrI_type = "HDI",
                                n_intervals = 20, round_digits = 2,
                                reps = 100, use_cpp = TRUE,
                                return_full = FALSE, silent = FALSE) {
  
  # avoid "no visible binding" note
  n <- NULL
  
  # check inputs
  assert_vector_pos_int(N)
  assert_single_bounded(prevalence, inclusive_left = FALSE, inclusive_right = FALSE)
  assert_single_bounded(ICC)
  assert_bounded(alpha, inclusive_left = FALSE, inclusive_right = FALSE)
  assert_single_bounded(prior_prev_shape1, left = 1e-3, right = 1e3)
  assert_single_bounded(prior_prev_shape2, left = 1e-3, right = 1e3)
  assert_single_bounded(prior_ICC_shape1, left = 1e-3, right = 1e3)
  assert_single_bounded(prior_ICC_shape2, left = 1e-3, right = 1e3)
  assert_single_string(CrI_type)
  assert_in(CrI_type, c("ETI", "HDI"))
  assert_single_pos_int(n_intervals)
  assert_greq(n_intervals, 5)
  assert_single_pos_int(round_digits)
  assert_single_pos_int(reps)
  assert_single_logical(use_cpp)
  assert_single_logical(return_full)
  assert_single_logical(silent)
  
  # draw data for all simulations. This allows us to group identical datasets,
  # which can save time in some situations (e.g. low prevalence where 0 counts
  # are common)
  l_n <- list()
  for (i in 1:reps) {
    l_n[[i]] <- data.frame(N = N,
                           n = rbbinom_reparam(n_clust = length(N), N = N,
                                               p = prevalence, rho = ICC)) %>%
      arrange(N, n)
  }
  
  # group duplicates
  l_u <- unique(l_n)
  l_w <- tabulate(match(l_n, l_u))
  
  # make progress bar
  pb <- progress_estimated(length(l_u))
  
  # simulate
  sim_df <- data.frame(lower = rep(NA, length(l_u)),
                       upper = NA)
  for (i in seq_along(l_u)) {
    if (!silent) {
      update_progress(pb)
    }
    
    p_est <- get_prevalence(n = l_u[[i]]$n,
                            N = l_u[[i]]$N,
                            alpha = alpha,
                            prev_thresh = 0.05,
                            prior_prev_shape1 = prior_prev_shape1,
                            prior_prev_shape2 = prior_prev_shape2,
                            prior_ICC_shape1 = prior_ICC_shape1,
                            prior_ICC_shape2 = prior_ICC_shape2,
                            MAP_on = FALSE,
                            post_mean_on = FALSE,
                            post_median_on = FALSE,
                            post_CrI_on = TRUE,
                            post_thresh_on = FALSE,
                            post_full_on = FALSE,
                            CrI_type = CrI_type,
                            n_intervals = n_intervals,
                            round_digits = round_digits,
                            use_cpp = use_cpp)
    
    sim_df$lower[i] <- p_est$CrI_lower
    sim_df$upper[i] <- p_est$CrI_upper
  }
  
  # unroll identical datasets
  sim_df <- sim_df[rep(seq_along(l_w), times = l_w),]
  
  # option to return complete distribution
  if (return_full) {
    return(sim_df)
  }
  
  # calculate mean and 95% CI on limits
  m <- colMeans(sim_df)
  v <- apply(sim_df, 2, var)
  d <- 1.96*sqrt(v / reps)
  ret <- data.frame(estimate = m,
                    CI_2.5 = m - d,
                    CI_97.5 = m + d)
  
  return(ret)
}
