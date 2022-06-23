
#------------------------------------------------
#' @title Estimate prevalence of pfhrp2 deletions from clustered data
#'
#' @description Takes raw counts of pfhrp2 deletions in multiple clusters (e.g.
#'   clinics) and estimates the overall prevalence. Uses a standard analysis
#'   approach in which the intra-cluster correlation (ICC) is first estimated,
#'   and this is used to compute an "effective sample size". Binomial confidence
#'   intervals are then calculated using this effective sample size rather than
#'   the raw sample size.
#'
#' @details The ICC is estimated via the ICCbin package using the stabilised
#'   estimate (\code{stab}) method proposed by Tamura and Young (1987). The
#'   design effect, \eqn{D_{eff}}{Deff}, is then calculated via the formula:
#'   \deqn{D_{eff} = 1 + (m - 1)r}{Deff = 1 + (m - 1)*r}
#'   where \eqn{m} is the number of clusters and \eqn{r} is the ICC. The
#'   effective sample size, \eqn{m_{eff}}{meff} is calculated via the formula:
#'   \deqn{m_{eff} = \frac{n m}{D_{eff}}}{meff = n*m / Deff}
#'   Finally, binomial confidence intervals are calculated using the
#'   Clopper-Pearson (1934) interval at a significance level \code{alpha}
#'   provided by the user (two-tailed).
#'   
#'   In some cases it is not possible to estimate the ICC using the stabilised
#'   method above - for example when there is a single cluster or when all
#'   counts are zero. In these cases a default value of 0.05 is used, which
#'   corresponds to a design effect of 2.8 for the reference situation of 37
#'   samples per clinic as described in the WHO master protocol.
#'
#' @references
#' Tamura, R.N. and Young, S.S., 1987. A stabilized moment estimator for the
#' beta-binomial distribution. Biometrics, pp.813-824.
#' 
#' Clopper, C.J. and Pearson, E.S., 1934. The use of confidence or fiducial
#' limits illustrated in the case of the binomial. Biometrika, 26, 404–413. doi:
#' 10.2307/2331986.
#'
#' @param pos_samples number of "positive" samples per cluster.
#' @param total_samples total sample size per cluster.
#' @param alpha the significance level of the credible interval - for example,
#'   use \code{alpha = 0.05} for a 95\% interval.
#' @param fix_ICC_estimate if defined, fix the ICC at this value rather than estimating
#'   from data.
#'
#' @importFrom ICCbin iccbin
#' @export

estimate_prevalence <- function(pos_samples, total_samples, alpha = 0.05, fix_ICC_estimate = NULL) {
  
  # check inputs
  assert_vector_pos_int(pos_samples)
  assert_single_pos_int(total_samples)
  assert_single_bounded(alpha)
  if (!is.null(fix_ICC_estimate)) {
    assert_single_bounded(fix_ICC_estimate)
  }
  
  # get basic data properties
  n <- length(pos_samples)
  
  # estimate ICC. Replace with 0.05 if method fails (this corresponds to a Deff
  # of 2.8 for the WHO master protocol situation of 10 clinics with 37 samples
  # per clinic)
  ICC_est <- 0.05
  if (!is.null(fix_ICC_estimate)) {
    ICC_est <- fix_ICC_estimate
  } else if (n > 1) {
    # get data into format required by ICCbin package
    df_dat <- data.frame(cid = as.factor(c(rep(1:n, times = pos_samples), rep(1:n, times = total_samples - pos_samples))),
                         y = c(rep(1, sum(pos_samples)), rep(0, sum(total_samples - pos_samples))))
    
    # estimate ICC
    ICC_raw <- try(suppressWarnings(ICCbin::iccbin(cid = 1, y = 2, data = df_dat, method = "stab", ci.type = NULL)), silent = TRUE)
    if (class(ICC_raw) != "try-error") {
      ICC_est <- ICC_raw$estimates$ICC
      if (ICC_est == "-") {
        ICC_est <- 0.05
      }
    }
  }
  
  # get design effect and effective sample size
  Dest <- 1 + (total_samples - 1) * ICC_est
  m_eff <- n * total_samples / Dest
  
  # get Clopper-Pearson interval on prevalence using effective sample size
  p_est <- mean(pos_samples / total_samples)
  ret <- ClopperPearson(n_success = p_est * m_eff, n_total = m_eff, alpha = alpha)
  
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
#'  \item Simulate repeatedly from the function \code{rbbinom_reparam()}
#'  using known values (e.g. a known "true" prevalence).
#'  \item Construct confidence intervals (CIs) on the prevalence from simualted
#'  data using \code{estimate_prevalence}.
#'  \item Make a decision as to whether prevalence is above the threshold,
#'  below, or inconclusive based on CIs.
#'  \item Count the number of simulations for which the correct conclusion is
#'  reached. This gives an estimate of empirical power, and upper and lower 95%
#'  binomial CIs on the power are produced using the method of Clopper and
#'  Pearson (1934).
#' }
#' 
#' @param clusters number of clusters.
#' @param total_samples total sample size per cluster.
#' @param prevalence assumed true prevalence of pfhrp2 deletions. Input as
#'   proportion between 0 and 1.
#' @param ICC,Deff assumed true intra-cluster correlation or design effect. Only
#'   one of these must be defined, the other must be \code{NULL}.
#' @param prevalence_threshold threshold used in decision-making. Input as
#'   proportion between 0 and 1.
#' @param alpha the significance level of the credible interval - for example,
#'   use \code{alpha = 0.05} for a 95\% interval.
#' @param reps number of times to repeat simulation per parameter combination.
#' @param fix_ICC_estimate if defined, fix the estimated ICC at this value
#'   rather than estimating from data. Note, this does not change the ICC used
#'   in simulation, which is still taken from \code{ICC}.
#'
#' @references
#' Clopper, C.J. and Pearson, E.S., 1934. The use of confidence or fiducial
#' limits illustrated in the case of the binomial. Biometrika, 26, 404–413. doi:
#' 10.2307/2331986.
#'
#' @importFrom rlang .data
#' @export

estimate_power <- function(clusters, total_samples, prevalence, ICC = NULL,
                           Deff = NULL, prevalence_threshold = 0.05, alpha = 0.05,
                           reps = 1e2, fix_ICC_estimate = NULL) {
  
  # check inputs
  assert_vector_pos_int(clusters)
  assert_vector_pos_int(total_samples)
  assert_vector_bounded(prevalence)
  if (!is.null(ICC)) {
    assert_vector_bounded(ICC)
  }
  if (!is.null(ICC)) {
    assert_vector_pos(Deff)
    assert_greq(Deff, 1.0)
  }
  assert_single_bounded(prevalence_threshold)
  assert_single_bounded(alpha)
  assert_single_pos_int(reps)
  if (!is.null(fix_ICC_estimate)) {
    assert_single_bounded(fix_ICC_estimate)
  }
  
  # one of ICC and Deff must be NULL
  if ((is.null(ICC) && is.null(Deff)) || (!is.null(ICC) && !is.null(Deff))) {
    stop("one (and only one) of ICC and Deff must be NULL")
  }
  
  # get ICC from Deff if needed
  if (is.null(ICC)) {
    ICC <- (Deff - 1) / (total_samples - 1)
  }
  
  # define dataframe holding parameter combinations, and that will eventually
  # hold results
  df_ret <- tidyr::expand_grid(clusters = clusters,
                               total_samples = total_samples,
                               prevalence = prevalence,
                               ICC = ICC) %>%
    dplyr::mutate(Deff = 1 + (clusters - 1) * ICC,
                  n_correct = 0)
  
  # loop through parameter combinations
  for (i in 1:nrow(df_ret)) {
    for (j in 1:reps) {
      
      # simulate samples from Beta-binomial model
      pos_samples <- rbbinom_reparam(n = df_ret$clusters[i],
                                     m = df_ret$total_samples[i],
                                     p = df_ret$prevalence[i],
                                     rho = df_ret$ICC[i])
      
      # estimate prevalence
      prev_CI <- estimate_prevalence(pos_samples = pos_samples,
                                     total_samples = df_ret$total_samples[i],
                                     alpha = alpha,
                                     fix_ICC_estimate = fix_ICC_estimate)
      
      # make call as to whether above threshold, below, or inconclusive 
      if (df_ret$prevalence[i] < prevalence_threshold) {
        if (prev_CI[2] < prevalence_threshold) {
          df_ret$n_correct[i] <- df_ret$n_correct[i] + 1
        }
      } else {
        if (prev_CI[1] > prevalence_threshold) {
          df_ret$n_correct[i] <- df_ret$n_correct[i] + 1
        }
      }
      
    }
  }
  
  # add power estimate and binomial CIs on power
  power_CI <- mapply(function(i) {
    ClopperPearson(n_success = df_ret$n_correct[i], n_total = reps, alpha = 0.05)
  }, 1:nrow(df_ret)) %>% t() %>% as.data.frame()
  df_ret <- df_ret %>%
    dplyr::mutate(power = .data$n_correct / reps * 100,
                  lower = power_CI$lower * 100,
                  upper = power_CI$upper * 100)
  
  return(df_ret)
}

