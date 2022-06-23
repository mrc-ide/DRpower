
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
#'
#' @importFrom ICCbin iccbin
#' @export

estimate_prevalence <- function(pos_samples, total_samples, alpha = 0.05) {
  
  # get basic data properties
  n <- length(pos_samples)
  
  # estimate ICC. Replace with 0.05 if method fails (this corresponds to a Deff
  # of 2.8 for the WHO master protocol situation of 10 clinics with 37 samples
  # per clinic)
  ICC_est <- 0.05
  if (n > 1) {
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
