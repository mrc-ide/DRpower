
#------------------------------------------------
# joint probability of data multiplied by priors on p and rho
#' @importFrom stats dbeta
#' @noRd

loglike_joint <- function(k, m, p, rho,
                          prior_p_shape1 = 1, prior_p_shape2 = 1,
                          prior_rho_shape1 = 1, prior_rho_shape2 = 1) {
  
  dbeta(p, shape1 = prior_p_shape1, shape2 = prior_p_shape2, log = TRUE) +
    dbeta(rho, shape1 = prior_rho_shape1, shape2 = prior_rho_shape2, log = TRUE) +
    sum(dbbinom_reparam(k = k, m = m, p = p, rho = rho, log_on = TRUE))
}

#------------------------------------------------
# loglike_joint vectorized in terms of p
#' @noRd

loglike_joint_Vp <- Vectorize(loglike_joint, vectorize.args = "p")

#------------------------------------------------
# loglike_joint vectorized in terms of rho
#' @noRd

loglike_joint_Vrho <- Vectorize(loglike_joint, vectorize.args = "rho")

#------------------------------------------------
# marginal log-likelihood of rho, integrated over p via trapezoidal rule. This
# is a brute force method that is slow and inaccurate compared with more elegant
# approaches, but has the advantage of being very simple and so can be used to
# sanity check other methods
#' @noRd

loglike_marginal_rho_trap <- function(k, m, rho, p_breaks = seq(0, 1, l = 201),
                                      prior_p_shape1 = 1, prior_p_shape2 = 1,
                                      prior_rho_shape1 = 1, prior_rho_shape2 = 1) {
  
  # evaluate log-likelihood
  ll <- loglike_joint_Vp(k = k, m = m, p = p_breaks, rho = rho,
                         prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                         prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
  
  # special case if all -Inf
  if (all(ll == -Inf)) {
    return(-Inf)
  }
  
  # subtract max log-likelihood and exponentiate
  ll_max <- max(ll)
  integrand <- exp(ll - ll_max)
  
  # calculate AUC via trapezoidal rule
  trap <- sum(0.5 * (integrand[-length(integrand)] + integrand[-1]) * diff(p_breaks))
  
  # add back in constant and return
  ret <- ll_max + log(trap)
  return(ret)
}

#------------------------------------------------
# loglike_marginal_rho_trap vectorized in terms of rho
#' @noRd

loglike_marginal_rho_trap_Vrho <- Vectorize(loglike_marginal_rho_trap, vectorize.args = "rho")

#------------------------------------------------
# marginal log-likelihood of p, integrated over rho via trapezoidal rule
#' @noRd

loglike_marginal_p_trap <- function(k, m, p, rho_breaks = seq(0, 1, l = 201),
                                    prior_p_shape1 = 1, prior_p_shape2 = 1,
                                    prior_rho_shape1 = 1, prior_rho_shape2 = 1) {
  
  # evaluate log-likelihood
  ll <- loglike_joint_Vrho(k = k, m = m, p = p, rho = rho_breaks,
                           prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                           prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
  
  # special case if all -Inf
  if (all(ll == -Inf)) {
    return(-Inf)
  }
  
  # subtract max log-likelihood and exponentiate
  ll_max <- max(ll)
  integrand <- exp(ll - ll_max)
  
  # calculate AUC via trapezoidal rule
  trap <- sum(0.5 * (integrand[-length(integrand)] + integrand[-1]) * diff(rho_breaks))
  
  # add back in constant and return
  ret <- ll_max + log(trap)
  return(ret)
}

#------------------------------------------------
# loglike_marginal_p_trap vectorized in terms of p
#' @noRd

loglike_marginal_p_trap_Vp <- Vectorize(loglike_marginal_p_trap, vectorize.args = "p")

#------------------------------------------------
# get credible intervals for rho via trapezoidal rule. This is a brute force
# method that is slow and inaccurate compared with more elegant approaches, but
# has the advantage of being very simple and so can be used to check other
# methods
#' @noRd

get_CrI_rho_trap <- function(k, m, p_breaks = seq(0, 1, l = 201), rho_breaks = seq(0, 1, l = 201),
                             alpha = 0.05, prior_p_shape1 = 1, prior_p_shape2 = 1,
                             prior_rho_shape1 = 1, prior_rho_shape2 = 1) {
  
  # evaluate marginal log-likelihood over range of rho values
  ll <- loglike_marginal_rho_trap_Vrho(k = k, m = m, rho = rho_breaks, p_breaks = p_breaks,
                                       prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                                       prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
  
  # rescale and calculate cumulative integral via trapezoidal rule
  integrand <- exp(ll - max(ll))
  trap <- c(0, cumsum(0.5 * (integrand[-length(integrand)] + integrand[-1]) * diff(rho_breaks)))
  trap_max <- trap[length(trap)]
  
  # get lower and upper CrIs
  w_lower <- which((trap / trap_max) > (alpha / 2))[1]
  CrI_lower <- rho_breaks[w_lower]
  
  w_upper <- which((trap / trap_max) > (1 - alpha / 2))[1]
  CrI_upper <- rho_breaks[w_upper]
  
  # return
  return(c(lower = CrI_lower, upper = CrI_upper))
}

#------------------------------------------------
# get credible intervals for p via trapezoidal rule
#' @noRd

get_CrI_p_trap <- function(k, m, p_breaks = seq(0, 1, l = 201), rho_breaks = seq(0, 1, l = 201),
                           alpha = 0.05, prior_p_shape1 = 1, prior_p_shape2 = 1,
                           prior_rho_shape1 = 1, prior_rho_shape2 = 1) {
  
  # evaluate marginal log-likelihood over range of rho values
  ll <- loglike_marginal_p_trap_Vp(k = k, m = m, p = p_breaks, rho_breaks = rho_breaks,
                                   prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                                   prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
  
  # rescale and calculate cumulative integral via trapezoidal rule
  integrand <- exp(ll - max(ll))
  trap <- c(0, cumsum(0.5 * (integrand[-length(integrand)] + integrand[-1]) * diff(rho_breaks)))
  trap_max <- trap[length(trap)]
  
  # get lower and upper CrIs
  w_lower <- which((trap / trap_max) > (alpha / 2))[1]
  CrI_lower <- p_breaks[w_lower]
  
  w_upper <- which((trap / trap_max) > (1 - alpha / 2))[1]
  CrI_upper <- p_breaks[w_upper]
  
  # return
  return(c(lower = CrI_lower, upper = CrI_upper))
}


##############################################
#                                            #
#   QUADRATURE-BASED METHODS                 #
#                                            #
##############################################

#------------------------------------------------
# marginal log-likelihood of rho, integrated over p by Gaussian quadrature on a
# range of sub-intervals
#' @importFrom statmod gauss.quad
#' @noRd

loglike_marginal_rho_Gauss <- function(k, m, rho, n_intervals = 10, n_nodes = 5,
                                       precision_limit = 6*log(10),
                                       prior_p_shape1 = 1, prior_p_shape2 = 1,
                                       prior_rho_shape1 = 1, prior_rho_shape2 = 1) {
  
  # get maximum of distribution
  ml <- optim(0.5, function(p) {
    -loglike_joint(k, m, p = p, rho = rho,
                   prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                   prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
  }, lower = 0, upper = 1, method = "Brent")
  
  # get lower and upper bounds at which distribution is a factor
  # exp(precision_limit) smaller than the maximum we just found. We will only
  # integrate over this interval as it contains nearly all the probability mass
  bound_lower <- optim(0, function(p) {
    ll <- loglike_joint(k, m, p = p, rho = rho,
                        prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                        prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
    abs(ll + ml$value + precision_limit)
  }, lower = 0, upper = ml$par, method = "Brent")
  
  bound_upper <- optim(1, function(p) {
    ll <- loglike_joint(k, m, p = p, rho = rho,
                        prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                        prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
    abs(ll + ml$value + precision_limit)
  }, lower = ml$par, upper = 1, method = "Brent")
  
  # uncomment these lines to plot the distribution of p and get total area by
  # trapezoidal rule as a sanity check
  #rho <- 0.1
  #p_vec <- seq(bound_lower$par, bound_upper$par, l = 1001)
  #z <- loglike_joint_Vp(k = k, m = m, p = p_vec, rho = rho,
  #                      prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
  #                      prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
  #plot(p_vec, exp(z))
  #sum(0.5*(exp(z[-1]) + exp(z[-length(z)])) * diff(p_vec))
  
  # get GQ nodes and weights for the standard [-1,1] interval, and transform to
  # apply to each of our new sub-intervals
  gq <- gauss.quad(n_nodes)
  int_width <- (bound_upper$par - bound_lower$par) / n_intervals
  int_midpoints <- seq(bound_lower$par + 0.5*int_width, bound_upper$par - 0.5*int_width, l = n_intervals)
  nodes_trans <- 0.5*int_width*gq$nodes + rep(int_midpoints, each = n_nodes)
  weights_trans <- rep(int_width/2 * gq$weights, times = n_intervals)
  
  # calculate loglikelihood and take out largest value as scaling factor
  ll <- loglike_joint_Vp(k = k, m = m, p = nodes_trans, rho = rho,
                         prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                         prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
  ll_max <- max(ll)
  
  # perform quadrature, take logs, and put scaling factor back in
  ret <- ll_max + log(sum(weights_trans * exp(ll - ll_max)))
  
  return(ret)
}

#------------------------------------------------
# loglike_marginal_rho_Gauss vectorized in terms of rho
#' @noRd

loglike_marginal_rho_Gauss_Vrho <- Vectorize(loglike_marginal_rho_Gauss, vectorize.args = "rho")

#------------------------------------------------
# marginal log-likelihood of p, integrated over rho by Gaussian quadrature on a
# range of sub-intervals
#' @importFrom statmod gauss.quad
#' @noRd

loglike_marginal_p_Gauss <- function(k, m, p, n_intervals = 10, n_nodes = 5,
                                     precision_limit = 6*log(10),
                                     prior_p_shape1 = 1, prior_p_shape2 = 1,
                                     prior_rho_shape1 = 1, prior_rho_shape2 = 1) {
  
  # if m == 1 then likelihood simplifies to a Bernoulli in p. rho does not
  # feature, therefore marginalising over prior is equivalent to marginalising
  # over the prior
  if (m == 1) {
    ret <- sum(k == 1) * log(p) + sum(k == 0) * log(1 - p)
    return(ret)
  }
  
  # get maximum of distribution
  ml <- optim(0.5, function(rho) {
    -loglike_joint(k, m, p = p, rho = rho,
                   prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                   prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
  }, lower = 0, upper = 1, method = "Brent")
  
  # get lower and upper bounds at which distribution is a factor
  # exp(precision_limit) smaller than the maximum we just found. We will only
  # integrate over this interval as it contains nearly all the probability mass
  bound_lower <- optim(0, function(rho) {
    ll <- loglike_joint(k, m, p = p, rho = rho,
                        prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                        prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
    abs(ll + ml$value + precision_limit)
  }, lower = 0, upper = ml$par, method = "Brent")
  
  bound_upper <- optim(1, function(rho) {
    ll <- loglike_joint(k, m, p = p, rho = rho,
                        prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                        prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
    abs(ll + ml$value + precision_limit)
  }, lower = ml$par, upper = 1, method = "Brent")
  
  # uncomment these lines to plot the distribution of rho and get total area by
  # trapezoidal rule as a sanity check
  #rho_vec <- seq(bound_lower$par, bound_upper$par, l = 1001)
  #z <- loglike_joint_Vrho(k = k, m = m, p = p, rho = rho_vec,
  #                        prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
  #                        prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
  #plot(rho_vec, exp(z))
  #sum(0.5*(exp(z[-1]) + exp(z[-length(z)])) * diff(rho_vec))
  
  # get GQ nodes and weights for the standard [-1,1] interval, and transform to
  # apply to each of our new sub-intervals
  gq <- gauss.quad(n_nodes)
  int_width <- (bound_upper$par - bound_lower$par) / n_intervals
  int_midpoints <- seq(bound_lower$par + 0.5*int_width, bound_upper$par - 0.5*int_width, l = n_intervals)
  nodes_trans <- 0.5*int_width*gq$nodes + rep(int_midpoints, each = n_nodes)
  weights_trans <- rep(int_width/2 * gq$weights, times = n_intervals)
  
  # calculate loglikelihood and take out largest value as scaling factor
  ll <- loglike_joint_Vrho(k = k, m = m, p = p, rho = nodes_trans,
                           prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                           prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
  ll_max <- max(ll)
  
  # perform quadrature, take logs, and put scaling factor back in
  ret <- ll_max + log(sum(weights_trans * exp(ll - ll_max)))
  
  return(ret)
}

#------------------------------------------------
# loglike_marginal_p_Gauss vectorized in terms of p
#' @noRd

loglike_marginal_p_Gauss_Vp <- Vectorize(loglike_marginal_p_Gauss, vectorize.args = "p")

#------------------------------------------------
#' @title Get credible intervals for the intra-cluster correlation coefficient
#'   (ICC)
#'
#' @description Produces lower and upper credible intervals on the intra-cluster
#'   correlation coefficient (ICC). By default these are 95\% credible
#'   intervals, although the user can choose any significance level by altering
#'   the \code{alpha} input value.
#'
#' @details There are two unknown quantities in the DRpower model: the
#'   prevalence and the ICC. This function integrates out the prevalence over a
#'   prior distribution (we 'marginalise' out the prevalence) to give a
#'   distribution just in terms of the ICC. Then it returns the credible
#'   intervals of this distribution at a specified significance level.
#'
#' @param pos_samples number of "positive" samples per cluster.
#' @param total_samples total sample size per cluster.
#' @param alpha the significance level of the credible interval - for example,
#'   use \code{alpha = 0.05} for a 95\% interval.
#' @param prior_prev_shape1,prior_prev_shape2,prior_ICC_shape1,prior_ICC_shape2
#'   parameters that dictate the shape of the priors on prevalence and the ICC.
#'   Increasing the first shape parameter (e.g. \code{prior_p_shape1}) pushes
#'   the distribution to the right, increasing the second shape parameter (e.g.
#'   \code{prior_p_shape2}) pushes the distribution to the left. Increasing both
#'   shape parameters squeezes the distribution and makes it narrower.
#'
#' @importFrom stats optim qbeta
#' @export
#' @examples
#' # define three clusters with different number of observed positive counts.
#' # Try the default 95% CrI and a more stringent significance level
#' get_credible_ICC(pos_samples = c(2, 5, 4), total_samples = 10)
#' get_credible_ICC(pos_samples = c(2, 5, 4), total_samples = 10, alpha = 0.01)

get_credible_ICC <- function(pos_samples, total_samples, alpha = 0.05,
                             prior_prev_shape1 = 1.0, prior_prev_shape2 = 1.0,
                             prior_ICC_shape1 = 1.0, prior_ICC_shape2 = 1.0) {
  
  # if total_samples == 1 then likelihood becomes independent of rho, therefore
  # the posterior equals the prior and we can return CrIs exactly
  if (total_samples == 1) {
    CrI_lower <- qbeta(p = alpha / 2, shape1 = prior_ICC_shape1, shape2 = prior_ICC_shape2)
    CrI_upper <- qbeta(p = 1 - alpha / 2, shape1 = prior_ICC_shape1, shape2 = prior_ICC_shape2)
    return(c(lower = CrI_lower, upper = CrI_upper))
  }
  
  # define a series of arguments that ordinarily would be part of the function
  # declaration, but I have made the design choice to put them here instead to
  # keep things as simple as possible for the user. Advanced users may want to
  # make their own version of this function and fiddle with these arguments
  precision_limit <- 6*log(10)
  n_intervals <- 40
  
  # get maximum of distribution
  ml <- optim(0.5, function(rho) {
    -loglike_marginal_rho_Gauss(k = pos_samples, m = total_samples, rho = rho,
                                prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
  }, lower = 0, upper = 1, method = "Brent")
  
  # get lower and upper bounds at which distribution is a factor
  # exp(precision_limit) smaller than the maximum we just found. We will only
  # integrate over this interval as it contains nearly all the probability mass
  bound_lower <- optim(0, function(rho) {
    ll <- loglike_marginal_rho_Gauss(k = pos_samples, m = total_samples, rho = rho,
                                     prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                     prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
    abs(ll + ml$value + precision_limit)
  }, lower = 0, upper = ml$par, method = "Brent")
  
  bound_upper <- optim(1, function(rho) {
    ll <- loglike_marginal_rho_Gauss(k = pos_samples, m = total_samples, rho = rho,
                                     prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                     prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
    abs(ll + ml$value + precision_limit)
  }, lower = ml$par, upper = 1, method = "Brent")
  
  # define interval width
  node_interval <- (bound_upper$par - bound_lower$par) / n_intervals
  
  # define node points making up left, right, and midpoint of intervals
  node_pos <- seq(bound_lower$par, bound_upper$par, l = n_intervals + 1)
  node_left <- node_pos[-length(node_pos)]
  node_right <- node_pos[-1]
  node_mids <- (node_left + node_right) / 2
  
  # define integrand at these points. Note that ml$value is added back into this
  # expression to avoid underflow issues. This means the integral we are
  # computing is not the true marginal likelihood, but rather is multiplied by
  # an arbitrary scalar. This does not matter for our purposes of constructing
  # CrIs, as they will be in the same positions.
  f_node <- exp(loglike_marginal_rho_Gauss_Vrho(k = pos_samples, m = total_samples, rho = node_pos,
                                                prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                                prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2) + ml$value)
  f_left <- f_node[-length(f_node)]
  f_right <- f_node[-1]
  f_mids <- exp(loglike_marginal_rho_Gauss_Vrho(k = pos_samples, m = total_samples, rho = node_mids,
                                                prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                                prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2) + ml$value)
  
  # uncomment to plot distribution of rho
  #rho_vec <- seq(bound_lower$par, bound_upper$par, l = 1001)
  #z <- loglike_marginal_rho_Gauss_Vrho(k = pos_samples, m = total_samples, rho = rho_vec,
  #                                     prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
  #                                     prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
  #plot(rho_vec, exp(z + ml$value), type = 'l')
  
  # uncomment to overlay nodes and Simpson's rule approximation
  #points(node_pos, f_node, col = 3, pch = 20)
  #z2 <- get_Simpsons_curve(rho_vec, node_left, node_right, f_left, f_mids, f_right)
  #lines(rho_vec, z2, lty = 2, col = 3)
  
  # define area and cumulative area of each interval
  interval_area <- node_interval / 6 * (f_left + 4*f_mids + f_right)
  interval_cumsum <- cumsum(interval_area)
  
  # get target area of integral for lower and upper CrI
  target_area_lower <- alpha / 2 * sum(interval_area)
  target_area_upper <- (1 - alpha / 2) * sum(interval_area)
  
  # find which intervals these target areas fall within
  which_interval_lower <- which(interval_cumsum  > target_area_lower)[1]
  which_interval_upper <- which(interval_cumsum  > target_area_upper)[1]
  
  # find the remaining area in each interval
  area_remaining_lower <- target_area_lower - c(0, interval_cumsum)[which_interval_lower]
  area_remaining_upper <- target_area_upper - c(0, interval_cumsum)[which_interval_upper]
  
  # solve for lower and upper CrIs
  CrI_lower <- solve_Simpsons_area(a = node_left[which_interval_lower],
                                   b = node_right[which_interval_lower],
                                   fa = f_left[which_interval_lower],
                                   fm = f_mids[which_interval_lower],
                                   fb = f_right[which_interval_lower],
                                   target_area = area_remaining_lower)
  
  CrI_upper <- solve_Simpsons_area(a = node_left[which_interval_upper],
                                   b = node_right[which_interval_upper],
                                   fa = f_left[which_interval_upper],
                                   fm = f_mids[which_interval_upper],
                                   fb = f_right[which_interval_upper],
                                   target_area = area_remaining_upper)
  
  # return
  return(c(lower = CrI_lower, upper = CrI_upper))
}

#------------------------------------------------
#' @title Get credible intervals for the prevalence over clusters
#'
#' @description Produces lower and upper credible intervals on the prevalence
#'   from clustered counts. By default these are 95\% credible intervals,
#'   although the user can choose any significance level by altering the
#'   \code{alpha} input value.
#'
#' @details There are two unknown quantities in the DRpower model: the
#'   prevalence and the ICC. This function integrates out the ICC over a prior
#'   distribution (we 'marginalise' out the ICC) to give a distribution just in
#'   terms of the prevalence. Then it returns the credible intervals of this
#'   distribution at a specified significance level.
#'
#' @param pos_samples number of "positive" samples per cluster.
#' @param total_samples total sample size per cluster.
#' @param alpha the significance level of the credible interval - for example,
#'   use \code{alpha = 0.05} for a 95\% interval.
#' @param prior_prev_shape1,prior_prev_shape2,prior_ICC_shape1,prior_ICC_shape2
#'   parameters that dictate the shape of the priors on prevalence and the ICC.
#'   Increasing the first shape parameter (e.g. \code{prior_p_shape1}) pushes
#'   the distribution to the right, increasing the second shape parameter (e.g.
#'   \code{prior_p_shape2}) pushes the distribution to the left. Increasing both
#'   shape parameters squeezes the distribution and makes it narrower.
#' @param debug_on if \code{TRUE} then produces a plot of the distribution of
#'   the prevalence for use in debugging (for advanced users only). Open circles
#'   are the Gaussian quadrature estimates, the red solid line is the
#'   brute-force approach via trapezoidal rule, and green closed circles and
#'   lines are the nodes and quadratic interpolations used in Simpson's rule -
#'   these are used to produce final credible interval estimates.
#'
#' @importFrom stats optim
#' @importFrom graphics lines points
#' @export
#' @examples
#' # define three clusters with different number of observed positive counts.
#' # Try the default 95% CrI and a more stringent significance level
#' get_credible_prevalence(pos_samples = c(2, 5, 4), total_samples = 10)
#' get_credible_prevalence(pos_samples = c(2, 5, 4), total_samples = 10, alpha = 0.01)

get_credible_prevalence <- function(pos_samples, total_samples, alpha = 0.05,
                                    prior_prev_shape1 = 1.0, prior_prev_shape2 = 1.0,
                                    prior_ICC_shape1 = 1.0, prior_ICC_shape2 = 1.0,
                                    debug_on = FALSE) {
  
  # define a series of arguments that ordinarily would be part of the function
  # declaration but I have made the design choice to put them here instead to
  # keep things as simple as possible for the user. Advanced users may want to
  # make their own version of this function and fiddle with these arguments
  precision_limit <- 6*log(10)
  n_intervals <- 40
  
  # get maximum of distribution
  ml <- optim(0.5, function(p) {
    -loglike_marginal_p_Gauss(k = pos_samples, m = total_samples, p = p,
                              prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                              prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
  }, lower = 0, upper = 1, method = "Brent")
  
  # get lower and upper bounds at which distribution is a factor
  # exp(precision_limit) smaller than the maximum we just found. We will only
  # integrate over this interval as it contains nearly all the probability mass
  bound_lower <- optim(0, function(p) {
    ll <- loglike_marginal_p_Gauss(k = pos_samples, m = total_samples, p = p,
                                   prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                   prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
    abs(ll + ml$value + precision_limit)
  }, lower = 0, upper = ml$par, method = "Brent")
  
  bound_upper <- optim(1, function(p) {
    ll <- loglike_marginal_p_Gauss(k = pos_samples, m = total_samples, p = p,
                                   prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                   prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
    abs(ll + ml$value + precision_limit)
  }, lower = ml$par, upper = 1, method = "Brent")
  
  # define interval width
  node_interval <- (bound_upper$par - bound_lower$par) / n_intervals
  
  # define node points making up left, right, and midpoint of intervals
  node_pos <- seq(bound_lower$par, bound_upper$par, l = n_intervals + 1)
  node_left <- node_pos[-length(node_pos)]
  node_right <- node_pos[-1]
  node_mids <- (node_left + node_right) / 2
  
  # define integrand at these points. Note that ml$value is added back into this
  # expression to avoid underflow issues. This means the integral we are
  # computing is not the true marginal likelihood, but rather is multiplied by
  # an arbitrary scalar. This does not matter for our purposes of constructing
  # CrIs, as they will be in the same positions.
  f_node <- exp(loglike_marginal_p_Gauss_Vp(k = pos_samples, m = total_samples, p = node_pos,
                                            prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                            prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2) + ml$value)
  f_left <- f_node[-length(f_node)]
  f_right <- f_node[-1]
  f_mids <- exp(loglike_marginal_p_Gauss_Vp(k = pos_samples, m = total_samples, p = node_mids,
                                            prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                            prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2) + ml$value)
  
  # plot distribution of p
  if (debug_on) {
    
    # get marginal distribution of p
    p_vec <- seq(bound_lower$par, bound_upper$par, l = 1001)
    ll_marginal_Gauss <- loglike_marginal_p_Gauss_Vp(k = pos_samples, m = total_samples, p = p_vec,
                                                     prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                                     prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
    ll_marginal_trap <- loglike_marginal_p_trap_Vp(k = pos_samples, m = total_samples, p = p_vec, rho_breaks = seq(0, 1, l = 201),
                                                   prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                                   prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
    
    # get interpolation by Simpson's rule
    ll_Simpsons <- get_Simpsons_curve(p_vec, node_left, node_right, f_left, f_mids, f_right)
    
    # plot
    plot(p_vec, exp(ll_marginal_Gauss + ml$value))
    lines(p_vec, exp(ll_marginal_trap + ml$value), type = 'l', col = 2, lwd = 2)
    points(node_pos, f_node, col = 3, pch = 20)
    lines(p_vec, ll_Simpsons, lty = 2, col = 3, lwd = 2)
  }
  
  # define area and cumulative area of each interval
  interval_area <- node_interval / 6 * (f_left + 4*f_mids + f_right)
  interval_cumsum <- cumsum(interval_area)
  
  # get target area of integral for lower and upper CrI
  target_area_lower <- alpha / 2 * sum(interval_area)
  target_area_upper <- (1 - alpha / 2) * sum(interval_area)
  
  # find which intervals these target areas fall within
  which_interval_lower <- which(interval_cumsum  > target_area_lower)[1]
  which_interval_upper <- which(interval_cumsum  > target_area_upper)[1]
  
  # find the remaining area in each interval
  area_remaining_lower <- target_area_lower - c(0, interval_cumsum)[which_interval_lower]
  area_remaining_upper <- target_area_upper - c(0, interval_cumsum)[which_interval_upper]
  
  # solve for lower and upper CrIs
  CrI_lower <- solve_Simpsons_area(a = node_left[which_interval_lower],
                                   b = node_right[which_interval_lower],
                                   fa = f_left[which_interval_lower],
                                   fm = f_mids[which_interval_lower],
                                   fb = f_right[which_interval_lower],
                                   target_area = area_remaining_lower)
  
  CrI_upper <- solve_Simpsons_area(a = node_left[which_interval_upper],
                                   b = node_right[which_interval_upper],
                                   fa = f_left[which_interval_upper],
                                   fm = f_mids[which_interval_upper],
                                   fb = f_right[which_interval_upper],
                                   target_area = area_remaining_upper)
  
  # return
  return(c(lower = CrI_lower, upper = CrI_upper))
}
