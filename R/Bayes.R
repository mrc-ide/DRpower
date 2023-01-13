
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
# loglike_joint taking p as a vector
#' @noRd

loglike_joint_Vp <- function(k, m, p, rho,
                             prior_p_shape1 = 1, prior_p_shape2 = 1,
                             prior_rho_shape1 = 1, prior_rho_shape2 = 1) {
  n_clust <- length(k)
  
  # if the same value of k is present many times (for the same corresponding
  # value of m) then we should only compute the likelihood once and raise to the
  # appropriate power. Do this by determining the unique values of k and the
  # number of times they are present (the weight). Note, this only applies when
  # m is a single value over all clusters, as here there is a decent chance of
  # duplication
  if (length(m) == 1) {
    k_unique <- sort(unique(k))
    k_weights <- as.vector(table(k))
  } else {
    k_unique <- k
    k_weights <- rep(1, length(k))
  }
  
  # main use-case: rho between 0 and 1
  if ((rho != 0) && (rho != 1)) {
    
    # beta-binomial distribution
    alpha <- p * (1 - rho) / rho
    beta <- (1 - p) * (1 - rho) / rho
    ret_mat <- mapply(function(i) {
      k_weights[i] * extraDistr::dbbinom(x = k_unique[i], size = m, alpha = alpha, beta = beta, log = TRUE)
    }, seq_along(k_unique))
    
    # sum over rows
    ret <- rowSums(ret_mat)
    
    # deal with special cases of p=0 and p=1
    ret[p == 0] <- ifelse(all(k == 0), 0, -Inf)
    ret[p == 1] <- ifelse(all(k == m), 0, -Inf)
    
    # apply priors
    ret <- ret + dbeta(p, shape1 = prior_p_shape1, shape2 = prior_p_shape2, log = TRUE) +
      dbeta(rho, shape1 = prior_rho_shape1, shape2 = prior_rho_shape2, log = TRUE)
    
    return(ret)
    
  } else {
    
    if (rho == 0) {
      # ordinary binomial distribution
      ret_mat <- mapply(function(i) {
        k_weights[i] * dbinom(x = k_unique[i], size = m, prob = p, log = TRUE)
      }, seq_along(k_unique))
      
      ret <- rowSums(ret_mat)
      return(ret)
      
    } else {
      # clusters must have prevalence 0 or 1
      if (any((k != 0) & (k != m))) {
        ret <- rep(-Inf, length(p))
        return(ret)
      } else{
        n1 <- sum(k == m)
        ret <- dbinom(x = n1, size = n_clust, prob = p, log = TRUE)
        return(ret)
      }
    }
    
  }
}

#------------------------------------------------
# loglike_joint taking rho as a vector
#' @importFrom extraDistr dbbinom
#' @importFrom stats dbinom
#' @noRd

loglike_joint_Vrho <- function(k, m, p, rho,
                                prior_p_shape1 = 1, prior_p_shape2 = 1,
                                prior_rho_shape1 = 1, prior_rho_shape2 = 1) {
  n_clust <- length(k)
  
  # if the same value of k is present many times (for the same corresponding
  # value of m) then we should only compute the likelihood once and raise to the
  # appropriate power. Do this by determining the unique values of k and the
  # number of times they are present (the weight). Note, this only applies when
  # m is a single value over all clusters, as here there is a decent chance of
  # duplication
  if (length(m) == 1) {
    k_unique <- sort(unique(k))
    k_weights <- as.vector(table(k))
  } else {
    k_unique <- k
    k_weights <- rep(1, length(k))
  }
  
  # main use-case: p between 0 and 1
  if (p != 0 && p != 1) {
    
    # beta-binomial distribution
    alpha <- p * (1 - rho) / rho
    beta <- (1 - p) * (1 - rho) / rho
    ret_mat <- mapply(function(i) {
      k_weights[i] * extraDistr::dbbinom(x = k_unique[i], size = m, alpha = alpha, beta = beta, log = TRUE)
    }, seq_along(k_unique))
    
    # deal with special cases of rho=0 and rho=1
    ret_mat[rho == 0,] <- k_weights * dbinom(x = k_unique, size = m, prob = p, log = TRUE)
    ret_mat[rho == 1,] <- -Inf
    if (any(k_unique == 0)) {
      ret_mat[rho == 1, k_unique == 0] <- k_weights[k_unique == 0] * log(1 - p)
    }
    if (any(k_unique == m)) {
      ret_mat[rho == 1, k == m] <- k_weights[k_unique == m] * log(p)
    }
    
    # sum log-likelihood and apply prior
    ret <- rowSums(ret_mat) + dbeta(p, shape1 = prior_p_shape1, shape2 = prior_p_shape2, log = TRUE) +
      dbeta(rho, shape1 = prior_rho_shape1, shape2 = prior_rho_shape2, log = TRUE)
    
    return(ret)
    
  } else {
    
    # deal with special case of p=0 and p=1
    ret <- rep(-Inf, length(rho))
    if (p == 0) {
      if (all(k == 0)) {
        ret <- rep(0, length(rho))
      }
    } else {
      if (all(k == m)) {
        ret <- rep(0, length(rho))
      }
    }
    return(ret)
    
  }
}

###############################################
#                                             #
#   TRAPEZOIDAL-BASED METHODS                 #
#                                             #
###############################################

#------------------------------------------------
# marginal log-likelihood of rho, integrated over p via trapezoidal rule. This
# is a brute force method that is slow and inaccurate compared with more elegant
# approaches, but has the advantage of being very simple and so can be used to
# check other methods
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
# marginal log-likelihood of p, integrated over rho via trapezoidal rule (the
# same caveats apply as for loglike_marginal_rho_trap())
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
# get credible intervals for p via trapezoidal rule (the same caveats apply as
# for get_CrI_rho_trap())
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
#' @importFrom graphics abline
#' @noRd

loglike_marginal_rho_Gauss <- function(k, m, rho, n_intervals = 10, n_nodes = 5,
                                       precision_limit = 6*log(10),
                                       prior_p_shape1 = 1, prior_p_shape2 = 1,
                                       prior_rho_shape1 = 1, prior_rho_shape2 = 1,
                                       debug_on = FALSE) {
  
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
  
  # get lower and upper bounds
  bounds <- c(bound_lower$par, bound_upper$par)
  
  # get GQ nodes and weights for the standard [-1,1] interval, and transform to
  # apply to each of our new sub-intervals
  gq <- statmod::gauss.quad(n_nodes)
  int_width <- diff(bounds) / n_intervals
  int_midpoints <- seq(bounds[1] + 0.5*int_width, bounds[2] - 0.5*int_width, l = n_intervals)
  nodes_trans <- 0.5*int_width*gq$nodes + rep(int_midpoints, each = n_nodes)
  weights_trans <- rep(int_width/2 * gq$weights, times = n_intervals)
  
  # calculate loglikelihood and take out largest value as scaling factor
  ll <- loglike_joint_Vp(k = k, m = m, p = nodes_trans, rho = rho,
                         prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                         prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
  ll_max <- max(ll)
  
  # perform quadrature, take logs, and put scaling factor back in
  ret <- ll_max + log(sum(weights_trans * exp(ll - ll_max)))
  
  # debug distribution of p
  if (debug_on) {
    
    # get conditional distribution of p on a fine grid
    p_vec <- seq(bounds[1], bounds[2], l = 1001)
    ll_grid <- loglike_joint_Vp(k = k, m = m, p = p_vec, rho = rho,
                                prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                                prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
    
    # plot
    plot(p_vec, exp(ll_grid), type = "l", lwd = 3)
    points(nodes_trans, exp(ll), col = 3, pch = 20)
    abline(v = seq(bounds[1], bounds[2], l = n_intervals + 1), col = 2, lt = 2)
  }
  
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
#' @importFrom graphics abline
#' @noRd

loglike_marginal_p_Gauss <- function(k, m, p, n_intervals = 10, n_nodes = 5,
                                     precision_limit = 6*log(10),
                                     prior_p_shape1 = 1, prior_p_shape2 = 1,
                                     prior_rho_shape1 = 1, prior_rho_shape2 = 1,
                                     debug_on = FALSE) {
  
  # if m == 1 then likelihood simplifies to a Bernoulli in p. rho does not
  # feature, therefore marginalising over joint is equivalent to marginalising
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
  
  # get lower and upper bounds
  bounds <- c(bound_lower$par, bound_upper$par)
  
  # get GQ nodes and weights for the standard [-1,1] interval, and transform to
  # apply to each of our new sub-intervals
  gq <- statmod::gauss.quad(n_nodes)
  int_width <- diff(bounds) / n_intervals
  int_midpoints <- seq(bounds[1] + 0.5*int_width, bounds[2] - 0.5*int_width, l = n_intervals)
  nodes_trans <- 0.5*int_width*gq$nodes + rep(int_midpoints, each = n_nodes)
  weights_trans <- rep(int_width/2 * gq$weights, times = n_intervals)
  
  # calculate loglikelihood and take out largest value as scaling factor
  ll <- loglike_joint_Vrho(k = k, m = m, p = p, rho = nodes_trans,
                           prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                           prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
  ll_max <- max(ll)
  
  # perform quadrature, take logs, and put scaling factor back in
  ret <- ll_max + log(sum(weights_trans * exp(ll - ll_max)))
  
  # debug distribution of rho
  if (debug_on) {
    
    # get conditional distribution of rho on a fine grid
    rho_vec <- seq(bounds[1], bounds[2], l = 1001)
    ll_grid <- loglike_joint_Vrho(k = k, m = m, p = p, rho = rho_vec,
                                  prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                                  prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
    
    # plot
    plot(rho_vec, exp(ll_grid), type = "l", lwd = 3)
    points(nodes_trans, exp(ll), col = 3, pch = 20)
    abline(v = seq(bounds[1], bounds[2], l = n_intervals + 1), col = 2, lt = 2)
  }
  
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
#'   correlation coefficient (ICC) from clustered counts. By default these are
#'   95\% credible intervals, although the significance level can be altered by
#'   changing the \code{alpha} input value.
#'
#' @details There are two unknown quantities in the DRpower model: the
#'   prevalence and the ICC. This function integrates out the prevalence over a
#'   prior distribution to give the marginal distribution of the ICC. Then it
#'   returns the credible intervals of this distribution at a specified
#'   significance level.
#'
#' @param n,N the numerator (\code{n}) and denominator (\code{N}) per cluster.
#' @param alpha the significance level of the credible interval - for example,
#'   use \code{alpha = 0.05} for a 95\% interval.
#' @param prior_prev_shape1,prior_prev_shape2,prior_ICC_shape1,prior_ICC_shape2
#'   parameters that dictate the shape of the priors on prevalence and the ICC.
#'   Increasing the first shape parameter (e.g. \code{prior_p_shape1}) pushes
#'   the distribution towards 1, increasing the second shape parameter (e.g.
#'   \code{prior_p_shape2}) pushes the distribution towards 0. Increasing both
#'   shape parameters squeezes the distribution and makes it narrower.
#' @param debug_on For use in debugging, for advanced users only. If \code{TRUE}
#'   then produces a plot of the posterior distribution comparing various
#'   computational methods and approximations that are used internally. The
#'   black solid line is the distribution of prevalence marginalised over rho
#'   via trapezoidal rule on a fine grid (see \code{debug_grid} argument). This
#'   can be used as a sanity check, as for a fine enough grid this must tend
#'   towards the correct answer. The dashed red line is the distribution of
#'   prevalence marginalised using Gaussian quadrature (GQ) instead of
#'   trapezoidal rule. In the full method, GQ is only used at a smaller number
#'   of points indicated by green circles. The dashed green line is the
#'   quadratic interpolation of these points via Simpson's rule. The final
#'   credible interval estimates are calculated from the dashed green line. If
#'   the method is working correctly we should see a dashed green and red line
#'   that overlays a solid black line so closely that this almost looks like a
#'   black border.
#' @param debug_grid The number of equally spaced intervals used when performing
#'   marginalisation via trapezoidal rule in debugging (see \code{debug_on}
#'   argument).
#'
#' @importFrom stats optim qbeta
#' @export
#' @examples
#' # define three clusters with different number of observed positive counts.
#' # Try the default 95% CrI and a more stringent significance level
#' get_credible_ICC(n = c(2, 5, 4), N = 10)
#' get_credible_ICC(n = c(2, 5, 4), N = 10, alpha = 0.01)

get_credible_ICC <- function(n, N, alpha = 0.05,
                             prior_prev_shape1 = 1.0, prior_prev_shape2 = 1.0,
                             prior_ICC_shape1 = 1.0, prior_ICC_shape2 = 1.0,
                             debug_on = FALSE, debug_grid = 1e3) {
  
  # define a series of arguments that ordinarily would be part of the function
  # declaration, but I have made the design choice to put them here instead to
  # keep things as simple as possible for the user. Advanced users may want to
  # make their own version of this function and fiddle with these arguments
  precision_limit <- 6*log(10)
  n_intervals <- 40
  
  # check inputs
  assert_vector_pos_int(n)
  assert_single_pos_int(N)
  assert_single_bounded(alpha)
  assert_single_pos(prior_prev_shape1, zero_allowed = FALSE)
  assert_single_pos(prior_prev_shape2, zero_allowed = FALSE)
  assert_single_pos(prior_ICC_shape1, zero_allowed = FALSE)
  assert_single_pos(prior_ICC_shape2, zero_allowed = FALSE)
  assert_single_logical(debug_on)
  assert_single_pos_int(debug_grid, zero_allowed = FALSE)
  
  # if N == 1 then likelihood becomes independent of rho, therefore
  # the posterior equals the prior and we can return CrIs exactly
  if (N == 1) {
    CrI_lower <- qbeta(p = alpha / 2, shape1 = prior_ICC_shape1, shape2 = prior_ICC_shape2)
    CrI_upper <- qbeta(p = 1 - alpha / 2, shape1 = prior_ICC_shape1, shape2 = prior_ICC_shape2)
    return(c(lower = CrI_lower, upper = CrI_upper))
  }
  
  # get maximum of distribution
  ml <- optim(0.5, function(rho) {
    -loglike_marginal_rho_Gauss(k = n, m = N, rho = rho,
                                prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
  }, lower = 0, upper = 1, method = "Brent")
  
  # get lower and upper bounds at which distribution is a factor
  # exp(precision_limit) smaller than the maximum we just found. We will only
  # integrate over this interval as it contains nearly all the probability mass
  bound_lower <- optim(0, function(rho) {
    ll <- loglike_marginal_rho_Gauss(k = n, m = N, rho = rho,
                                     prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                     prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
    abs(ll + ml$value + precision_limit)
  }, lower = 0, upper = ml$par, method = "Brent")
  
  bound_upper <- optim(1, function(rho) {
    ll <- loglike_marginal_rho_Gauss(k = n, m = N, rho = rho,
                                     prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                     prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
    abs(ll + ml$value + precision_limit)
  }, lower = ml$par, upper = 1, method = "Brent")
  
  # define bounds
  bounds <- c(bound_lower$par, bound_upper$par)
  
  # define interval width
  node_interval <- diff(bounds) / n_intervals
  
  # define node points making up left, right, and midpoint of intervals
  node_pos <- seq(bounds[1], bounds[2], l = n_intervals + 1)
  node_left <- node_pos[-length(node_pos)]
  node_right <- node_pos[-1]
  node_mids <- (node_left + node_right) / 2
  
  # define integrand at these points. Note that ml$value is added back into this
  # expression to avoid underflow issues. This means the integral we are
  # computing is not the true marginal likelihood, but rather is multiplied by
  # an arbitrary scalar. This does not matter for our purposes of constructing
  # CrIs, as they will be in the same positions.
  f_node <- exp(loglike_marginal_rho_Gauss_Vrho(k = n, m = N, rho = node_pos,
                                                prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                                prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2) + ml$value)
  f_left <- f_node[-length(f_node)]
  f_right <- f_node[-1]
  f_mids <- exp(loglike_marginal_rho_Gauss_Vrho(k = n, m = N, rho = node_mids,
                                                prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                                prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2) + ml$value)
  
  # debug distribution of rho
  if (debug_on) {
    
    #code to debug internal GQ integration function
    # loglike_marginal_rho_Gauss(k = n, m = N, rho = 0.2,
    #                            prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
    #                            prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2,
    #                            debug_on = TRUE, precision_limit = 6*log(10))
    
    # get posterior distribution of rho using marginalisation via trapezoidal rule
    rho_vec <- seq(bound_lower$par, bound_upper$par, l = 201)
    ll_marginal_trap <- loglike_marginal_rho_trap_Vrho(k = n, m = N, rho = rho_vec, p_breaks = seq(0, 1, l = debug_grid),
                                                       prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                                       prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
    
    # get the same distribution using marginalisation via Gaussian quadrature
    ll_marginal_Gauss <- loglike_marginal_rho_Gauss_Vrho(k = n, m = N, rho = rho_vec,
                                                         prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                                         prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
    
    # get interpolation by Simpson's rule
    ll_Simpsons <- get_Simpsons_curve(rho_vec, node_left, node_right, f_left, f_mids, f_right)
    
    # plot
    plot(rho_vec, exp(ll_marginal_trap + ml$value), type = "l", lwd = 3)
    lines(rho_vec, exp(ll_marginal_Gauss + ml$value), col = 2)
    points(node_pos, f_node, col = 3, pch = 20)
    lines(rho_vec, ll_Simpsons, lty = 2, col = 3, lwd = 1)
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

#------------------------------------------------
#' @title Get credible intervals for the prevalence over clusters
#'
#' @description Produces lower and upper credible intervals on the prevalence
#'   from clustered counts. By default these are 95\% credible intervals,
#'   although the significance level can be altered by changing the \code{alpha}
#'   input value.
#'
#' @details There are two unknown quantities in the DRpower model: the
#'   prevalence and the ICC. This function integrates out the ICC over a prior
#'   distribution to give the marginal distribution of the prevalence. Then it
#'   returns the credible intervals of this distribution at a specified
#'   significance level.
#'
#' @param n,N the numerator (\code{n}) and denominator (\code{N}) per cluster.
#' @param alpha the significance level of the credible interval - for example,
#'   use \code{alpha = 0.05} for a 95\% interval.
#' @param prior_prev_shape1,prior_prev_shape2,prior_ICC_shape1,prior_ICC_shape2
#'   parameters that dictate the shape of the priors on prevalence and the ICC.
#'   Increasing the first shape parameter (e.g. \code{prior_p_shape1}) pushes
#'   the distribution towards 1, increasing the second shape parameter (e.g.
#'   \code{prior_p_shape2}) pushes the distribution towards 0. Increasing both
#'   shape parameters squeezes the distribution and makes it narrower.
#' @param debug_on For use in debugging, for advanced users only. If \code{TRUE}
#'   then produces a plot of the posterior distribution comparing various
#'   computational methods and approximations that are used internally. The
#'   black solid line is the distribution of prevalence marginalised over rho
#'   via trapezoidal rule on a fine grid (see \code{debug_grid} argument). This
#'   can be used as a sanity check, as for a fine enough grid this must tend
#'   towards the correct answer. The dashed red line is the distribution of
#'   prevalence marginalised using Gaussian quadrature (GQ) instead of
#'   trapezoidal rule. In the full method, GQ is only used at a smaller number
#'   of points indicated by green circles. The dashed green line is the
#'   quadratic interpolation of these points via Simpson's rule. The final
#'   credible interval estimates are calculated from the dashed green line. If
#'   the method is working correctly we should see a dashed green and red line
#'   that overlays a solid black line so closely that this almost looks like a
#'   black border.
#' @param debug_grid The number of equally spaced intervals used when performing
#'   marginalisation via trapezoidal rule in debugging (see \code{debug_on}
#'   argument).
#'
#' @importFrom stats optim
#' @importFrom graphics lines points
#' @export
#' @examples
#' # define three clusters with different number of observed positive counts.
#' # Try the default 95% CrI and a more stringent significance level
#' get_credible_prevalence(n = c(2, 5, 4), N = 10)
#' get_credible_prevalence(n = c(2, 5, 4), N = 10, alpha = 0.01)

get_credible_prevalence <- function(n, N, alpha = 0.05,
                                    prior_prev_shape1 = 1.0, prior_prev_shape2 = 1.0,
                                    prior_ICC_shape1 = 1.0, prior_ICC_shape2 = 1.0,
                                    debug_on = FALSE, debug_grid = 1e3) {
  
  # define a series of arguments that ordinarily would be part of the function
  # declaration but I have made the design choice to put them here instead to
  # keep things as simple as possible for the user. Advanced users may want to
  # make their own version of this function and fiddle with these arguments
  precision_limit <- 6*log(10)
  n_intervals <- 40
  
  # check inputs
  assert_vector_pos_int(n)
  assert_single_pos_int(N)
  assert_single_bounded(alpha)
  assert_single_pos(prior_prev_shape1, zero_allowed = FALSE)
  assert_single_pos(prior_prev_shape2, zero_allowed = FALSE)
  assert_single_pos(prior_ICC_shape1, zero_allowed = FALSE)
  assert_single_pos(prior_ICC_shape2, zero_allowed = FALSE)
  assert_single_logical(debug_on)
  assert_single_pos_int(debug_grid, zero_allowed = FALSE)
  
  # get maximum of distribution
  ml <- optim(0.5, function(p) {
    -loglike_marginal_p_Gauss(k = n, m = N, p = p,
                              prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                              prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
  }, lower = 0, upper = 1, method = "Brent")
  
  # get lower and upper bounds at which distribution is a factor
  # exp(precision_limit) smaller than the maximum we just found. We will only
  # integrate over this interval as it contains nearly all the probability mass
  bound_lower <- optim(0, function(p) {
    ll <- loglike_marginal_p_Gauss(k = n, m = N, p = p,
                                   prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                   prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
    abs(ll + ml$value + precision_limit)
  }, lower = 0, upper = ml$par, method = "Brent")
  
  bound_upper <- optim(1, function(p) {
    ll <- loglike_marginal_p_Gauss(k = n, m = N, p = p,
                                   prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                   prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
    abs(ll + ml$value + precision_limit)
  }, lower = ml$par, upper = 1, method = "Brent")
  
  # define bounds
  bounds <- c(bound_lower$par, bound_upper$par)
  
  # define interval width
  node_interval <- diff(bounds) / n_intervals
  
  # define node points making up left, right, and midpoint of intervals
  node_pos <- seq(bounds[1], bounds[2], l = n_intervals + 1)
  node_left <- node_pos[-length(node_pos)]
  node_right <- node_pos[-1]
  node_mids <- (node_left + node_right) / 2
  
  # define integrand at these points. Note that ml$value is added back into this
  # expression to avoid underflow issues. This means the integral we are
  # computing is not the true marginal likelihood, but rather is multiplied by
  # an arbitrary scalar. This does not matter for our purposes of constructing
  # CrIs, as they will be in the same positions.
  f_node <- exp(loglike_marginal_p_Gauss_Vp(k = n, m = N, p = node_pos,
                                            prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                            prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2) + ml$value)
  f_left <- f_node[-length(f_node)]
  f_right <- f_node[-1]
  f_mids <- exp(loglike_marginal_p_Gauss_Vp(k = n, m = N, p = node_mids,
                                            prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                            prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2) + ml$value)
  
  # debug distribution of p
  if (debug_on) {
    
    # code to debug internal GQ integration function
    # loglike_marginal_p_Gauss(k = n, m = N, p = 0.1,
    #                            prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
    #                            prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2,
    #                            debug_on = TRUE, precision_limit = 6*log(10))
    
    # get posterior distribution of rho using marginalisation via trapezoidal rule
    p_vec <- seq(bounds[1], bounds[2], l = 201)
    ll_marginal_trap <- loglike_marginal_p_trap_Vp(k = n, m = N, p = p_vec, rho_breaks = seq(0, 1, l = debug_grid),
                                                   prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                                   prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
    
    # get the same distribution using marginalisation via Gaussian quadrature
    ll_marginal_Gauss <- loglike_marginal_p_Gauss_Vp(k = n, m = N, p = p_vec,
                                                     prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                                                     prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
    
    # get interpolation by Simpson's rule
    ll_Simpsons <- get_Simpsons_curve(p_vec, node_left, node_right, f_left, f_mids, f_right)
    
    # plot
    plot(p_vec, exp(ll_marginal_trap + ml$value), type = "l", lwd = 3)
    lines(p_vec, exp(ll_marginal_Gauss + ml$value), col = 2)
    points(node_pos, f_node, col = 3, pch = 20)
    lines(p_vec, ll_Simpsons, lty = 2, col = 3, lwd = 1)
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
  
  #print(c(which_interval_lower, which_interval_upper))
  #print(interval_cumsum)
  #print(c(target_area_lower, target_area_upper))
  #print(c(area_remaining_lower, area_remaining_upper))
  
  
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

########################################
#                                      #
#   Rcpp-BASED METHODS                 #
#                                      #
########################################

#------------------------------------------------
#' @title Rcpp implementation of \code{get_credible_prevalence()}
#'
#' @description Equivalent to \code{get_credible_prevalence()} in functionality,
#'   but implemented in Rcpp for speed.
#'
#' @inheritParams get_credible_prevalence
#'
#' @export

get_credible_prevalence_fast <- function(n, N, alpha = 0.05,
                                         prior_prev_shape1 = 1.0, prior_prev_shape2 = 1.0,
                                         prior_ICC_shape1 = 1.0, prior_ICC_shape2 = 1.0) {
  
  # check inputs
  assert_vector_pos_int(n)
  assert_single_pos_int(N)
  assert_single_bounded(alpha)
  assert_single_pos(prior_prev_shape1, zero_allowed = FALSE)
  assert_single_pos(prior_prev_shape2, zero_allowed = FALSE)
  assert_single_pos(prior_ICC_shape1, zero_allowed = FALSE)
  assert_single_pos(prior_ICC_shape2, zero_allowed = FALSE)
  
  # define number of intervals and GQ nodes used in integration
  GQ_intervals <- 10
  GQ_nodes <- 5
  gq <- statmod::gauss.quad(GQ_nodes)
  
  # get arguments into list
  args_params <- list(n = n, N = rep(N, length(n)), alpha = alpha,
                      prior_prev_shape1 = prior_prev_shape1,
                      prior_prev_shape2 = prior_prev_shape2,
                      prior_ICC_shape1 = prior_ICC_shape1,
                      prior_ICC_shape2 = prior_ICC_shape2,
                      GQ_intervals = GQ_intervals,
                      GQ_nodes = GQ_nodes,
                      GQ_node_pos = 0.5*(gq$nodes + 1),
                      GQ_node_weights = gq$weights)
  
  # R functions to pass to C++
  args_functions <- list(solve_Simpsons_area = solve_Simpsons_area)
  
  # run efficient C++ function
  output_raw <- get_credible_prevalence_cpp(args_params, args_functions)
  
  # return
  names(output_raw) <- c("lower", "upper")
  return(output_raw)
}
