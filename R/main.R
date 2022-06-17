
#------------------------------------------------
# apply Simpson's rule to interpolate betwee a series of node x- and y-values.
# New values are calcualted at x_new, which is allowed to project beyond the
# range of the nodes, in which case values are extrapolated from the last
# quadratic fit (not recommended for anything other than very short distances)
#' @noRd

get_Simpsons_curve <- function(x_new, node_left, node_right, f_left, f_mid, f_right) {
  
  # get half-widths and midpoints of intervals
  h <- (node_right - node_left) / 2
  node_mid <- (node_left + node_right) / 2
  
  # define a series of coefficients. The equation of the quadratic line
  # interpolating between the three points of interest for each interval is:
  # a1*x^2 + a2*x + a3
  c1 <- 0.5 * f_left / h^2
  c2 <- -f_mid / h^2
  c3 <- 0.5 * f_right / h^2
  
  a1 <- c1 + c2 + c3
  a2 <- -c1*(node_mid + node_right) - c2*(node_left + node_right) - c3*(node_mid + node_left)
  a3 <- (c1 * node_mid * node_right) + (c2 * node_left * node_right) + (c3 * node_mid * node_left)
  
  # find which interval each of the x_new points corresponds to
  w <- as.numeric(cut(x_new, breaks = c(node_left[1], node_right)))
  
  # if none of the x_new are inside any interval then error
  if (all(is.na(w))) {
    stop("none of x_new inside given intervals")
  }
  
  # for x_new values outside intervals, extrapolate to the left and right
  if (any(is.na(w))) {
    first_good_value <- which(!is.na(w))[1]
    w[1:first_good_value] <- 1
  }
  if (any(is.na(w))) {
    last_good_value <- which(is.na(w))[1] - 1
    w[last_good_value:length(w)] <- w[last_good_value]
  }
  
  # a1*x^2 + a2*x + a3
  ret <- a1[w]*x_new^2 + a2[w]*x_new + a3[w]
  
  return(ret)
}

#------------------------------------------------
# solve Simpson's rule over the interval [a, b] (with corresponging y-values
# [fa, fb] and y-value at the midpoint equal to fm) to find the x-value at which
# the area under the curve equals target_area
#' @noRd

solve_Simpsons_area <- function(a, b, fa, fm, fb, target_area) {
  
  # get half-width and midpoint of interval
  h <- (b - a) / 2
  m <- (a + b) / 2
  
  # define a series of coefficients. The equation of the quadratic line
  # interpolating between the three points of interest is:
  # a1*x^2 + a2*x + a3
  c1 <- 0.5 * fa / h^2
  c2 <- -fm / h^2
  c3 <- 0.5 * fb / h^2
  
  a1 <- c1 + c2 + c3
  a2 <- -c1*(m + b) - c2*(a + b) - c3*(m + a)
  a3 <- (c1 * m * b) + (c2 * a * b) + (c3 * m * a)
  
  # define a new series of coefficients for the integral of this quadratic
  # equation from 0 to z. The solution can be written:
  # b1*z^3 + b2*z^2 + b3*z + b4
  b1 <- a1 / 3
  b2 <- a2 / 2
  b3 <- a3
  b4 <- -b1*a^3 - b2*a^2 - b3*a
  
  # solve for the z value at which the area under curve equals the target area
  roots <- polyroot(c(b4 - target_area, b3, b2, b1))
  
  # there should be exactly one root within the interval
  w <- which((Re(roots) > a) & (Re(roots) < b))
  if (length(w) != 1) {
    stop("could not find single root to cubic in Simpson's rule within the defined interval")
  }
  ret <- Re(roots)[w]
  
  return(ret)
}

#------------------------------------------------
# reparameterisation of the beta-binomial distribution in terms of a mean (p)
# and an intra-cluster correlation coefficient (rho). Deals with special cases
# that simplify to the binomial or bernoulli distributions
#' @importFrom extraDistr dbbinom
#' @importFrom stats dbinom
#' @noRd

dbbinom_reparam <- function(k, m, p, rho, log_on = TRUE) {
  if ((rho == 0) || (m == 1)) {
    # simplifies to binomial distribution
    ret <- dbinom(x = k, size = m, prob = p, log = TRUE)
    
  } else if (rho == 1) {
    # likelihood still positive whenever k == 0 or k == m
    n <- length(k)
    ret <- rep(-Inf, n)
    ret[k == 0] <- log(1 - p)
    ret[k == m] <- log(p)
    
  } else {
    if (p == 0) {
      # likelihood 1 whenever k == 0
      n <- length(k)
      ret <- rep(-Inf, n)
      ret[k == 0] <- 0
      
    } else if (p == 1) {
      # likelihood 1 whenever k == m
      n <- length(k)
      ret <- rep(-Inf, n)
      ret[k == m] <- 0
      
    } else {
      # beta-binomial distribution
      alpha <- p * (1 - rho) / rho
      beta <- (1 - p) * (1 - rho) / rho
      ret <- extraDistr::dbbinom(x = k, size = m, alpha = alpha, beta = beta, log = TRUE)
      
    }
  }
  if (!log_on) {
    ret <- exp(ret)
  }
  return(ret)
}

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
#'
#' @description Get lower and upper credible intervals on the intra-cluster
#'   correlation coefficient (ICC).
#'
#' @param k counts of positive samples.
#' @param m sample size per cluster.
#' @param alpha the significance level of the credible interval - for sample,
#'   use \code{alpha = 0.05} for 95.
#' @param prior_p_shape1,prior_p_shape2 parameters of the Beta prior on the
#'   prevalence.
#' @param prior_rho_shape1,prior_rho_shape2 parameters of the Beta prior on the
#'   ICC.
#'
#' @importFrom stats optim qbeta
#' @export
#' @examples
#' # TODO
#' print("foo")

get_credible_ICC <- function(k, m, alpha = 0.05,
                             prior_p_shape1 = 1.0, prior_p_shape2 = 1.0,
                             prior_rho_shape1 = 1.0, prior_rho_shape2 = 1.0) {
  
  # if m == 1 then likelihood becomes independent of rho, therefore the
  # posterior equals the prior and we can return CrIs exactly
  if (m == 1) {
    CrI_lower <- qbeta(p = alpha / 2, shape1 = prior_rho_shape1, shape2 = prior_rho_shape2)
    CrI_upper <- qbeta(p = 1 - alpha / 2, shape1 = prior_rho_shape1, shape2 = prior_rho_shape2)
    return(c(lower = CrI_lower, upper = CrI_upper))
  }
  
  # define a series of arguments that ordinarily would be part of the function
  # declaration. I have made the design choice to put them here instead to keep
  # things as simple as possible for the user. Advanced users may want to make
  # their own version of this function and fiddle with these arguments
  precision_limit <- 6*log(10)
  n_intervals <- 40
  
  # get maximum of distribution
  ml <- optim(0.5, function(rho) {
    -loglike_marginal_rho_Gauss(k = k, m = m, rho = rho,
                                prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                                prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
  }, lower = 0, upper = 1, method = "Brent")
  
  # get lower and upper bounds at which distribution is a factor
  # exp(precision_limit) smaller than the maximum we just found. We will only
  # integrate over this interval as it contains nearly all the probability mass
  bound_lower <- optim(0, function(rho) {
    ll <- loglike_marginal_rho_Gauss(k, m, rho = rho,
                                     prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                                     prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
    abs(ll + ml$value + precision_limit)
  }, lower = 0, upper = ml$par, method = "Brent")
  
  bound_upper <- optim(1, function(rho) {
    ll <- loglike_marginal_rho_Gauss(k = k, m = m, rho = rho,
                                     prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                                     prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
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
  f_node <- exp(loglike_marginal_rho_Gauss_Vrho(k = k, m = m, rho = node_pos,
                                                prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                                                prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2) + ml$value)
  f_left <- f_node[-length(f_node)]
  f_right <- f_node[-1]
  f_mids <- exp(loglike_marginal_rho_Gauss_Vrho(k = k, m = m, rho = node_mids,
                                                prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                                                prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2) + ml$value)
  
  # uncomment to plot distribution of rho
  #rho_vec <- seq(bound_lower$par, bound_upper$par, l = 1001)
  #z <- loglike_marginal_rho_Gauss_Vrho(k = k, m = m, rho = rho_vec,
  #                                     prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
  #                                     prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
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
#' @title Get credible intervals for prevalence
#'
#' @description Get lower and upper credible intervals on the prevalence over
#'   all clusters.
#'
#' @param k counts of positive samples.
#' @param m sample size per cluster.
#' @param alpha the significance level of the credible interval - for sample,
#'   use \code{alpha = 0.05} for 95.
#' @param prior_p_shape1,prior_p_shape2 parameters of the Beta prior on the
#'   prevalence.
#' @param prior_rho_shape1,prior_rho_shape2 parameters of the Beta prior on the
#'   ICC.
#'
#' @importFrom stats optim
#' @export
#' @examples
#' # TODO
#' print("foo")

get_credible_prevalence <- function(k, m, alpha = 0.05,
                                    prior_p_shape1 = 1.0, prior_p_shape2 = 1.0,
                                    prior_rho_shape1 = 1.0, prior_rho_shape2 = 1.0) {
  
  # define a series of arguments that ordinarily would be part of the function
  # declaration. I have made the design choice to put them here instead to keep
  # things as simple as possible for the user. Advanced users may want to make
  # their own version of this function and fiddle with these arguments
  precision_limit <- 6*log(10)
  n_intervals <- 40
  
  # get maximum of distribution
  ml <- optim(0.5, function(p) {
    -loglike_marginal_p_Gauss(k = k, m = m, p = p,
                              prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                              prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
  }, lower = 0, upper = 1, method = "Brent")
  
  # get lower and upper bounds at which distribution is a factor
  # exp(precision_limit) smaller than the maximum we just found. We will only
  # integrate over this interval as it contains nearly all the probability mass
  bound_lower <- optim(0, function(p) {
    ll <- loglike_marginal_p_Gauss(k = k, m = m, p = p,
                                   prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                                   prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
    abs(ll + ml$value + precision_limit)
  }, lower = 0, upper = ml$par, method = "Brent")
  
  bound_upper <- optim(1, function(p) {
    ll <- loglike_marginal_p_Gauss(k = k, m = m, p = p,
                                   prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                                   prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
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
  f_node <- exp(loglike_marginal_p_Gauss_Vp(k = k, m = m, p = node_pos,
                                            prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                                            prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2) + ml$value)
  f_left <- f_node[-length(f_node)]
  f_right <- f_node[-1]
  f_mids <- exp(loglike_marginal_p_Gauss_Vp(k = k, m = m, p = node_mids,
                                            prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
                                            prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2) + ml$value)
  
  # uncomment to plot distribution of p
  #p_vec <- seq(bound_lower$par, bound_upper$par, l = 1001)
  #z <- loglike_marginal_p_Gauss_Vp(k = k, m = m, p = p_vec,
  #                                 prior_p_shape1 = prior_p_shape1, prior_p_shape2 = prior_p_shape2,
  #                                 prior_rho_shape1 = prior_rho_shape1, prior_rho_shape2 = prior_rho_shape2)
  #plot(p_vec, exp(z + ml$value), type = 'l')
  
  # uncomment to overlay nodes and Simpson's rule approximation
  #points(node_pos, f_node, col = 3, pch = 20)
  #z2 <- get_Simpsons_curve(p_vec, node_left, node_right, f_left, f_mids, f_right)
  #lines(p_vec, z2, lty = 2, col = 3)
  
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

