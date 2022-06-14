
#------------------------------------------------
# reparameterisation of the beta-binomial distribution in terms of a mean (p)
# and an intra-cluster correlation coefficient (rho). Deals with special cases
# in which rho is 0 or 1, corresponding to the binomial and bernoulli
# distributions, respectively.
#' @importFrom extraDistr dbbinom
#' @importFrom stats dbinom
#' @noRd

dbbinom_reparam <- function(k, m, p, rho, log_on = TRUE) {
  if (rho == 0) {
    # simplifies to binomial distribution
    ret <- dbinom(x = k, size = m, prob = p, log = TRUE)
    
  } else if (rho == 1) {
    # likelihood still finite whenever k == 0 or k == m
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

loglike_joint <- function(k, m, p, rho) {
  dbeta(p, shape1 = 1, shape2 = 1, log = TRUE) +
    dbeta(rho, shape1 = 1, shape2 = 1, log = TRUE) +
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
# marginal log-likelihood of rho, integrated over p via trapezoidal rule
#' @noRd

loglike_marginal_rho_trap <- function(k, m, rho, p_breaks = seq(0, 1, l = 201)) {
  
  # evaluate log-likelihood
  ll <- loglike_joint_Vp(k = k, m = m, p = p_breaks, rho = rho)
  
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

loglike_marginal_p_trap <- function(k, m, p, rho_breaks = seq(0, 1, l = 201)) {
  
  # evaluate log-likelihood
  ll <- loglike_joint_Vrho(k = k, m = m, p = p, rho = rho_breaks)
  
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
# get credible intervals for rho via trapezoidal rule
#' @noRd

get_CrI_rho_trap <- function(k, m, p_breaks = seq(0, 1, l = 201), rho_breaks = seq(0, 1, l = 201),
                             alpha = 0.05) {
  
  # evaluate marginal log-likelihood over range of rho values
  ll <- loglike_marginal_rho_trap_Vrho(k = k, m = m, rho = rho_breaks, p_breaks = p_breaks)
  
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
                           alpha = 0.05) {
  
  # evaluate marginal log-likelihood over range of rho values
  ll <- loglike_marginal_p_trap_Vp(k = k, m = m, p = p_breaks, rho_breaks = rho_breaks)
  
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

#------------------------------------------------
# marginal log-likelihood of rho, integrated over p via gaussian quadrature
#' @importFrom statmod gauss.quad.prob
#' @noRd

loglike_marginal_rho_Gauss <- function(k, m, rho, p_nodes = 20) {
  
  # get Gaussian quadrature nodes and weights
  gq <- statmod::gauss.quad.prob(p_nodes)
  
  # evaluate log-likelihood at Gaussian quadrature node locations
  ll <- loglike_joint_Vp(k = k, m = m, p = gq$nodes, rho = rho)
  
  # special case if all -Inf
  if (all(ll == -Inf)) {
    return(-Inf)
  }
  
  # subtract max log-likelihood and exponentiate
  ll_max <- max(ll)
  integrand <- exp(ll - ll_max)
  
  # calculate AUC by Gaussian quadrature method
  gq_res <- sum(gq$weights * integrand)
  
  # add back in constant and return
  ret <- ll_max + log(gq_res)
  return(ret)
}

#------------------------------------------------
# loglike_marginal_rho_Gauss vectorized in terms of rho
#' @noRd

loglike_marginal_rho_Gauss_Vrho <- Vectorize(loglike_marginal_rho_Gauss, vectorize.args = "rho")

#------------------------------------------------
# marginal log-likelihood of p, integrated over rho via gaussian quadrature
#' @importFrom statmod gauss.quad.prob
#' @noRd

loglike_marginal_p_Gauss <- function(k, m, p, rho_nodes = 20) {
  
  # get Gaussian quadrature nodes and weights
  gq <- statmod::gauss.quad.prob(rho_nodes)
  
  # evaluate log-likelihood at Gaussian quadrature node locations
  ll <- loglike_joint_Vrho(k = k, m = m, p = p, rho = gq$nodes)
  
  # special case if all -Inf
  if (all(ll == -Inf)) {
    return(-Inf)
  }
  
  # subtract max log-likelihood and exponentiate
  ll_max <- max(ll)
  integrand <- exp(ll - ll_max)
  
  # calculate AUC by Gaussian quadrature method
  gq_res <- sum(gq$weights * integrand)
  
  # add back in constant and return
  ret <- ll_max + log(gq_res)
  return(ret)
}

#------------------------------------------------
# loglike_marginal_p_Gauss vectorized in terms of p
#' @noRd

loglike_marginal_p_Gauss_Vp <- Vectorize(loglike_marginal_p_Gauss, vectorize.args = "p")

#------------------------------------------------
# integrate likelihood up to rho_upper by Gaussian quadrature
#' @importFrom statmod gauss.quad.prob
#' @noRd

loglike_rho_upper_Gauss <- function(k, m, rho_upper, p_nodes = 20, rho_nodes = 20, alpha = 0.05) {
  
  # get Gaussian quadrature nodes and weights
  gq <- statmod::gauss.quad.prob(rho_nodes, u = rho_upper)
  
  # evaluate marginal log-likelihood at nodes
  ll <- loglike_marginal_rho_Gauss_Vrho(k = k, m = m, rho = gq$nodes, p_nodes = p_nodes)
  ll_max <- max(ll)
  
  # rescale and integrate up to limit by Gaussian quadrature
  integrand <- exp(ll - ll_max)
  gq_res <- rho_upper * sum(gq$weights * integrand)
  
  # add back in constant and return
  ret <- ll_max + log(gq_res)
  return(ret)
}

#------------------------------------------------
# loglike_rho_upper_Gauss vectorized in terms of rho_upper
#' @noRd

loglike_rho_upper_Gauss_Vrho <- Vectorize(loglike_rho_upper_Gauss, vectorize.args = "rho_upper")

#------------------------------------------------
# integrate likelihood up to p_upper by Gaussian quadrature
#' @importFrom statmod gauss.quad.prob
#' @noRd

loglike_p_upper_Gauss <- function(k, m, p_upper, p_nodes = 20, rho_nodes = 20, alpha = 0.05) {
  
  # get Gaussian quadrature nodes and weights
  gq <- statmod::gauss.quad.prob(p_nodes, u = p_upper)
  
  # evaluate marginal log-likelihood at nodes
  ll <- loglike_marginal_p_Gauss_Vp(k = k, m = m, p = gq$nodes, rho_nodes = rho_nodes)
  ll_max <- max(ll)
  
  # rescale and integrate up to limit by Gaussian quadrature
  integrand <- exp(ll - ll_max)
  gq_res <- p_upper * sum(gq$weights * integrand)
  
  # add back in constant and return
  ret <- ll_max + log(gq_res)
  return(ret)
}

#------------------------------------------------
# loglike_p_upper_Gauss vectorized in terms of p
#' @noRd

loglike_p_upper_Gauss_Vp <- Vectorize(loglike_p_upper_Gauss, vectorize.args = "p_upper")

#------------------------------------------------
#' @title Get credible intervals for the intra-cluster correlation coefficient
#'
#' @description Get lower and upper credible intervals on the intra-cluster
#'   correlation coefficient (ICC).
#'
#' @param k counts of positive samples.
#' @param m sample size per cluster.
#' @param alpha the significance level of the credible interval - for sample,
#'   use \code{alpha = 0.05} for 95
#' @param p_nodes,rho_nodes number of nodes to use in the Gaussian quadrature
#'   approximate integral method.
#'
#' @importFrom stats optim
#' @export
#' @examples
#' # TODO
#' print("foo")

get_credible_ICC <- function(k, m, alpha = 0.05, p_nodes = 20, rho_nodes = 20) {
  
  # get total AUC
  int_total <- loglike_rho_upper_Gauss(k, m, rho_upper = 1)
  
  # get lower and upper CrIs
  CrI_lower <- optim(0.5, fn = function(rho_upper) {
    tmp <- loglike_rho_upper_Gauss(k = k, m = m, rho_upper = rho_upper, p_nodes = p_nodes, rho_nodes = rho_nodes) - int_total
    abs(exp(tmp) - (alpha / 2))
  }, lower = 0, upper = 1, method = "Brent")$par
  
  CrI_upper <- optim(0.5, fn = function(rho_upper) {
    tmp <- loglike_rho_upper_Gauss(k = k, m = m, rho_upper = rho_upper, p_nodes = p_nodes, rho_nodes = rho_nodes) - int_total
    abs(exp(tmp) - (1 - alpha / 2))
  }, lower = 0, upper = 1, method = "Brent")$par
  
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
#'   use \code{alpha = 0.05} for 95
#' @param p_nodes,rho_nodes number of nodes to use in the Gaussian quadrature
#'   approximate integral method.
#'
#' @importFrom stats optim
#' @export
#' @examples
#' # TODO
#' print("foo")

get_credible_prevalence <- function(k, m, alpha = 0.05, p_nodes = 20, rho_nodes = 20) {
  
  # get total AUC
  int_total <- loglike_p_upper_Gauss(k, m, p_upper = 1)
  
  # get lower and upper CrIs
  CrI_lower <- optim(0.5, fn = function(p_upper) {
    tmp <- loglike_p_upper_Gauss_Vp(k = k, m = m, p_upper = p_upper, p_nodes = p_nodes, rho_nodes = rho_nodes) - int_total
    abs(exp(tmp) - (alpha / 2))
  }, lower = 0, upper = 1, method = "Brent")$par
  
  CrI_upper <- optim(0.5, fn = function(p_upper) {
    tmp <- loglike_p_upper_Gauss_Vp(k = k, m = m, p_upper = p_upper, p_nodes = p_nodes, rho_nodes = rho_nodes) - int_total
    abs(exp(tmp) - (1 - alpha / 2))
  }, lower = 0, upper = 1, method = "Brent")$par
  
  # return
  return(c(lower = CrI_lower, upper = CrI_upper))
}

