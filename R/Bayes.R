
#------------------------------------------------
# reparameterisation of the beta-binomial distribution in terms of a mean (p)
# and an intra-cluster correlation coefficient (rho). Deals with special cases
# that simplify to the binomial or bernoulli distributions
#' @importFrom extraDistr dbbinom
#' @importFrom stats dbinom
#' @noRd

dbbinom_reparam <- function(n, N, p, rho, log_on = TRUE) {
  
  # count number of clusters
  n_clust <- length(n)
  
  if (rho == 0) {
    # simplifies to binomial distribution
    ret <- dbinom(x = n, size = N, prob = p, log = TRUE)
    
  } else if (rho == 1) {
    # likelihood still positive whenever n == 0 or n == N
    ret <- rep(-Inf, n_clust)
    ret[n == 0] <- log(1 - p)
    ret[n == N] <- log(p)
    
  } else {
    if (p == 0) {
      # likelihood 1 whenever n == 0
      ret <- rep(-Inf, n_clust)
      ret[n == 0] <- 0
      
    } else if (p == 1) {
      # likelihood 1 whenever n == N
      ret <- rep(-Inf, n_clust)
      ret[n == N] <- 0
      
    } else {
      # beta-binomial distribution
      alpha <- p * (1 - rho) / rho
      beta <- (1 - p) * (1 - rho) / rho
      ret <- extraDistr::dbbinom(x = n, size = N, alpha = alpha, beta = beta, log = TRUE)
      
    }
  }
  if (!log_on) {
    ret <- exp(ret)
  }
  return(ret)
}

#------------------------------------------------
# draw from reparameterisation of the beta-binomial distribution in terms of a
# mean (p) and an intra-cluster correlation coefficient (rho). Deals with
# special cases that simplify to the binomial distribution
#' @importFrom extraDistr rbbinom
#' @importFrom stats rbinom
#' @noRd

rbbinom_reparam <- function(n_clust, N, p, rho) {
  if (rho == 0) {
    # simplifies to binomial distribution
    ret <- rbinom(n = n_clust, size = N, prob = p)
    
  } else {
    # beta-binomial distribution
    alpha <- p * (1 - rho) / rho
    beta <- (1 - p) * (1 - rho) / rho
    ret <- extraDistr::rbbinom(n = n_clust, size = N, alpha = alpha, beta = beta)
    
  }
  return(ret)
}

#------------------------------------------------
# joint probability of data multiplied by priors on p and rho
#' @importFrom stats dbeta
#' @noRd

loglike_joint <- function(n, N, p, rho,
                          prior_p_shape1 = 1, prior_p_shape2 = 1,
                          prior_rho_shape1 = 1, prior_rho_shape2 = 1) {
  
  dbeta(p, shape1 = prior_p_shape1, shape2 = prior_p_shape2, log = TRUE) +
    dbeta(rho, shape1 = prior_rho_shape1, shape2 = prior_rho_shape2, log = TRUE) +
    sum(dbbinom_reparam(n = n, N = N, p = p, rho = rho, log_on = TRUE))
}

#------------------------------------------------
# marginal log-likelihood of rho, integrated over p by adaptive quadrature
#' @noRd

loglike_marginal_rho <- function(n, N, rho, n_intervals = 40,
                                 prior_p_shape1 = 1, prior_p_shape2 = 1,
                                 prior_rho_shape1 = 1, prior_rho_shape2 = 1,
                                 debug_on = FALSE) {
  
  
  # integrate over rho via adaptive quadrature
  df_quad <- adaptive_quadrature(f1 = function(p) {
    loglike_joint(n = n, N = N, p = p, rho = rho, prior_p_shape1 = prior_p_shape1,
                  prior_p_shape2 = prior_p_shape2, prior_rho_shape1 = prior_rho_shape1,
                  prior_rho_shape2 = prior_rho_shape2)
  }, n_intervals = n_intervals, left = 0, right = 1, debug_on = debug_on)
  
  # sum area over intervals
  log_area <- df_quad$log_area_S
  if (all(log_area == -Inf)) {
    log_area_sum <- -Inf
  } else {
    log_area_sum <- max(log_area) + log(sum(exp(log_area - max(log_area))))
  }
  
  return(log_area_sum)
}

#------------------------------------------------
# marginal log-likelihood of p, integrated over rho by adaptive quadrature
#' @noRd

loglike_marginal_p <- function(n, N, p, n_intervals = 40,
                               prior_p_shape1 = 1, prior_p_shape2 = 1,
                               prior_rho_shape1 = 1, prior_rho_shape2 = 1,
                               debug_on = FALSE) {
  
  # integrate over rho via adaptive quadrature
  df_quad <- adaptive_quadrature(f1 = function(rho) {
    loglike_joint(n = n, N = N, p = p, rho = rho, prior_p_shape1 = prior_p_shape1,
                  prior_p_shape2 = prior_p_shape2, prior_rho_shape1 = prior_rho_shape1,
                  prior_rho_shape2 = prior_rho_shape2)
  }, n_intervals = n_intervals, left = 0, right = 1, debug_on = debug_on)
  
  # sum area over intervals
  log_area <- df_quad$log_area_S
  if (all(log_area == -Inf)) {
    log_area_sum <- -Inf
  } else {
    log_area_sum <- max(log_area) + log(sum(exp(log_area - max(log_area))))
  }
  
  return(log_area_sum)
}

#------------------------------------------------
##' @name get_posterior
##' @rdname get_posterior
##'
#' @title Estimate prevalence and intra-cluster correlation from raw counts
#'
#' @description Takes raw counts of the number of positive samples per cluster
#'   (numerator) and the number of tested samples per cluster (denominator) and
#'   returns posterior estimates of the prevalence and the intra-cluster
#'   correlation coefficient (ICC).
#'
#' @details There are two unknown quantities in the DRpower model: the
#'   prevalence and the ICC. Thes following functions integrate over a prior on
#'   one quantity to give the marginal posterior distribution of the other.
#'   Possible outputs include the posterior mean, median, credible interval
#'   (CrI), probability of being above a threshold, and the full posterior
#'   distribution. For speed, distributions are approximated using an adaptive
#'   quadrature approach, in which the full distribution is split into intervals
#'   of differing width and each intervals is approximated using Simpson's rule.
#'   The number of intervals used in quadrature can be increased for more
#'   accurate results, at the cost of slower speed.
#'
#' @param n,N the numerator (\code{n}) and denominator (\code{N}) per cluster
#'   (vectors).
#' @param alpha the significance level of the credible interval - for example,
#'   use \code{alpha = 0.05} for a 95\% interval. See also \code{CrI_type}
#'   argument for how this is calculated.
#' @param prev_thresh return the probability that the prevalence is above this
#'   threshold. Can be a vector, in which case the return object contains one
#'   value for each input.
#' @param prior_prev_shape1,prior_prev_shape2,prior_ICC_shape1,prior_ICC_shape2
#'   parameters that dictate the shape of the Beta priors on prevalence and the
#'   ICC. Increasing the first shape parameter (e.g. \code{prior_prev_shape1})
#'   pushes the distribution towards 1, increasing the second shape parameter
#'   (e.g. \code{prior_prev_shape2}) pushes the distribution towards 0.
#'   Increasing both shape parameters squeezes the distribution towards the
#'   centre and therefore makes it narrower.
#' @param return_type a list object specifying which outputs to produce (only a
#'   subset of outputs is active by default to reduce computation time and to
#'   simplify things for the basic user). The following options are available:
#'   \itemize{
#'     \item \code{mean_on}: if \code{TRUE} then return the posterior mean.
#'     \item \code{median_on}: if \code{TRUE} then return the posterior median.
#'     \item \code{CrI_on}: if \code{TRUE} then return the posterior credible
#'     interval at significance level \code{alpha}. See \code{CrI_type} argument
#'     for how this is calculated.
#'     \item \code{thresh_on}: if \code{TRUE} then return the posterior
#'     probability of being above the threshold specified by \code{prev_thresh}.
#'     \item \code{full_on}: if \code{TRUE} then return the full posterior
#'     distribution (as approximated using the adaptive quadrature approach) in
#'     0.1\% intervals from 0\% to 100\%.
#'   }
#' @param CrI_type which method to use when computing credible intervals, with
#'   options "ETI" (equal-tailed interval) and "HDI" (high-density interval).
#'   The ETI searches a distance \code{alpha/2} from either side of the [0,1]
#'   interval. The HDI method returns the narrowest interval that subtends a
#'   proportion \code{1-alpha} of the distribution. The HDI method is used by
#'   default.
#' @param n_intervals the number of intervals used in the adaptive quadrature
#'   method. Increasing this value gives a more accurate representation of the
#'   true posterior, but comes at the cost of reduced speed.
#' @param debug_on for use in debugging. If \code{TRUE} and if \code{use_cpp =
#'   FALSE} then produces a plot of the posterior distribution evaluated by
#'   brute force (black) overlaid with the adaptive quadrature approximation
#'   (blue) for direct comparison. If \code{use_cpp = TRUE} then only the
#'   approximation is plotted. If the approximate distribution does not agree
#'   closely with the brute force solution then consider increasing the value of
#'   \code{n_intervals}. Note that producing this plot can be very slow, and so
#'   this option should be turned off when not needed.
#' @param use_cpp if \code{TRUE} then use an Rcpp implementation of the adaptive
#'   quadrature approach that is much faster and therefore useful when running
#'   the function a large number of times.
NULL

##' @rdname get_posterior
#' @importFrom stats optim
#' @importFrom graphics lines points abline
#' @export

get_prevalence <- function(n, N, alpha = 0.05, prev_thresh = 0.05,
                           prior_prev_shape1 = 1.0, prior_prev_shape2 = 1.0,
                           prior_ICC_shape1 = 1.0, prior_ICC_shape2 = 3.0,
                           return_type = list(mean_on = FALSE,
                                              median_on = TRUE,
                                              CrI_on = TRUE,
                                              thresh_on = TRUE,
                                              full_on = FALSE),
                           CrI_type = "HDI", n_intervals = 20,
                           debug_on = FALSE, use_cpp = TRUE) {
  
  # avoid "no visible binding" note
  dummy <- NULL
  
  # check inputs
  assert_vector_pos_int(n)
  assert_vector_pos_int(N)
  if (length(N) == 1) {
    N <- rep(N, length(n))
  }
  assert_same_length(n, N, message = "N must be either a single value (all clusters the same size) or a vector with the same length as n")
  assert_single_bounded(alpha, inclusive_left = FALSE, inclusive_right = FALSE)
  assert_vector_bounded(prev_thresh)
  assert_single_bounded(prior_prev_shape1, left = 1, right = 1e3)
  assert_single_bounded(prior_prev_shape2, left = 1, right = 1e3)
  assert_single_bounded(prior_ICC_shape1, left = 1, right = 1e3)
  assert_single_bounded(prior_ICC_shape2, left = 1, right = 1e3)
  assert_list(return_type, message = "return_type must be a list")
  assert_in(c("mean_on", "median_on", "CrI_on", "thresh_on", "full_on"), names(return_type), name_y = "names(return_type)")
  assert_single_string(CrI_type)
  assert_in(CrI_type, c("HDI", "ETI"))
  assert_pos_int(n_intervals)
  assert_greq(n_intervals, 5)
  assert_single_logical(debug_on)
  assert_single_logical(use_cpp)
  
  # split based on C++ vs. R
  if (use_cpp) {
    
    # get arguments into list
    args_params <- list(n = n,
                        N = N,
                        prior_prev_shape1 = prior_prev_shape1,
                        prior_prev_shape2 = prior_prev_shape2,
                        prior_ICC_shape1 = prior_ICC_shape1,
                        prior_ICC_shape2 = prior_ICC_shape2,
                        n_intervals = n_intervals)
    
    # run efficient C++ function
    output_raw <- get_prevalence_cpp(args_params)
    df_quad <- as.data.frame(output_raw)
    
    # debug distribution
    if (debug_on) {
      
      # normalise and add coefficients
      df_norm <- normalise_quadrature(df_quad)
      
      # calculate interpolated curve from Simpson's rule
      x <- seq(0, 1, l = 201)
      z <- findInterval(x, vec = df_norm$x0)
      fx <- df_norm$A[z]*x^2 + df_norm$B[z]*x + df_norm$C[z]
      
      # produce plot
      plot(x, fx, type = 'l')
      points(df_norm$x0, exp(df_norm$log_y0), pch = 20, cex = 0.75)
      abline(h = 0, lty = 3)
    }
    
  } else {
    
    # get distribution of p via adaptive quadrature
    df_quad <- adaptive_quadrature(f1 = function(p) {
      loglike_marginal_p(n = n, N = N, p = p, n_intervals = n_intervals,
                         prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                         prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
    }, n_intervals = n_intervals, left = 0, right = 1, debug_on = debug_on)
    
  }
  
  # normalise and add coefficients
  df_norm <- normalise_quadrature(df_quad)
  
  # initialise return data.frame
  ret <- data.frame(dummy = NA)
  
  # solve for posterior mean
  if (return_type$mean_on) {
    ret$mean = NA
  }
  
  # solve for posterior median
  if (return_type$median_on) {
    ret$median <- qquad(df_norm, q = 0.5)
  }
  
  # solve for lower and upper CrIs
  if (return_type$CrI_on) {
    if (CrI_type == "ETI") {
      ret$lower <- qquad(df_norm, q = alpha / 2)
      ret$upper <- qquad(df_norm, q = 1 - alpha / 2)
    } else if (CrI_type == "HDI") {
      HDI <- get_HDI(df_norm, alpha = alpha)
      ret$lower <- HDI["lower"]
      ret$upper <- HDI["upper"]
    }
  }
  
  # solve for prob above threshold
  if (return_type$thresh_on) {
    prob_above_threshold <- rep(NA, length(prev_thresh))
    for (i in seq_along(prev_thresh)) {
      prob_above_threshold[i] <- 1 - pquad(df_norm, p = prev_thresh[i])
    }
    prob_above_threshold[prob_above_threshold > 1] <- 1
    if (length(prev_thresh) == 1) {
      ret$prob_above_threshold <- prob_above_threshold
    } else {
      ret$prob_above_threshold <- I(list(prob_above_threshold))
    }
  }
  
  # get full posterior interpolated curve from Simpson's rule
  if (return_type$full_on) {
    x <- seq(0, 1, l = 1001)
    z <- findInterval(x, vec = df_norm$x0)
    y <- df_norm$A[z]*x^2 + df_norm$B[z]*x + df_norm$C[z]
    ret$full <- I(list(y))
  }
  
  # finalise output
  ret <- dplyr::select(ret, -dummy)
  row.names(ret) <- NULL
  
  # return
  return(ret)
}

##' @rdname get_posterior
#' @importFrom stats optim qbeta
#' @importFrom graphics lines points abline
#' @export

get_ICC <- function(n, N, alpha = 0.05,
                    prior_prev_shape1 = 1.0, prior_prev_shape2 = 1.0,
                    prior_ICC_shape1 = 1.0, prior_ICC_shape2 = 1.0,
                    return_type = list(mean_on = FALSE,
                                       median_on = TRUE,
                                       CrI_on = TRUE,
                                       full_on = FALSE),
                    CrI_type = "HDI", n_intervals = 20,
                    debug_on = FALSE, use_cpp = TRUE) {
  
  # avoid "no visible binding" note
  dummy <- NULL
  
  # check inputs
  assert_vector_pos_int(n)
  assert_vector_pos_int(N)
  if (length(N) == 1) {
    N <- rep(N, length(n))
  }
  assert_same_length(n, N, message = "N must be either a single value (all clusters the same size) or a vector with the same length as n")
  assert_single_bounded(prior_prev_shape1, left = 1, right = 1e3)
  assert_single_bounded(prior_prev_shape2, left = 1, right = 1e3)
  assert_single_bounded(prior_ICC_shape1, left = 1, right = 1e3)
  assert_single_bounded(prior_ICC_shape2, left = 1, right = 1e3)
  assert_list(return_type, message = "return_type must be a list")
  assert_in(c("mean_on", "median_on", "CrI_on", "full_on"), names(return_type), name_y = "names(return_type)")
  assert_single_string(CrI_type)
  assert_in(CrI_type, c("HDI", "ETI"))
  assert_pos_int(n_intervals)
  assert_greq(n_intervals, 5)
  assert_single_logical(debug_on)
  assert_single_logical(use_cpp)
  
  # split based on C++ vs. R
  if (use_cpp) {
    
    # get arguments into list
    args_params <- list(n = n,
                        N = N,
                        prior_prev_shape1 = prior_prev_shape1,
                        prior_prev_shape2 = prior_prev_shape2,
                        prior_ICC_shape1 = prior_ICC_shape1,
                        prior_ICC_shape2 = prior_ICC_shape2,
                        n_intervals = n_intervals)
    
    # run efficient C++ function
    output_raw <- get_ICC_cpp(args_params)
    df_quad <- as.data.frame(output_raw)
    
    # debug distribution
    if (debug_on) {
      
      # normalise and add coefficients
      df_norm <- normalise_quadrature(df_quad)
      
      # calculate interpolated curve from Simpson's rule
      x <- seq(0, 1, l = 201)
      z <- findInterval(x, vec = df_norm$x0)
      fx <- df_norm$A[z]*x^2 + df_norm$B[z]*x + df_norm$C[z]
      
      # produce plot
      plot(x, fx, type = 'l')
      points(df_norm$x0, exp(df_norm$log_y0), pch = 20, cex = 0.75)
      abline(h = 0, lty = 3)
    }
    
  } else {
    
    # get distribution of rho via adaptive quadrature
    df_quad <- adaptive_quadrature(f1 = function(rho) {
      loglike_marginal_rho(n = n, N = N, rho = rho, n_intervals = n_intervals,
                           prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                           prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
    }, n_intervals = n_intervals, left = 0, right = 1, debug_on = debug_on)
    
  }
  
  # normalise and add coefficients
  df_norm <- normalise_quadrature(df_quad)
  
  # initialise return data.frame
  ret <- data.frame(dummy = NA)
  
  # solve for posterior mean
  if (return_type$mean_on) {
    ret$mean = NA
  }
  
  # solve for posterior median
  if (return_type$median_on) {
    ret$median <- qquad(df_norm, q = 0.5)
  }
  
  # solve for lower and upper CrIs
  if (return_type$CrI_on) {
    if (CrI_type == "ETI") {
      ret$lower <- qquad(df_norm, q = alpha / 2)
      ret$upper <- qquad(df_norm, q = 1 - alpha / 2)
    } else if (CrI_type == "HDI") {
      HDI <- get_HDI(df_norm, alpha = alpha)
      ret$lower <- HDI["lower"]
      ret$upper <- HDI["upper"]
    }
  }
  
  # get full posterior interpolated curve from Simpson's rule
  if (return_type$full_on) {
    x <- seq(0, 1, l = 1001)
    z <- findInterval(x, vec = df_norm$x0)
    y <- df_norm$A[z]*x^2 + df_norm$B[z]*x + df_norm$C[z]
    ret$full <- I(list(y))
  }
  
  # finalise output
  ret <- dplyr::select(ret, -dummy)
  row.names(ret) <- NULL
  
  # return
  return(ret)
}

#------------------------------------------------
#' @title Get posterior distribution of both prevalence and the ICC on a grid
#'
#' @description Get posterior distribution of both prevalence and the ICC on a
#'   grid.
#'
#' @param n,N the numerator (\code{n}) and denominator (\code{N}) per cluster.
#' @param prior_prev_shape1,prior_prev_shape2,prior_ICC_shape1,prior_ICC_shape2
#'   parameters that dictate the shape of the priors on prevalence and the ICC.
#'   Increasing the first shape parameter (e.g. \code{prior_p_shape1}) pushes
#'   the distribution towards 1, increasing the second shape parameter (e.g.
#'   \code{prior_p_shape2}) pushes the distribution towards 0. Increasing both
#'   shape parameters squeezes the distribution and makes it narrower.
#' @param prev_cells,ICC_cells the number of cells in the grid in each
#'   dimension.
#'  
#' @export

get_joint_grid <- function(n, N, 
                           prior_prev_shape1 = 1.0, prior_prev_shape2 = 1.0,
                           prior_ICC_shape1 = 1.0, prior_ICC_shape2 = 1.0,
                           prev_cells = 64, ICC_cells = 64) {
  
  # check inputs
  assert_vector_pos_int(n)
  assert_vector_pos_int(N)
  if (length(N) == 1) {
    N <- rep(N, length(n))
  }
  assert_single_bounded(prior_prev_shape1, left = 1, right = 1e3)
  assert_single_bounded(prior_prev_shape2, left = 1, right = 1e3)
  assert_single_bounded(prior_ICC_shape1, left = 1, right = 1e3)
  assert_single_bounded(prior_ICC_shape2, left = 1, right = 1e3)
  assert_single_pos_int(prev_cells, zero_allowed = FALSE)
  assert_single_pos_int(ICC_cells, zero_allowed = FALSE)
  
  # make parameter grids
  p_vec <- seq(0, 1, l = prev_cells)
  rho_vec <- seq(0, 1, l = ICC_cells)
  p <- matrix(rep(p_vec, each = length(rho_vec)), length(rho_vec))
  rho <- matrix(rep(rho_vec, length(p_vec)), length(rho_vec))
  
  # define implied beta-binomial parameters
  alpha <- p*(1/rho - 1)
  beta <- (1 - p)*(1/rho - 1)
  
  # evaluate beta-binomial likelihood
  ll <- 0
  for (i in seq_along(n)) {
    ll <- ll + extraDistr::dbbinom(x = n[i], size = N[i], alpha = alpha, beta = beta, log = TRUE)
  }
  ll <- matrix(ll, nrow = length(rho_vec))
  
  # special cases
  w <- which(rho == 0)
  ll[w] <- 0
  for (i in seq_along(n)) {
    ll[w] <- ll[w] + dbinom(x = n[i], size = N[i], prob = p[w], log = TRUE)
  }
  w <- which(rho == 1)
  ll[w] <- 0
  for (i in seq_along(n)) {
    if (n[i] == 0) {
      ll[w] <- ll[w] + log(1 - p[w])
    } else if (n[i] == N[i]) {
      ll[w] <- ll[w] + log(p[w])
    } else {
      ll[w] <- -Inf
    }
  }
  w <- which(p == 0)
  if (all(n == 0)) {
    ll[w] <- 0
  } else {
    ll[w] <- -Inf
  }
  w <- which(p == 1)
  if (all(n == N)) {
    ll[w] <- 0
  } else {
    ll[w] <- -Inf
  }
  
  # normalise to largest value and exponentiate
  ret <- exp(ll - max(ll))
  ret <- ret / sum(ret)
  
  # return
  return(ret)
}
