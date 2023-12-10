#------------------------------------------------
#' @title Check that DRpower package has loaded successfully
#'
#' @description Simple function to check that DRpower package has loaded 
#'   successfully. Prints "DRpower loaded successfully!" if so.
#'
#' @export

check_DRpower_loaded <- function() {
  message("DRpower version 1.0.2 loaded successfully!")
}

#------------------------------------------------
# reparameterisation of the beta-binomial distribution in terms of a mean (p)
# and an intra-cluster correlation coefficient (rho). The shape parameters of
# this distribution are alpha = p*(1/rho - 1) and beta = (1 - p)*(1/rho - 1).
# The implied beta distribution has mean p and variance p*(1-p)*rho. Deals with
# special cases that simplify to the binomial or bernoulli distributions
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
    # perfect correlation within clusters. Likelihood finite whenever n == 0 or
    # n == N
    ret <- rep(-Inf, n_clust)
    ret[n == 0] <- log(1 - p)
    ret[n == N] <- log(p)
    
  } else {
    if (p == 0) {
      # no positives allowed irrespective of ICC. Likelihood 1 whenever n == 0
      ret <- rep(-Inf, n_clust)
      ret[n == 0] <- 0
      
    } else if (p == 1) {
      # no negatives allowed irrespective of ICC. Likelihood 1 whenever n == N
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
# draw from reparameterisation of the beta-binomial distribution (see
# dbbinom_reparam()). Deals with special cases that simplify to the binomial or
# bernoulli distributions
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
# joint probability of data and rho, integrated over the prior on p by adaptive
# quadrature
#' @noRd

loglike_joint_rho <- function(n, N, rho, n_intervals = 40,
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
# joint probability of data and p, integrated over the prior on rho by adaptive
# quadrature
#' @noRd

loglike_joint_p <- function(n, N, p, rho_fixed = NULL, n_intervals = 40,
                            prior_p_shape1 = 1, prior_p_shape2 = 1,
                            prior_rho_shape1 = 1, prior_rho_shape2 = 1,
                            debug_on = FALSE) {
  
  # simple evaluation if using fixed rho (no need to integrate)
  if (!is.null(rho_fixed)) {
    ret <- loglike_joint(n = n, N = N, p = p, rho = rho_fixed, prior_p_shape1 = prior_p_shape1,
                         prior_p_shape2 = prior_p_shape2, prior_rho_shape1 = prior_rho_shape1,
                         prior_rho_shape2 = prior_rho_shape2)
    return(ret)
  }
  
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
#'   returns posterior estimates of the prevalence and intra-cluster correlation
#'   coefficient (ICC).
#'
#' @details There are two unknown quantities in the DRpower model: the
#'   prevalence and the intra-cluster correlation (ICC). These functions
#'   integrate over a prior on one quantity to arrive at the marginal posterior
#'   distribution of the other. Possible outputs include the maximum \emph{a
#'   posteriori} (MAP) estimate, the posterior mean, posterior median, credible
#'   interval (CrI), probability of being above a set threshold, and the full
#'   posterior distribution. For speed, distributions are approximated using an
#'   adaptive quadrature approach in which the full distribution is split into
#'   intervals and each intervals is approximated using Simpson's rule. The
#'   number of intervals used in quadrature can be increased for more accurate
#'   results at the cost of slower speed.
#'
#' @param n,N the numerator (\code{n}) and denominator (\code{N}) per cluster.
#'   These are both integer vectors.
#' @param alpha the significance level of the credible interval - for example,
#'   use \code{alpha = 0.05} for a 95\% interval. See also \code{CrI_type}
#'   argument for how this is calculated.
#' @param prev_thresh the prevalence threshold that we are comparing against.
#'   Can be a vector, in which case the return object contains one value for
#'   each input.
#' @param ICC normally this should be set to \code{NULL} (the default), in which
#'   case the ICC is estimated from the data. However, a fixed value can be
#'   entered here, in which case this overrides the use of the prior
#'   distribution as specified by \code{prior_ICC_shape1} and
#'   \code{prior_ICC_shape2}.
#' @param prior_prev_shape1,prior_prev_shape2,prior_ICC_shape1,prior_ICC_shape2
#'   parameters that dictate the shape of the Beta priors on prevalence and the
#'   ICC. See the \href{https://en.wikipedia.org/wiki/Beta_distribution}{Wikipedia page on the Beta
#'   distribution} for more detail. The default values of these parameters were chosen
#'   based on an
#'   \href{https://mrc-ide.github.io/DRpower/articles/historical_analysis.html}{analysis
#'   of historical pfhrp2/3 studies}, although this does not guarantee that they
#'   will be suitable in all settings.
#' @param MAP_on,post_mean_on,post_median_on,post_CrI_on,post_thresh_on,post_full_on a
#'   series of boolean values specifying which outputs to produce. The options are:
#'   \itemize{
#'     \item \code{MAP_on}: if \code{TRUE} then return the maximum \emph{a
#'     posteriori}.
#'     \item \code{post_mean_on}: if \code{TRUE} then return the posterior mean.
#'     \item \code{post_median_on}: if \code{TRUE} then return the posterior
#'     median.
#'     \item \code{post_CrI_on}: if \code{TRUE} then return the posterior
#'     credible interval at significance level \code{alpha}. See \code{CrI_type}
#'     argument for how this is calculated.
#'     \item \code{post_thresh_on}: if \code{TRUE} then return the posterior
#'     probability of being above the threshold(s) specified by
#'     \code{prev_thresh}.
#'     \item \code{post_full_on}: if \code{TRUE} then return the full posterior
#'     distribution, produced using the adaptive quadrature approach, at breaks
#'     specified by \code{post_full_breaks}.
#'   }
#' @param post_full_breaks a vector of breaks at which to evaluate the full
#'   posterior distribution (only if \code{post_full_on = TRUE}). Defaults to
#'   0.1\% intervals from 0\% to 100\%.
#' @param CrI_type which method to use when computing credible intervals.
#'   Options are "ETI" (equal-tailed interval) or "HDI" (high-density interval).
#'   The ETI searches a distance \code{alpha/2} from either side of the [0,1]
#'   interval. The HDI method returns the narrowest interval that subtends a
#'   proportion \code{1-alpha} of the distribution. The HDI method is used by
#'   default as it guarantees that the MAP estimate is within the credible
#'   interval, which is not always the case for the ETI.
#' @param n_intervals the number of intervals used in the adaptive quadrature
#'   method. Increasing this value gives a more accurate representation of the
#'   true posterior, but comes at the cost of reduced speed.
#' @param round_digits the number of digits after the decimal point that are
#'   used when reporting estimates. This is to simplify results and to avoid
#'   giving the false impression of extreme precision.
#' @param use_cpp if \code{TRUE} (the default) then use an Rcpp implementation
#'   of the adaptive quadrature approach that is much faster than the base R
#'   method.
#' @param silent if \code{TRUE} then suppress all console output.
NULL

##' @rdname get_posterior
#' @importFrom stats optim
#' @importFrom graphics lines points abline
#' @examples
#' # basic example of estimating prevalence and
#' # ICC from observed counts
#' sample_size <- c(80, 110, 120)
#' deletions <- c(3, 5, 6)
#' 
#' get_prevalence(n = deletions, N = sample_size)
#' get_ICC(n = deletions, N = sample_size)
#' 
#' @export

get_prevalence <- function(n, N, alpha = 0.05, prev_thresh = 0.05, ICC = NULL,
                           prior_prev_shape1 = 1.0, prior_prev_shape2 = 1.0,
                           prior_ICC_shape1 = 1.0, prior_ICC_shape2 = 9.0,
                           MAP_on = TRUE, post_mean_on = FALSE, post_median_on = FALSE,
                           post_CrI_on = TRUE, post_thresh_on = TRUE,
                           post_full_on = FALSE, post_full_breaks = seq(0, 1, l = 1001),
                           CrI_type = "HDI", n_intervals = 20, round_digits = 2,
                           use_cpp = TRUE, silent = FALSE) {
  
  # avoid "no visible binding" note
  dummy <- A <- B <- C <- x0 <- x1 <- post_mean <- NULL
  
  # check inputs
  assert_vector_pos_int(n)
  assert_vector_pos_int(N)
  if (length(N) == 1) {
    N <- rep(N, length(n))
  }
  assert_same_length(n, N, message = "N must be either a single value (all clusters the same size) or a vector with the same length as n")
  assert_single_bounded(alpha, inclusive_left = FALSE, inclusive_right = FALSE)
  assert_vector_bounded(prev_thresh)
  if (!is.null(ICC)) {
    if (!silent) {
      message(sprintf("Note: the ICC is set to the fixed value %s. This overrides the use of the prior on ICC. Use silent = TRUE to suppress this message.", ICC))
    }
    assert_single_bounded(ICC)
  }
  assert_single_bounded(prior_prev_shape1, left = 1e-3, right = 1e3)
  assert_single_bounded(prior_prev_shape2, left = 1e-3, right = 1e3)
  assert_single_bounded(prior_ICC_shape1, left = 1e-3, right = 1e3)
  assert_single_bounded(prior_ICC_shape2, left = 1e-3, right = 1e3)
  assert_single_logical(MAP_on)
  assert_single_logical(post_mean_on)
  assert_single_logical(post_median_on)
  assert_single_logical(post_CrI_on)
  assert_single_logical(post_thresh_on)
  assert_single_logical(post_full_on)
  assert_vector_bounded(post_full_breaks)
  assert_increasing(post_full_breaks)
  assert_single_string(CrI_type)
  assert_in(CrI_type, c("HDI", "ETI"))
  assert_pos_int(n_intervals)
  assert_greq(n_intervals, 5)
  assert_single_pos_int(round_digits, zero_allowed = FALSE)
  assert_single_logical(use_cpp)
  assert_single_logical(silent)
  
  # split based on C++ vs. R
  if (use_cpp) {
    
    # get arguments into list
    args_params <- list(n = n,
                        N = N,
                        ICC = ifelse(is.null(ICC), -1, ICC),
                        prior_prev_shape1 = prior_prev_shape1,
                        prior_prev_shape2 = prior_prev_shape2,
                        prior_ICC_shape1 = prior_ICC_shape1,
                        prior_ICC_shape2 = prior_ICC_shape2,
                        n_intervals = n_intervals)
    
    # run efficient C++ function
    output_raw <- get_prevalence_cpp(args_params)
    df_quad <- as.data.frame(output_raw)
    
    # debug distribution
    if (FALSE) {
      
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
    # NB, set debug_on = TRUE to plot distribution
    df_quad <- adaptive_quadrature(f1 = function(p) {
      loglike_joint_p(n = n, N = N, p = p, rho_fixed = ICC, n_intervals = n_intervals,
                      prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                      prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
    }, n_intervals = n_intervals, left = 0, right = 1, debug_on = FALSE)
    
  }
  
  # normalise and add coefficients
  df_norm <- normalise_quadrature(df_quad)
  
  # initialise return data.frame
  ret <- data.frame(dummy = NA)
  
  # solve for maximum a posteriori (MAP)
  if (MAP_on) {
    MAP <- get_max_x(df_norm)
    ret$MAP <- round(MAP * 100, round_digits)
  }
  
  # solve for posterior mean
  if (post_mean_on) {
    post_mean = df_norm %>%
      dplyr::mutate(post_mean = 1/4*A*(x1^4 - x0^4) + 1/3*B*(x1^3 - x0^3) + 1/2*C*(x1^2 - x0^2)) %>%
      dplyr::pull(post_mean) %>%
      sum()
    ret$post_mean = round(post_mean * 100, round_digits)
  }
  
  # solve for posterior median
  if (post_median_on) {
    post_median <- qquad(df_norm, q = 0.5)
    ret$post_median <- round(post_median * 100,  round_digits)
  }
  
  # solve for lower and upper CrIs
  if (post_CrI_on) {
    if (CrI_type == "ETI") {
      CrI_lower <- qquad(df_norm, q = alpha / 2)
      CrI_upper <- qquad(df_norm, q = 1 - alpha / 2)
    } else if (CrI_type == "HDI") {
      HDI <- get_HDI(df_norm, alpha = alpha)
      CrI_lower <- HDI["lower"]
      CrI_upper <- HDI["upper"]
    }
    ret$CrI_lower <- round(CrI_lower * 100, round_digits)
    ret$CrI_upper <- round(CrI_upper * 100, round_digits)
  }
  
  # solve for prob above threshold
  if (post_thresh_on) {
    prob_above_threshold <- rep(NA, length(prev_thresh))
    for (i in seq_along(prev_thresh)) {
      prob_above_threshold[i] <- 1 - pquad(df_norm, p = prev_thresh[i])
    }
    prob_above_threshold[prob_above_threshold > 1] <- 1
    prob_above_threshold <- round(prob_above_threshold, 4)
    if (length(prev_thresh) == 1) {
      ret$prob_above_threshold <- prob_above_threshold
    } else {
      ret$prob_above_threshold <- I(list(prob_above_threshold))
    }
  }
  
  # get full posterior interpolated curve from Simpson's rule
  if (post_full_on) {
    x <- post_full_breaks
    z <- findInterval(x, vec = df_norm$x0)
    y <- df_norm$A[z]*x^2 + df_norm$B[z]*x + df_norm$C[z]
    y[y < 0] <- 0
    ret$post_full <- I(list(y))
  }
  
  # finalise output
  ret <- dplyr::select(ret, -dummy)
  row.names(ret) <- NULL
  
  # return
  return(ret)
}

# NB, this form of documentation means that both get_prevalence() and get_ICC()
# fall within the same single help page
##' @rdname get_posterior
#' @importFrom stats optim qbeta
#' @importFrom graphics lines points abline
#' @export

get_ICC <- function(n, N, alpha = 0.05,
                    prior_prev_shape1 = 1.0, prior_prev_shape2 = 1.0,
                    prior_ICC_shape1 = 1.0, prior_ICC_shape2 = 9.0,
                    MAP_on = TRUE, post_mean_on = FALSE, post_median_on = FALSE,
                    post_CrI_on = TRUE, post_full_on = FALSE,
                    post_full_breaks = seq(0, 1, l = 1001), CrI_type = "HDI",
                    n_intervals = 20, round_digits = 4, use_cpp = TRUE) {
  
  # avoid "no visible binding" note
  dummy <- A <- B <- C <- x0 <- x1 <- post_mean <- NULL
  
  # check inputs
  assert_vector_pos_int(n)
  assert_vector_pos_int(N)
  if (length(N) == 1) {
    N <- rep(N, length(n))
  }
  assert_same_length(n, N, message = "N must be either a single value (all clusters the same size) or a vector with the same length as n")
  assert_single_bounded(prior_prev_shape1, left = 1e-3, right = 1e3)
  assert_single_bounded(prior_prev_shape2, left = 1e-3, right = 1e3)
  assert_single_bounded(prior_ICC_shape1, left = 1e-3, right = 1e3)
  assert_single_bounded(prior_ICC_shape2, left = 1e-3, right = 1e3)
  assert_single_logical(MAP_on)
  assert_single_logical(post_mean_on)
  assert_single_logical(post_median_on)
  assert_single_logical(post_CrI_on)
  assert_single_logical(post_full_on)
  assert_vector_bounded(post_full_breaks)
  assert_increasing(post_full_breaks)
  assert_single_string(CrI_type)
  assert_in(CrI_type, c("HDI", "ETI"))
  assert_pos_int(n_intervals)
  assert_greq(n_intervals, 5)
  assert_single_pos_int(round_digits, zero_allowed = FALSE)
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
    if (FALSE) {
      
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
    # NB, set debug_on = TRUE to plot distribution
    df_quad <- adaptive_quadrature(f1 = function(rho) {
      loglike_joint_rho(n = n, N = N, rho = rho, n_intervals = n_intervals,
                        prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                        prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
    }, n_intervals = n_intervals, left = 0, right = 1, debug_on = FALSE)
    
  }
  
  # normalise and add coefficients
  df_norm <- normalise_quadrature(df_quad)
  
  # initialise return data.frame
  ret <- data.frame(dummy = NA)
  
  # solve for maximum a posteriori (MAP)
  if (MAP_on) {
    MAP <- get_max_x(df_norm)
    ret$MAP <- round(MAP, round_digits)
  }
  
  # solve for posterior mean
  if (post_mean_on) {
    post_mean = df_norm %>%
      dplyr::mutate(post_mean = 1/4*A*(x1^4 - x0^4) + 1/3*B*(x1^3 - x0^3) + 1/2*C*(x1^2 - x0^2)) %>%
      dplyr::pull(post_mean) %>%
      sum()
    ret$post_mean = round(post_mean, round_digits)
  }
  
  # solve for posterior median
  if (post_median_on) {
    post_median <- qquad(df_norm, q = 0.5)
    ret$post_median <- round(post_median, round_digits)
  }
  
  # solve for lower and upper CrIs
  if (post_CrI_on) {
    if (CrI_type == "ETI") {
      CrI_lower <- qquad(df_norm, q = alpha / 2)
      CrI_upper <- qquad(df_norm, q = 1 - alpha / 2)
    } else if (CrI_type == "HDI") {
      HDI <- get_HDI(df_norm, alpha = alpha)
      CrI_lower <- HDI["lower"]
      CrI_upper <- HDI["upper"]
    }
    ret$CrI_lower <- round(CrI_lower, round_digits)
    ret$CrI_upper <- round(CrI_upper, round_digits)
  }
  
  # get full posterior interpolated curve from Simpson's rule
  if (post_full_on) {
    x <- post_full_breaks
    z <- findInterval(x, vec = df_norm$x0)
    y <- df_norm$A[z]*x^2 + df_norm$B[z]*x + df_norm$C[z]
    y[y < 0] <- 0
    ret$post_full <- I(list(y))
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
#'   grid. Prevalence is returned in columns, ICC in rows. See also
#'   \code{plot_joint} for how to make a contour plot from this grid.
#'
#' @inheritParams get_posterior
#' @param prev_breaks,ICC_breaks the values at which to evaluate the posterior
#'   in both dimensions. Prevalence is returned in columns, ICC in rows.
#'
#' @examples
#' get_joint(n = c(5, 2, 9), N = c(100, 80, 120))
#' 
#' @export

get_joint <- function(n, N, 
                      prior_prev_shape1 = 1.0, prior_prev_shape2 = 1.0,
                      prior_ICC_shape1 = 1.0, prior_ICC_shape2 = 9.0,
                      prev_breaks = seq(0, 1, 0.01), ICC_breaks = seq(0, 1, 0.01)) {
  
  # check inputs
  assert_vector_pos_int(n)
  assert_vector_pos_int(N)
  if (length(N) == 1) {
    N <- rep(N, length(n))
  }
  assert_single_bounded(prior_prev_shape1, left = 1e-3, right = 1e3)
  assert_single_bounded(prior_prev_shape2, left = 1e-3, right = 1e3)
  assert_single_bounded(prior_ICC_shape1, left = 1e-3, right = 1e3)
  assert_single_bounded(prior_ICC_shape2, left = 1e-3, right = 1e3)
  assert_bounded(prev_breaks)
  assert_increasing(prev_breaks)
  assert_bounded(ICC_breaks)
  assert_increasing(ICC_breaks)
  
  # make parameter grids
  p <- t(matrix(prev_breaks, nrow = length(prev_breaks), ncol = length(ICC_breaks)))
  rho <- matrix(ICC_breaks, nrow = length(ICC_breaks), ncol = length(prev_breaks))
  
  # define implied beta-binomial parameters
  alpha <- p*(1/rho - 1)
  beta <- (1 - p)*(1/rho - 1)
  
  # evaluate beta-binomial likelihood. NB, we use a new implementation here
  # rather than relying on the dbbinom_reparam() function because we need the
  # likelihood to accept matrix input for both alpha and beta parameters and to
  # be summed over all clusters
  ll <- 0
  for (i in seq_along(n)) {
    ll <- ll + extraDistr::dbbinom(x = n[i], size = N[i], alpha = alpha, beta = beta, log = TRUE)
  }
  ll <- matrix(ll, nrow = length(ICC_breaks))
  
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
  
  # incorporate priors
  ll <- ll + dbeta(p, prior_prev_shape1, prior_prev_shape2, log = TRUE) +
    dbeta(rho, prior_ICC_shape1, prior_ICC_shape2, log = TRUE)
  
  # normalise to largest value and exponentiate (to avoid underflow), then
  # renormalise
  ret <- exp(ll - max(ll))
  ret <- ret / sum(ret)
  
  # return
  return(ret)
}

#------------------------------------------------
#' @title Estimate power when testing prevalence against a threshold
#'
#' @description Estimates power when conducting a clustered prevalence survey
#'   and comparing against a set threshold. Estimates power empirically via
#'   repeated simulation. Returns an estimate of the power, along with lower and
#'   upper 95\% confidence interval of this estimate.
#' 
#' @details
#' Estimates power using the following approach:
#' \enumerate{
#'  \item Simulate data via the function \code{rbbinom_reparam()} using known
#'  values (e.g. a known "true" prevalence and intra-cluster correlation).
#'  \item Analyse data using \code{get_prevalence()} to determine the
#'  probability of being above \code{prev_thresh}.
#'  \item If this probability is above \code{rejection_threshold} then
#'  reject the null hypothesis. Encode this as a single correct conclusion.
#'  \item Repeat steps 1-3 many times. Count the number of simulations for which
#'  the correct conclusion is reached, and divide by the total number of
#'  simulations. This gives an estimate of empirical power, along with upper and
#'  lower 95\% binomial CIs via the method of Clopper and Pearson (1934).
#' }
#' Note that this function can be run even when \code{prevalence} is less than
#' \code{prev_thresh}, although in this case what is returned is not the power.
#' Power is defined as the probability of \emph{correctly} rejecting the null
#' hypothesis, whereas here we would be incorrectly rejecting the null.
#' Therefore, what we obtain in this case is an estimate of the false positive
#' rate.
#' 
#' @inheritParams get_posterior
#' @param N vector of the number of samples obtained from each cluster.
#' @param prevalence assumed true prevalence of pfhrp2/3 deletions as a
#'   proportion between 0 and 1. If a vector of two values is given here then
#'   prevalence is drawn uniformly from between these limits independently for
#'   each simulation. This allows power to be calculated for a composite
#'   hypothesis.
#' @param ICC assumed true intra-cluster correlation (ICC) as a value between 0
#'   and 1.
#' @param prev_thresh the threshold prevalence that we are testing against (5\%
#'   by default).
#' @param rejection_threshold the posterior probability of being above the
#'   prevalence threshold needs to be greater than \code{rejection_threshold} in
#'   order to reject the null hypothesis.
#' @param ICC_infer the value of the ICC assumed in the inference step. If we
#'   plan on estimating the ICC from our data, i.e. running
#'   \code{get_prevalence(ICC = NULL)} (the default), then we should also set
#'   \code{ICC=NULL} here (the default). However, if we plan on running
#'   \code{get_prevalence()} with ICC set to a known value then we should
#'   insert this value here as \code{ICC_infer}.
#' @param reps number of times to repeat simulation per parameter combination.
#' @param silent if \code{TRUE} then suppress all console output.
#'
#' @references
#' Clopper, C.J. and Pearson, E.S., 1934. The use of confidence or fiducial
#' limits illustrated in the case of the binomial. Biometrika, 26, 404–413. doi:
#' 10.2307/2331986.
#'
#' @examples
#' get_power_threshold(N = c(120, 90, 150), prevalence = 0.15, ICC = 0.1 , reps = 1e2)
#' 
#' @importFrom stats runif
#' @importFrom knitrProgressBar progress_estimated update_progress
#' @export

get_power_threshold <- function(N, prevalence = 0.10, ICC = 0.05,
                                prev_thresh = 0.05,
                                rejection_threshold = 0.95,
                                ICC_infer = NULL,
                                prior_prev_shape1 = 1.0, prior_prev_shape2 = 1.0,
                                prior_ICC_shape1 = 1.0, prior_ICC_shape2 = 9.0,
                                n_intervals = 20, round_digits = 2,
                                reps = 1e2, use_cpp = TRUE, silent = FALSE) {
  
  # avoid "no visible binding" note
  n <- NULL
  
  # check inputs
  assert_vector_pos_int(N)
  assert_bounded(prevalence)
  assert_in(length(prevalence), c(1, 2))
  assert_increasing(prevalence)
  assert_single_bounded(ICC)
  assert_vector_bounded(prev_thresh)
  assert_single_bounded(rejection_threshold)
  if (!is.null(ICC_infer)) {
    if (!silent) {
      message(sprintf("Note: this assumes the ICC will be set to the fixed value %s when analysing your data. This overrides the use of the prior on ICC. Use silent = TRUE to suppress this message.", ICC_infer))
    }
    assert_single_bounded(ICC_infer)
  }
  assert_single_bounded(prior_prev_shape1, left = 1e-3, right = 1e3)
  assert_single_bounded(prior_prev_shape2, left = 1e-3, right = 1e3)
  assert_single_bounded(prior_ICC_shape1, left = 1e-3, right = 1e3)
  assert_single_bounded(prior_ICC_shape2, left = 1e-3, right = 1e3)
  assert_single_pos_int(n_intervals)
  assert_greq(n_intervals, 5)
  assert_single_pos_int(round_digits)
  assert_single_pos_int(reps)
  assert_single_logical(use_cpp)
  assert_single_logical(silent)
  
  # prevalence has the option of drawing uniformly between limits
  if (length(prevalence) == 1) {
    prevalence <- rep(prevalence, 2)
  }
  
  # draw data for all simulations. This allows us to group identical datasets,
  # which can save time in some situations (e.g. low prevalence where 0 counts
  # are common)
  l_n <- list()
  for (i in 1:reps) {
    l_n[[i]] <- data.frame(N = N,
                           n = rbbinom_reparam(n_clust = length(N), N = N,
                                               p = runif(1, min = prevalence[1], max = prevalence[2]),
                                               rho = ICC)) %>%
      arrange(N, n)
  }
  
  # group duplicates
  l_u <- unique(l_n)
  l_w <- tabulate(match(l_n, l_u))
  
  # make progress bar
  pb <- knitrProgressBar::progress_estimated(length(l_u))
  
  # simulate
  sim_correct <- matrix(NA, nrow = length(l_u), ncol = length(prev_thresh))
  for (i in seq_along(l_u)) {
    if (!silent) {
      update_progress(pb)
    }
    
    p_est <- get_prevalence(n = l_u[[i]]$n,
                            N = l_u[[i]]$N,
                            alpha = 0.05,
                            prev_thresh = prev_thresh,
                            ICC = ICC_infer,
                            prior_prev_shape1 = prior_prev_shape1,
                            prior_prev_shape2 = prior_prev_shape2,
                            prior_ICC_shape1 = prior_ICC_shape1,
                            prior_ICC_shape2 = prior_ICC_shape2,
                            MAP_on = TRUE,
                            post_mean_on = FALSE,
                            post_median_on = FALSE,
                            post_CrI_on = FALSE,
                            post_thresh_on = TRUE,
                            post_full_on = FALSE,
                            CrI_type = "HDI",
                            n_intervals = n_intervals,
                            round_digits = 5,
                            use_cpp = use_cpp,
                            silent = TRUE)
    
    sim_correct[i,] <- (p_est$prob_above_threshold[[1]] > rejection_threshold)
  }
  
  # weighted sum
  n_correct <- colSums(sim_correct * l_w)
  
  # get 95% CIs on power
  power_CI <- ClopperPearson(n_success = n_correct, n_total = reps, alpha = 0.05)
  ret <- data.frame(prev_thresh = prev_thresh,
                    power = round(n_correct / reps * 100, round_digits),
                    lower = round(power_CI$lower * 100, round_digits),
                    upper = round(power_CI$upper * 100, round_digits))
  rownames(ret) <- NULL
  return(ret)
}

#------------------------------------------------
#' @title Calculate power when testing for presence of deletions
#'
#' @description Calculates power directly for the case of a clustered prevalence
#'   survey where the aim is to detect the presence of *any* deletions over all
#'   clusters. This design can be useful as a pilot study to identify priority
#'   regions where high deletion prevalence is likely. Note that we need to take
#'   account of intra-cluster correlation here, as a high ICC will make it more
#'   likely that we see zero deletions even when the prevalence is non-zero.
#' 
#' @param N vector of the number of samples obtained from each cluster.
#' @param prevalence assumed true prevalence of pfhrp2/3 deletions as a
#'   proportion between 0 and 1.
#' @param ICC assumed true intra-cluster correlation (ICC) between 0 and 1.
#'
#' @examples
#' get_power_presence(N = c(120, 90, 150), prevalence = 0.01, ICC = 0.1)
#' 
#' @export

get_power_presence <- function(N, prevalence = 0.01, ICC = 0.05) {
  
  # check inputs
  assert_vector_pos_int(N, zero_allowed = FALSE)
  assert_single_bounded(prevalence)
  assert_single_bounded(ICC)
  
  # calculate log-probability of seeing zero deletions, dealing with special
  # cases
  if (prevalence == 0) {
    log_prob_zero <- 0
    
  } else if (ICC == 0) {
    log_prob_zero <- sum(N)*log(1 - prevalence)
    
  } else if (ICC == 1) {
    log_prob_zero <- length(N)*log(1 - prevalence)
    
  } else {
    alpha <- prevalence*(1 / ICC - 1)
    beta <- (1 - prevalence)*(1 / ICC - 1)
    log_prob_zero <-  sum(lgamma(alpha + beta) - lgamma(beta) + lgamma(N + beta) - lgamma(N + alpha + beta))
    
  }
  
  # power is probability of not seeing zero deletions
  power <- 1e2*(1 - exp(log_prob_zero))
  
  return(power)
}

#------------------------------------------------
#' @title Get minimum sample size when testing for presence of deletions
#'
#' @description Calculates the minimum sample size required per cluster to
#'   achieve a certain power for the case of a clustered prevalence survey where
#'   the aim is to detect the presence of \emph{any} deletions over all clusters (see
#'   \code{?get_power_presence()}). Assumes the same sample size per cluster.
#' 
#' @param n_clust the number of clusters.
#' @param target_power the power we are aiming to achieve.
#' @param prevalence assumed true prevalence of pfhrp2/3 deletions as a
#'   proportion between 0 and 1.
#' @param ICC assumed true intra-cluster correlation (ICC) between 0 and 1.
#' @param N_max the maximum allowed sample size. Sample sizes are only explored
#'   up to this value, after which point an error is returned.
#'
#' @examples
#' get_sample_size_presence(n_clust = 5, prevalence = 0.01, ICC = 0.1)
#' 
#' @export

get_sample_size_presence <- function(n_clust, target_power = 0.8,
                                     prevalence = 0.01, ICC = 0.05,
                                     N_max = 2e3) {
  
  # check inputs
  assert_single_pos_int(n_clust, zero_allowed = FALSE)
  assert_single_bounded(target_power)
  assert_single_bounded(prevalence)
  assert_single_bounded(ICC)
  assert_single_pos_int(N_max, zero_allowed = FALSE)
  
  # calculate log-probability of seeing zero deletions, dealing with special
  # cases
  N <- 1:N_max
  if (prevalence == 0) {
    log_prob_zero <- 0
    
  } else if (ICC == 0) {
    log_prob_zero <- n_clust*N*log(1 - prevalence)
    
  } else if (ICC == 1) {
    log_prob_zero <- n_clust*log(1 - prevalence)
    
  } else {
    alpha <- prevalence*(1 / ICC - 1)
    beta <- (1 - prevalence)*(1 / ICC - 1)
    log_prob_zero <-  n_clust * (lgamma(alpha + beta) - lgamma(beta) + lgamma(N + beta) - lgamma(N + alpha + beta))
    
  }
  
  # power is probability of not seeing zero deletions
  power <- 1 - exp(log_prob_zero)
  
  # break if power not reached
  if (!any(power > target_power)) {
    stop(sprintf("target power not reached within sample size up to N_max = %s", N_max))
  }
  
  # get minimum sample size
  N_opt <- N[which(power > target_power)[1]]
  
  return(N_opt)
}

