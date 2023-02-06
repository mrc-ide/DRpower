
#------------------------------------------------
# reparameterisation of the beta-binomial distribution in terms of a mean (p)
# and an intra-cluster correlation coefficient (rho). Deals with special cases
# that simplify to the binomial or bernoulli distributions
#' @importFrom extraDistr dbbinom
#' @importFrom stats dbinom
#' @noRd

dbbinom_reparam <- function(k, m, p, rho, log_on = TRUE) {
  
  if (rho == 0) {
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
# draw from reparameterisation of the beta-binomial distribution in terms of a
# mean (p) and an intra-cluster correlation coefficient (rho). Deals with
# special cases that simplify to the binomial distribution
#' @importFrom extraDistr rbbinom
#' @importFrom stats rbinom
#' @noRd

rbbinom_reparam <- function(n, m, p, rho) {
  if (rho == 0) {
    # simplifies to binomial distribution
    ret <- rbinom(n = n, size = m, prob = p)
    
  } else {
    # beta-binomial distribution
    alpha <- p * (1 - rho) / rho
    beta <- (1 - p) * (1 - rho) / rho
    ret <- extraDistr::rbbinom(n = n, size = m, alpha = alpha, beta = beta)
    
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
    sum(dbbinom_reparam(k = n, m = N, p = p, rho = rho, log_on = TRUE))
}

#------------------------------------------------
# calculate log(area) via trapezoidal rule
#' @noRd

get_area_single_trap <- function(x0, x1, log_y0, log_y1) {
  if ((log_y0 == -Inf) & (log_y1 == -Inf)) {
    ret <- -Inf
  } else if (log_y0 > log_y1) {
    ret <- log_y0 + log(1 + exp(log_y1 - log_y0)) + log(x1 - x0) - log(2)
  } else {
    ret <- log_y1 + log(1 + exp(log_y0 - log_y1)) + log(x1 - x0) - log(2)
  }
  return(ret)
}

#------------------------------------------------
# calculate log(area) given two points and a midpoint via trapezoidal rule
#' @noRd

get_area_trap <- function(x0, x1, log_y0, log_ym, log_y1) {
  l1 <- get_area_single_trap(x0, 0.5*(x0 + x1), log_y0, log_ym)
  l2 <- get_area_single_trap(0.5*(x0 + x1), x1, log_ym, log_y1)
  if ((l1 == -Inf) & (l2 == -Inf)) {
    return(-Inf)
  } else if (l1 > l2) {
    ret <- l1 + log(1 + exp(l2 - l1))
  } else {
    ret <- l2 + log(1 + exp(l1 - l2))
  }
  return(ret)
}

#------------------------------------------------
# calculate log(area) via Simpson's rule
#' @noRd

get_area_Simp <- function(x0, x1, log_y0, log_ym, log_y1) {
  if ((log_y0 == -Inf) & (log_ym == -Inf) & (log_y1 == -Inf)) {
    ret <- -Inf
  } else {
    z <- c(log_y0, log(4) + log_ym, log_y1)
    ret <- max(z) + log(sum(exp(z - max(z)))) + log(x1 - x0) - log(6)
  }
  return(ret)
}

#------------------------------------------------
# calculate A coefficient in the Simpson's rule formula Ax^2 + Bx + C
#' @noRd

get_Simp_A <- function(x0, xm, x1, log_y0, log_ym, log_y1) {
  c1 <- xm - x0
  c2 <- x1 - xm
  (c1*(exp(log_y1) - exp(log_ym)) - c2*(exp(log_ym) - exp(log_y0))) / (c1*(x1^2 - xm^2) + c2*(x0^2 - xm^2))
}

#------------------------------------------------
# calculate B coefficient in the Simpson's rule formula Ax^2 + Bx + C
#' @noRd

get_Simp_B <- function(x0, xm, log_y0, log_ym, A) {
  (exp(log_ym) - exp(log_y0) + A*(x0^2 - xm^2)) / (xm - x0)
}

#------------------------------------------------
# calculate C coefficient in the Simpson's rule formula Ax^2 + Bx + C
#' @noRd

get_Simp_C <- function(x0, log_y0, A, B) {
  exp(log_y0) - A*x0^2 - B*x0
}

#------------------------------------------------
# get lowest point in Simpson's rule quadratic
#' @noRd

get_Simp_lowest <- function(x0, x1, A, B, C) {
  
  # get y-values at both ends of interval and point of inflection
  y0 <- A*x0^2 + B*x0 + C
  y1 <- A*x1^2 + B*x1 + C
  x_inflect <- -B / (2*A)
  y_inflect <- A*x_inflect^2 + B*x_inflect + C
  
  # get smallest y-value. This can be the point of inflection as long as this is
  # inside the range [x0, x1]
  y_min <- pmin(y0, y1)
  w <- which((x_inflect > x0) & (x_inflect < x1))
  if (any(w)) {
    y_min[w] <- y_inflect[w]
  }
  
  return(y_min)
}

#------------------------------------------------
# calculate gradient at midpoint and over range, and return absolute relative
# difference in gradient
#' @noRd

compare_grad <- function(x0, x1, log_y0, log_y1, log_ym, log_ydelta, delta) {
  if (!is.finite(log_y0) && !is.finite(log_ym) && !is.finite(log_y1)) {
    ret <- 1
  } else {
    grad_mid <- (exp(log_ydelta) - exp(log_ym)) / delta;
    grad_trap <- (exp(log_y1) - exp(log_y0)) / (x1 - x0);
    ret <- abs(grad_mid / grad_trap - 1.0);
  }
  return(ret)
}

#------------------------------------------------
# general function for performing integration via adaptive quadrature
# n_intervals must be 2 or greater
#' @noRd

adaptive_quadrature <- function(f1, n_intervals, left, right, debug_on = FALSE) {
  
  delta <- 1e-4
  
  # initialise empty dataframe for storing results
  df_quad <- data.frame(x0 = rep(NA, n_intervals), xm = NA, x1 = NA,
                        log_y0 = NA, log_ym = NA, log_y1 = NA,
                        log_area_trap = NA, log_area_Simp = NA,
                        log_area_diff = NA)
  
  # create first interval
  df_quad$x0[1] <- left
  df_quad$x1[1] <- right
  df_quad$xm[1] <- 0.5*(df_quad$x0[1] + df_quad$x1[1])
  df_quad$log_y0[1] <- f1(df_quad$x0[1])
  df_quad$log_y1[1] <- f1(df_quad$x1[1])
  df_quad$log_ym[1] <- f1(df_quad$xm[1])
  df_quad$log_area_trap[1] <- get_area_trap(df_quad$x0[1], df_quad$x1[1], df_quad$log_y0[1], df_quad$log_ym[1], df_quad$log_y1[1])
  df_quad$log_area_Simp[1] <- get_area_Simp(df_quad$x0[1], df_quad$x1[1], df_quad$log_y0[1], df_quad$log_ym[1], df_quad$log_y1[1])
  df_quad$log_area_diff[1] <- 0
  
  # loop through remaining n_intervals
  for (i in 2:n_intervals) {
    
    # find largest discrepancy in area
    w <- which.max(df_quad$log_area_diff)
    
    # create new entry and modify existing
    df_quad$x0[i] <- df_quad$xm[w]
    df_quad$x1[i] <- df_quad$x1[w]
    df_quad$xm[i] <- 0.5*(df_quad$x0[i] + df_quad$x1[i])
    df_quad$log_y0[i] <- df_quad$log_ym[w]
    df_quad$log_y1[i] <- df_quad$log_y1[w]
    df_quad$log_ym[i] <- f1(df_quad$xm[i])
    df_quad$log_area_trap[i] <- get_area_trap(df_quad$x0[i], df_quad$x1[i], df_quad$log_y0[i], df_quad$log_ym[i], df_quad$log_y1[i])
    df_quad$log_area_Simp[i] <- get_area_Simp(df_quad$x0[i], df_quad$x1[i], df_quad$log_y0[i], df_quad$log_ym[i], df_quad$log_y1[i])
    log_ydelta = f1(df_quad$xm[i] + delta);
    grad_diff <- compare_grad(df_quad$x0[i], df_quad$x1[i], df_quad$log_y0[i], df_quad$log_y1[i], df_quad$log_ym[i], log_ydelta, delta);
    df_quad$log_area_diff[i] <- grad_diff * abs(exp(df_quad$log_area_trap[i]) - exp(df_quad$log_area_Simp[i]))
    
    df_quad$x1[w] = df_quad$xm[w]
    df_quad$xm[w] = 0.5*(df_quad$x0[w] + df_quad$x1[w])
    df_quad$log_y1[w] <- df_quad$log_ym[w]
    df_quad$log_ym[w] <- f1(df_quad$xm[w])
    df_quad$log_area_trap[w] <- get_area_trap(df_quad$x0[w], df_quad$x1[w], df_quad$log_y0[w], df_quad$log_ym[w], df_quad$log_y1[w])
    df_quad$log_area_Simp[w] <- get_area_Simp(df_quad$x0[w], df_quad$x1[w], df_quad$log_y0[w], df_quad$log_ym[w], df_quad$log_y1[w])
    log_ydelta = f1(df_quad$xm[w] + delta);
    grad_diff <- compare_grad(df_quad$x0[w], df_quad$x1[w], df_quad$log_y0[w], df_quad$log_y1[w], df_quad$log_ym[w], log_ydelta, delta);
    df_quad$log_area_diff[w] <- grad_diff * abs(exp(df_quad$log_area_trap[w]) - exp(df_quad$log_area_Simp[w]))
  }
  
  # debug distribution
  if (debug_on) {
    
    # get in order of x0
    df_order <- df_quad[order(df_quad$x0),]
    
    # calculate Simpson's rule coefficients
    df_order$A <- get_Simp_A(df_order$x0, df_order$xm, df_order$x1, df_order$log_y0, df_order$log_ym, df_order$log_y1)
    df_order$B <- get_Simp_B(df_order$x0, df_order$xm, df_order$log_y0, df_order$log_ym, df_order$A)
    df_order$C <- get_Simp_C(df_order$x0, df_order$log_y0, df_order$A, df_order$B)
    
    # create sequence of x-values spanning total range and evaluate function by
    # brute force over this range
    x <- seq(left, right, l = 201)
    fx <- rep(NA, length(x))
    for (i in seq_along(x)) {
      fx[i] <- exp(f1(x[i]))
    }
    
    # calculate interpolated curve from Simpson's rule
    z <- findInterval(x, vec = df_order$x0)
    fx2 <- df_order$A[z]*x^2 + df_order$B[z]*x + df_order$C[z]
    
    # produce plot
    plot(x, fx, type = 'l', lwd = 2)
    lines(x, fx2, col = 4, lty = 2, lwd = 2)
    points(df_quad$x0, exp(df_quad$log_y0), pch = 20, cex = 0.75)
    points(df_quad$xm, exp(df_quad$log_ym), pch = 20, cex = 0.75, col = 4)
    segments(df_quad$x0, exp(df_quad$log_y0), df_quad$xm, exp(df_quad$log_ym), col = 2)
    segments(df_quad$xm, exp(df_quad$log_ym), df_quad$x1, exp(df_quad$log_y1), col = 2)
    abline(h = 0, lty = 3)
  }
  
  return(df_quad)
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
                             n_intervals = 20, debug_on = FALSE) {
  
  # check inputs
  assert_vector_pos_int(n)
  assert_vector_pos_int(N)
  if (length(N) == 1) {
    N <- rep(N, length(n))
  }
  assert_same_length(n, N, message = "N must be either a single value (all clusters the same size) or a vector with the same length as n")
  assert_single_bounded(alpha, inclusive_left = FALSE, inclusive_right = FALSE)
  assert_single_bounded(prior_prev_shape1, left = 1, right = 1e3)
  assert_single_bounded(prior_prev_shape2, left = 1, right = 1e3)
  assert_single_bounded(prior_ICC_shape1, left = 1, right = 1e3)
  assert_single_bounded(prior_ICC_shape2, left = 1, right = 1e3)
  assert_vector_pos_int(n_intervals)
  assert_greq(n_intervals, 5)
  assert_single_logical(debug_on)
  
  # approximate distribution of rho via adaptive quadrature
  df_quad <- adaptive_quadrature(f1 = function(rho) {
    loglike_marginal_rho(n = n, N = N, rho = rho, n_intervals = n_intervals,
                         prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                         prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
  }, n_intervals = n_intervals, left = 0, right = 1, debug_on = debug_on)
  
  # reorder in terms of increasing rho
  df_quad <- df_quad[order(df_quad$x0),]
  
  # normalise interval areas
  log_area_Simp <- df_quad$log_area_Simp
  log_area_sum <- max(log_area_Simp) + log(sum(exp(log_area_Simp - max(log_area_Simp))))
  area <- exp(log_area_Simp - log_area_sum)
  area_cs <- cumsum(area)
  
  # solve for lower and upper CrIs
  w1 <- which(area_cs > alpha / 2)[1]
  area_remaining1 <- alpha / 2 - (area_cs[w1] - area[w1])
  CrI_lower <- solve_Simpsons_area(a = df_quad$x0[w1], b = df_quad$x1[w1],
                                   fa = exp(df_quad$log_y0[w1] - log_area_sum),
                                   fm = exp(df_quad$log_ym[w1] - log_area_sum),
                                   fb = exp(df_quad$log_y1[w1] - log_area_sum),
                                   target_area = area_remaining1)
  
  w2 <- which(area_cs > (1 - alpha / 2))[1]
  area_remaining2 <- area_cs[w2] - (1 - alpha / 2)
  CrI_upper <- solve_Simpsons_area(a = df_quad$x0[w2], b = df_quad$x1[w2],
                                   fa = exp(df_quad$log_y0[w2] - log_area_sum),
                                   fm = exp(df_quad$log_ym[w2] - log_area_sum),
                                   fb = exp(df_quad$log_y1[w2] - log_area_sum),
                                   target_area = area_remaining2)
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

get_credible_prevalence <- function(n, N, alpha = 0.05,
                                    prior_prev_shape1 = 1.0, prior_prev_shape2 = 1.0,
                                    prior_ICC_shape1 = 1.0, prior_ICC_shape2 = 1.0,
                                    n_intervals = 20, debug_on = FALSE) {
  
  # check inputs
  assert_vector_pos_int(n)
  assert_vector_pos_int(N)
  if (length(N) == 1) {
    N <- rep(N, length(n))
  }
  assert_same_length(n, N, message = "N must be either a single value (all clusters the same size) or a vector with the same length as n")
  assert_single_bounded(alpha, inclusive_left = FALSE, inclusive_right = FALSE)
  assert_single_bounded(prior_prev_shape1, left = 1, right = 1e3)
  assert_single_bounded(prior_prev_shape2, left = 1, right = 1e3)
  assert_single_bounded(prior_ICC_shape1, left = 1, right = 1e3)
  assert_single_bounded(prior_ICC_shape2, left = 1, right = 1e3)
  assert_vector_pos_int(n_intervals)
  assert_greq(n_intervals, 5)
  assert_single_logical(debug_on)
  
  # approximate distribution of p via adaptive quadrature
  df_quad <- adaptive_quadrature(f1 = function(p) {
    loglike_marginal_p(n = n, N = N, p = p, n_intervals = n_intervals,
                       prior_p_shape1 = prior_prev_shape1, prior_p_shape2 = prior_prev_shape2,
                       prior_rho_shape1 = prior_ICC_shape1, prior_rho_shape2 = prior_ICC_shape2)
  }, n_intervals = n_intervals, left = 0, right = 1, debug_on = debug_on)
  
  # reorder in terms of increasing p
  df_quad <- df_quad[order(df_quad$x0),]
  
  # normalise interval areas
  log_area_Simp <- df_quad$log_area_Simp
  log_area_sum <- max(log_area_Simp) + log(sum(exp(log_area_Simp - max(log_area_Simp))))
  area <- exp(log_area_Simp - log_area_sum)
  area_cs <- cumsum(area)
  
  # solve for lower and upper CrIs
  w1 <- which(area_cs > alpha / 2)[1]
  area_remaining1 <- alpha / 2 - (area_cs[w1] - area[w1])
  CrI_lower <- solve_Simpsons_area(a = df_quad$x0[w1], b = df_quad$x1[w1],
                                   fa = exp(df_quad$log_y0[w1] - log_area_sum),
                                   fm = exp(df_quad$log_ym[w1] - log_area_sum),
                                   fb = exp(df_quad$log_y1[w1] - log_area_sum),
                                   target_area = area_remaining1)
  
  w2 <- which(area_cs > (1 - alpha / 2))[1]
  area_remaining2 <- area_cs[w2] - (1 - alpha / 2)
  CrI_upper <- solve_Simpsons_area(a = df_quad$x0[w2], b = df_quad$x1[w2],
                                   fa = exp(df_quad$log_y0[w2] - log_area_sum),
                                   fm = exp(df_quad$log_ym[w2] - log_area_sum),
                                   fb = exp(df_quad$log_y1[w2] - log_area_sum),
                                   target_area = area_remaining2)
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
#'   but implemented in Rcpp for increased speed.
#'
#' @inheritParams get_credible_prevalence
#'
#' @export

get_credible_prevalence_fast <- function(n, N, alpha = 0.05,
                                         prior_prev_shape1 = 1.0, prior_prev_shape2 = 1.0,
                                         prior_ICC_shape1 = 1.0, prior_ICC_shape2 = 1.0,
                                         n_intervals = 20, debug_on = FALSE) {
  
  # check inputs
  assert_vector_pos_int(n)
  assert_vector_pos_int(N)
  if (length(N) == 1) {
    N <- rep(N, length(n))
  }
  assert_same_length(n, N, message = "N must be either a single value (all clusters the same size) or a vector with the same length as n")
  assert_single_bounded(alpha, inclusive_left = FALSE, inclusive_right = FALSE)
  assert_single_bounded(prior_prev_shape1, left = 1, right = 1e3)
  assert_single_bounded(prior_prev_shape2, left = 1, right = 1e3)
  assert_single_bounded(prior_ICC_shape1, left = 1, right = 1e3)
  assert_single_bounded(prior_ICC_shape2, left = 1, right = 1e3)
  assert_vector_pos_int(n_intervals)
  assert_greq(n_intervals, 5)
  assert_single_logical(debug_on)
  
  # get arguments into list
  args_params <- list(n = n,
                      N = rep(N, length(n)),
                      alpha = alpha,
                      prior_prev_shape1 = prior_prev_shape1,
                      prior_prev_shape2 = prior_prev_shape2,
                      prior_ICC_shape1 = prior_ICC_shape1,
                      prior_ICC_shape2 = prior_ICC_shape2,
                      n_intervals = n_intervals)
  
  # run efficient C++ function
  output_raw <- get_credible_prevalence_cpp(args_params)
  df_quad <- as.data.frame(output_raw)
  
  # reorder in terms of increasing p
  df_quad <- df_quad[order(df_quad$x0),]
  
  # normalise interval areas
  log_area_Simp <- df_quad$log_area_Simp
  log_area_sum <- max(log_area_Simp) + log(sum(exp(log_area_Simp - max(log_area_Simp))))
  area <- exp(log_area_Simp - log_area_sum)
  area_cs <- cumsum(area)
  
  # solve for lower and upper CrIs
  w1 <- which(area_cs > alpha / 2)[1]
  area_remaining1 <- alpha / 2 - (area_cs[w1] - area[w1])
  CrI_lower <- solve_Simpsons_area(a = df_quad$x0[w1], b = df_quad$x1[w1],
                                   fa = exp(df_quad$log_y0[w1] - log_area_sum),
                                   fm = exp(df_quad$log_ym[w1] - log_area_sum),
                                   fb = exp(df_quad$log_y1[w1] - log_area_sum),
                                   target_area = area_remaining1)
  
  w2 <- which(area_cs > (1 - alpha / 2))[1]
  area_remaining2 <- area_cs[w2] - (1 - alpha / 2)
  CrI_upper <- solve_Simpsons_area(a = df_quad$x0[w2], b = df_quad$x1[w2],
                                   fa = exp(df_quad$log_y0[w2] - log_area_sum),
                                   fm = exp(df_quad$log_ym[w2] - log_area_sum),
                                   fb = exp(df_quad$log_y1[w2] - log_area_sum),
                                   target_area = area_remaining2)
  
  # debug distribution
  if (debug_on) {
    
    # get in order of x0
    df_order <- df_quad[order(df_quad$x0),]
    
    # calculate Simpson's rule coefficients
    df_order$A <- get_Simp_A(df_order$x0, df_order$xm, df_order$x1, df_order$log_y0, df_order$log_ym, df_order$log_y1)
    df_order$B <- get_Simp_B(df_order$x0, df_order$xm, df_order$log_y0, df_order$log_ym, df_order$A)
    df_order$C <- get_Simp_C(df_order$x0, df_order$log_y0, df_order$A, df_order$B)
    
    # calculate interpolated curve from Simpson's rule
    x <- seq(0, 1, l = 201)
    z <- findInterval(x, vec = df_order$x0)
    fx <- df_order$A[z]*x^2 + df_order$B[z]*x + df_order$C[z]
    
    # produce plot
    plot(x, fx, type = 'l', col = 4, lty = 2, lwd = 1)
    points(df_quad$x0, exp(df_quad$log_y0), pch = 20, cex = 0.75)
    points(df_quad$xm, exp(df_quad$log_ym), pch = 20, cex = 0.75, col = 4)
    segments(df_quad$x0, exp(df_quad$log_y0), df_quad$xm, exp(df_quad$log_ym), col = 2)
    segments(df_quad$xm, exp(df_quad$log_ym), df_quad$x1, exp(df_quad$log_y1), col = 2)
    abline(h = 0, lty = 3)
  }
  
  # return
  return(c(lower = CrI_lower, upper = CrI_upper))
}

#------------------------------------------------
#' @title Rcpp implementation of \code{get_credible_ICC()}
#'
#' @description Equivalent to \code{get_credible_ICC()} in functionality, but
#'   implemented in Rcpp for increased speed.
#'
#' @inheritParams get_credible_prevalence
#'
#' @export

get_credible_ICC_fast <- function(n, N, alpha = 0.05,
                                  prior_prev_shape1 = 1.0, prior_prev_shape2 = 1.0,
                                  prior_ICC_shape1 = 1.0, prior_ICC_shape2 = 1.0,
                                  n_intervals = 20, debug_on = FALSE) {
  
  # check inputs
  assert_vector_pos_int(n)
  assert_vector_pos_int(N)
  if (length(N) == 1) {
    N <- rep(N, length(n))
  }
  assert_same_length(n, N, message = "N must be either a single value (all clusters the same size) or a vector with the same length as n")
  assert_single_bounded(alpha, inclusive_left = FALSE, inclusive_right = FALSE)
  assert_single_bounded(prior_prev_shape1, left = 1, right = 1e3)
  assert_single_bounded(prior_prev_shape2, left = 1, right = 1e3)
  assert_single_bounded(prior_ICC_shape1, left = 1, right = 1e3)
  assert_single_bounded(prior_ICC_shape2, left = 1, right = 1e3)
  assert_vector_pos_int(n_intervals)
  assert_greq(n_intervals, 5)
  assert_single_logical(debug_on)
  
  # get arguments into list
  args_params <- list(n = n,
                      N = rep(N, length(n)),
                      alpha = alpha,
                      prior_prev_shape1 = prior_prev_shape1,
                      prior_prev_shape2 = prior_prev_shape2,
                      prior_ICC_shape1 = prior_ICC_shape1,
                      prior_ICC_shape2 = prior_ICC_shape2,
                      n_intervals = n_intervals)
  
  # run efficient C++ function
  output_raw <- get_credible_ICC_cpp(args_params)
  df_quad <- as.data.frame(output_raw)
  
  # reorder in terms of increasing p
  df_quad <- df_quad[order(df_quad$x0),]
  
  # normalise interval areas
  log_area_Simp <- df_quad$log_area_Simp
  log_area_sum <- max(log_area_Simp) + log(sum(exp(log_area_Simp - max(log_area_Simp))))
  area <- exp(log_area_Simp - log_area_sum)
  area_cs <- cumsum(area)
  
  # solve for lower and upper CrIs
  w1 <- which(area_cs > alpha / 2)[1]
  area_remaining1 <- alpha / 2 - (area_cs[w1] - area[w1])
  CrI_lower <- solve_Simpsons_area(a = df_quad$x0[w1], b = df_quad$x1[w1],
                                   fa = exp(df_quad$log_y0[w1] - log_area_sum),
                                   fm = exp(df_quad$log_ym[w1] - log_area_sum),
                                   fb = exp(df_quad$log_y1[w1] - log_area_sum),
                                   target_area = area_remaining1)
  
  w2 <- which(area_cs > (1 - alpha / 2))[1]
  area_remaining2 <- area_cs[w2] - (1 - alpha / 2)
  CrI_upper <- solve_Simpsons_area(a = df_quad$x0[w2], b = df_quad$x1[w2],
                                   fa = exp(df_quad$log_y0[w2] - log_area_sum),
                                   fm = exp(df_quad$log_ym[w2] - log_area_sum),
                                   fb = exp(df_quad$log_y1[w2] - log_area_sum),
                                   target_area = area_remaining2)
  
  # debug distribution
  if (debug_on) {
    
    # get in order of x0
    df_order <- df_quad[order(df_quad$x0),]
    
    # calculate Simpson's rule coefficients
    df_order$A <- get_Simp_A(df_order$x0, df_order$xm, df_order$x1, df_order$log_y0, df_order$log_ym, df_order$log_y1)
    df_order$B <- get_Simp_B(df_order$x0, df_order$xm, df_order$log_y0, df_order$log_ym, df_order$A)
    df_order$C <- get_Simp_C(df_order$x0, df_order$log_y0, df_order$A, df_order$B)
    
    # calculate interpolated curve from Simpson's rule
    x <- seq(0, 1, l = 201)
    z <- findInterval(x, vec = df_order$x0)
    fx <- df_order$A[z]*x^2 + df_order$B[z]*x + df_order$C[z]
    
    # produce plot
    plot(x, fx, type = 'l', col = 4, lty = 2, lwd = 1)
    points(df_quad$x0, exp(df_quad$log_y0), pch = 20, cex = 0.75)
    points(df_quad$xm, exp(df_quad$log_ym), pch = 20, cex = 0.75, col = 4)
    segments(df_quad$x0, exp(df_quad$log_y0), df_quad$xm, exp(df_quad$log_ym), col = 2)
    segments(df_quad$xm, exp(df_quad$log_ym), df_quad$x1, exp(df_quad$log_y1), col = 2)
    abline(h = 0, lty = 3)
  }
  
  # return
  return(c(lower = CrI_lower, upper = CrI_upper))
}
