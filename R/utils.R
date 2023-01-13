
#------------------------------------------------
#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL


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
# draw from reparameterisation of the beta-binomial distribution in terms of a
# mean (p) and an intra-cluster correlation coefficient (rho). Deals with
# special cases that simplify to the binomial distribution
#' @importFrom extraDistr rbbinom
#' @importFrom stats rbinom
#' @noRd

rbbinom_reparam <- function(n, m, p, rho) {
  if ((rho == 0) || (m == 1)) {
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
# apply Simpson's rule to interpolate between a series of node x- and y-values.
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
    stop("could not find single root to cubic expression in Simpson's rule within the defined interval")
  }
  ret <- Re(roots)[w]
  
  return(ret)
}


