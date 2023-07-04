
#------------------------------------------------
# take a vector x of log values and calculate log(sum(exp(x))) in an
# underflow-safe way
#' @noRd

log_sum <- function(x) {
  if (all(x == -Inf)) {
    return(-Inf)
  }
  mx <- max(x)
  mx + log(sum(exp(x - mx)))
}

#------------------------------------------------
# calculate A coefficient in the Simpson's rule formula Ax^2 + Bx + C.
# accepts vector inputs
#' @noRd

get_Simp_A <- function(x0, xm, x1, log_y0, log_ym, log_y1) {
  c1 <- xm - x0
  c2 <- x1 - xm
  (c1*(exp(log_y1) - exp(log_ym)) - c2*(exp(log_ym) - exp(log_y0))) / (c1*(x1^2 - xm^2) + c2*(x0^2 - xm^2))
}

#------------------------------------------------
# calculate B coefficient in the Simpson's rule formula Ax^2 + Bx + C.
# accepts vector inputs
#' @noRd

get_Simp_B <- function(x0, xm, log_y0, log_ym, A) {
  (exp(log_ym) - exp(log_y0) + A*(x0^2 - xm^2)) / (xm - x0)
}

#------------------------------------------------
# calculate C coefficient in the Simpson's rule formula Ax^2 + Bx + C.
# accepts vector inputs
#' @noRd

get_Simp_C <- function(x0, log_y0, A, B) {
  exp(log_y0) - A*x0^2 - B*x0
}

#------------------------------------------------
# evaluate Simpson's quadratic Ax^2 + Bx + C.
# accepts vector inputs
#' @noRd

get_Simp_y <- function(x, A, B, C) {
  A*x^2 + B*x + C
}

#------------------------------------------------
# solve Simpson's quadratic y = Ax^2 + Bx + C for x given y in a given x range.
# The first solution within the range is returned, meaning if there is a single
# solution within the range this is unique, otherwise the choice is arbitrary
# (hence, this should only be used when there is a single unique solution). If
# no solution then error
#' @noRd

get_Simp_x <- function(y, x0, x1, A, B, C) {
  descriminant <- B^2 - 4*A*(C - y)
  if (descriminant < 0) {
    stop("descriminant < 0")
  }
  x_prop <- (-B + c(-1, 1)*sqrt(descriminant)) / (2*A)
  w <- which((x_prop + 1e-10 >= x0) & (x_prop - 1e-10 <= x1))
  if (!any(w)) {
    stop("no solution to quadratic inside range")
  }
  x <- x_prop[w[1]]
  return(x)
}

#------------------------------------------------
# calculate log(area) of a single trapezium given x and log(y) coordinates
#' @noRd

get_trap_area <- function(x0, x1, log_y0, log_y1) {
  log_sum(c(log_y0, log_y1)) + log(x1 - x0) - log(2)
}

#------------------------------------------------
# calculate log(area) via trapezoidal rule given three x coordinates and three
# log(y) coordinates. Total area is the sum of two trapezia
#' @noRd

get_trap_area_double <- function(x0, xm, x1, log_y0, log_ym, log_y1) {
  l1 <- get_trap_area(x0, xm, log_y0, log_ym)
  l2 <- get_trap_area(xm, x1, log_ym, log_y1)
  log_sum(c(l1, l2))
}

#------------------------------------------------
# calculate log(area) via Simpson's rule given x corrdinates at ends and log(y)
# coordinates at ends and midpoint
#' @noRd

get_Simp_area <- function(x0, x1, log_y0, log_ym, log_y1) {
  log_sum(c(log_y0, log(4) + log_ym, log_y1)) + log(x1 - x0) - log(6)
}

#------------------------------------------------
# integrate Simpson's quadratic Ax^2 + Bx + C between x0 and x1. Direction of
# integration does not matter, meaning area returned is always positive.
# accepts vector inputs. 
#' @noRd

get_Simp_integral <- function(x0, x1, A, B, C) {
  ifelse(x1 > x0,
         1/3*A*(x1^3 - x0^3) + 1/2*B*(x1^2 - x0^2) + C*(x1 - x0),
         1/3*A*(x0^3 - x1^3) + 1/2*B*(x0^2 - x1^2) + C*(x0 - x1))
}

#------------------------------------------------
# solve Simpson's rule over the interval [x0, x1] to find the x-value at which
# the area under the curve equals target_area
#' @noRd

solve_Simp_area <- function(x0, x1, A, B, C, target_area) {
  
  # define a new series of coefficients for the integral of this quadratic
  # equation from 0 to z. The solution can be written:
  # A2*z^3 + B2*z^2 + C2*z + D2
  A2 <- A / 3
  B2 <- B / 2
  C2 <- C
  D2 <- -A2*x0^3 - B2*x0^2 - C2*x0
  
  # solve for the z value at which the area under curve equals the target area
  roots <- polyroot(c(D2 - target_area, C2, B2, A2))
  real_roots <- Re(roots)
  
  # most often there will be a single root within the interval
  if (any((real_roots >= x0) & (real_roots <= x1))) {
    w <- which((real_roots >= x0) & (real_roots <= x1))[1]
    ret <- real_roots[w]
    
  } else { # underflow can cause problems when the root is exactly at the border
    
    eq_x0 <- mapply(function(z) isTRUE(all.equal(z, x0)), real_roots)
    if (any(eq_x0)) {
      w <- which(eq_x0)[1]
      ret <- real_roots[w]
    } else {
      eq_x1 <- mapply(function(z) isTRUE(all.equal(z, x1)), real_roots)
      if (any(eq_x1)) {
        w <- which(eq_x1)[1]
        ret <- real_roots[w]
      } else {
        stop("function solve_Simp_area() unable to find roots to polynomial within defined interval")
      }
    }
  }
  
  return(ret)
}

#------------------------------------------------
# get lowest point in Simpson's rule quadratic between x0 and x1.
# accepts vector inputs
#' @noRd

get_Simp_lowest <- function(x0, x1, A, B, C) {
  
  # get y-values at both ends of interval and point of inflection
  y0 <- get_Simp_y(x0, A, B, C)
  y1 <- get_Simp_y(x1, A, B, C)
  x_inflect <- -B / (2*A)
  y_inflect <- get_Simp_y(x_inflect, A, B, C)
  
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
# given two log(area) values, calculate the log(area) of the difference between
# them
#' @noRd

get_log_area_diff <- function(log_area1, log_area2) {
  if ((log_area1 == -Inf) & (log_area2 == -Inf)) {
    ret <- -Inf
  } else if (log_area1 > log_area2) {
    ret <- log_area1 + log(1 - exp(log_area2 - log_area1))
  } else {
    ret <- log_area2 + log(1 - exp(log_area1 - log_area2))
  }
  return(ret)
}

#------------------------------------------------
# calculate gradient in two ways. First, between points at the ends of the range
# [x0, x1], with corresponding log(y) values. Second, between points at [xm, xm
# + delta] with corresponding log(y) values. Calculate both gradients as angles,
# and return the difference normalised such that identical angles return 0 and
# maximally different angles (180 degrees apart) return 1.
#' @noRd

get_grad_diff <- function(x0, x1, log_y0, log_y1, log_ym, log_yd, delta) {
  
  # get two angles between [-0.5, 0.5] and subtract to get relative difference,
  # then normalise to [0,1] scale where 0 represents no differece and 1
  # represents maximal difference
  theta1 <- atan((exp(log_y1) - exp(log_y0)) / (x1 - x0)) / pi
  theta2 <- atan((exp(log_yd) - exp(log_ym)) / (delta*(x1 - x0)/2)) / pi
  ret <- abs(theta1 - theta2)
  return(ret)
}

#------------------------------------------------
# combines values to produce total diff score
#' @noRd

get_total_diff <- function(log_area_diff, grad_diff) {
  exp(log_area_diff) * grad_diff
}

#------------------------------------------------
#' @title Break a function into intervals using adaptive quadrature
#'
#' @description Takes a function \code{f1(x)} that returns values in log space
#'   (e.g. a log-likelihood) and adaptively chooses x-values such that the final
#'   discretised curve is reasonably smooth. The adaptive method uses Simpson's
#'   rule, meaning the final curve can be approximated by a series of
#'   quadratics. Returns a data.frame with all the intervals used in quadrature.
#'   This method is only appropriate for simple, unimodal distributions and is
#'   not guaranteed to work for every possible distribution shape.
#'
#' @details The adaptive method works by breaking the curve into a series of
#'   intervals, and for each interval evaluating the area via both trapezoidal
#'   rule and Simpson's rule. If an interval has a large discrepancy in area
#'   between the two methods then this implies that the function is not very
#'   linear over this interval and so it is a good target for further
#'   subdivision. However, there is an important edge-case in which the midpoint
#'   of the interval happens to have a corresponding function value that is
#'   exactly at the midpoint of the y-range by chance, in which case the
#'   function would appear linear. For this reason, gradient information is also
#'   taken into account. The gradient is evaluated over the full interval and at
#'   the midpoint (searching a small distance delta from the midpoint). If the
#'   gradients are very different this also suggests the function is non-linear
#'   over the interval and so is a target for further subdivision. Both area-
#'   and gradient-based information are combined to produce a total score, and
#'   the interval with the highest score is split into two subintervals about
#'   the midpoint. This process is repeated \code{n_intervals} times. This
#'   method ensures that breakpoints are placed at areas where 1) the curve is
#'   changing rapidly, and 2) the area is large. Hence, it does a decent job of
#'   approximating the total area under the curve, and can also be used to
#'   calculate tails of the distribution reasonably accurately.
#'
#' @param f1 a function returning values in log space. Effort has been taken to
#'   ensure the quadrature method is relatively underflow-safe, working with log
#'   values wherever possible.
#' @param n_intervals the final number of intervals resulting from the
#'   quadrature process.
#' @param left,right the range of the interval to explore.
#' @param debug_on if TRUE then produces a plot comparing the approximating
#'   function with the true function evaluated over a grid.
#' @param delta the small distance used when calculating gradients near the
#'   midpoint. This value is used as a proportion of the interval, not as an
#'   absolute value, meaning values can never be proposed outside of the
#'   quadrature range.
#'
# do not document to keep package help simple
# @export
#' @noRd

adaptive_quadrature <- function(f1, n_intervals, left, right, debug_on = FALSE, delta = 1e-3) {
  
  # initialise empty dataframe for storing results
  df_quad <- data.frame(x0 = rep(NA, n_intervals), xm = NA, xd = NA, x1 = NA,
                        log_y0 = NA, log_ym = NA, log_yd = NA, log_y1 = NA,
                        log_area_trap = NA, log_area_Simp = NA,
                        log_area_diff = NA, grad_diff = NA, total_diff = NA)
  
  # create first interval over entire domain
  df_quad$x0[1] <- left
  df_quad$x1[1] <- right
  df_quad$xm[1] <- 0.5*(left + right)
  df_quad$xd[1] <- df_quad$xm[1] + (df_quad$x1[1] - df_quad$xm[1])*delta
  
  df_quad$log_y0[1] <- f1(df_quad$x0[1])
  df_quad$log_y1[1] <- f1(df_quad$x1[1])
  df_quad$log_ym[1] <- f1(df_quad$xm[1])
  df_quad$log_yd[1] <- f1(df_quad$xd[1])
  
  # no need to compute all quantities for first interval, as this will
  # definitely be picked first
  df_quad$total_diff[1] <- 1
  
  # loop through remaining n_intervals
  for (i in 2:n_intervals) {
    
    # find which interval contains largest diff
    w <- which.max(df_quad$total_diff)
    
    # create new entry from xm to x1
    df_quad$x0[i] <- df_quad$xm[w]
    df_quad$x1[i] <- df_quad$x1[w]
    df_quad$xm[i] <- 0.5*(df_quad$x0[i] + df_quad$x1[i])
    df_quad$xd[i] <- df_quad$xm[i] + (df_quad$x1[i] - df_quad$xm[i])*delta
    
    df_quad$log_y0[i] <- df_quad$log_ym[w]
    df_quad$log_y1[i] <- df_quad$log_y1[w]
    df_quad$log_ym[i] <- f1(df_quad$xm[i])
    df_quad$log_yd[i] <- f1(df_quad$xd[i])
    
    df_quad$log_area_trap[i] <- get_trap_area_double(df_quad$x0[i], df_quad$xm[i], df_quad$x1[i], df_quad$log_y0[i], df_quad$log_ym[i], df_quad$log_y1[i])
    df_quad$log_area_Simp[i] <- get_Simp_area(df_quad$x0[i], df_quad$x1[i], df_quad$log_y0[i], df_quad$log_ym[i], df_quad$log_y1[i])
    df_quad$log_area_diff[i] <- get_log_area_diff(df_quad$log_area_trap[i], df_quad$log_area_Simp[i])
    df_quad$grad_diff[i] <- get_grad_diff(df_quad$x0[i], df_quad$x1[i], df_quad$log_y0[i], df_quad$log_y1[i], df_quad$log_ym[i], df_quad$log_yd[i], delta)
    df_quad$total_diff[i] <- get_total_diff(df_quad$log_area_diff[i], df_quad$grad_diff[i])
    
    # modify existing entry to go from x0 to xm
    df_quad$x1[w] = df_quad$xm[w]
    df_quad$xm[w] = 0.5*(df_quad$x0[w] + df_quad$x1[w])
    df_quad$xd[w] <- df_quad$xm[w] + (df_quad$x1[w] - df_quad$xm[w])*delta
    
    df_quad$log_y1[w] <- df_quad$log_ym[w]
    df_quad$log_ym[w] <- f1(df_quad$xm[w])
    df_quad$log_yd[w] <- f1(df_quad$xd[w])
    
    df_quad$log_area_trap[w] <- get_trap_area_double(df_quad$x0[w], df_quad$xm[w], df_quad$x1[w], df_quad$log_y0[w], df_quad$log_ym[w], df_quad$log_y1[w])
    df_quad$log_area_Simp[w] <- get_Simp_area(df_quad$x0[w], df_quad$x1[w], df_quad$log_y0[w], df_quad$log_ym[w], df_quad$log_y1[w])
    df_quad$log_area_diff[w] <- get_log_area_diff(df_quad$log_area_trap[w], df_quad$log_area_Simp[w])
    df_quad$grad_diff[w] <- get_grad_diff(df_quad$x0[w], df_quad$x1[w], df_quad$log_y0[w], df_quad$log_y1[w], df_quad$log_ym[w], df_quad$log_yd[w], delta)
    df_quad$total_diff[w] <- get_total_diff(df_quad$log_area_diff[w], df_quad$grad_diff[w])
  }
  
  # debug distribution
  if (debug_on) {
    plot_quadrature(f1, df_quad)
  }
  
  return(df_quad)
}

#------------------------------------------------
#' @title Normalise an adaptive quadrature data.frame
#'
#' @description Takes data.frame output from the function
#'   \code{adaptive_quadrature()} and normalises the y-values so the curve
#'   integrates to 1 via Simpson's rule. Also reorders rows in terms of
#'   x-position, drops some columns, and appends some columns. The final columns
#'   consist of:
#'   \itemize{
#'     \item The x and log(y) positions of all intervals.
#'     \item The log(area) via Simpson's rule.
#'     \item The Simpson's rule coefficients of the quadratic Ax^2 + Bx + C for
#'     each interval.
#'   }
#'
#' @param df_quad a data.frame output from the function
#'   \code{adaptive_quadrature()}.
#'  
#' @importFrom dplyr mutate select arrange
#'  
# do not document to keep package help simple
# @export
#' @noRd

normalise_quadrature <- function(df_quad) {
  
  # avoid "no visible binding" messages
  x0 <- xm <- x1 <- A <- B <- C <- log_y0 <- log_ym <- log_y1 <- log_area_Simp <- NULL
  
  # check inputs
  assert_dataframe(df_quad)
  assert_in(c("x0", "xm", "x1", "log_y0", "log_ym", "log_y1", "log_area_Simp"), names(df_quad))
  
  # normalise by total area
  mx <- max(df_quad$log_area_Simp)
  log_area_sum <- mx + log(sum(exp(df_quad$log_area_Simp - mx)))
  ret <- df_quad %>%
    mutate(log_y0 = log_y0 - log_area_sum,
           log_ym = log_ym - log_area_sum,
           log_y1 = log_y1 - log_area_sum,
           log_area_Simp = log_area_Simp - log_area_sum)
  
  # calculate Simpson's rule coefficients, subset columns and reorder rows
  ret <- ret %>%
    dplyr::mutate(A = get_Simp_A(x0, xm, x1, log_y0, log_ym, log_y1),
                  B = get_Simp_B(x0, xm, log_y0, log_ym, A),
                  C = get_Simp_C(x0, log_y0, A, B)) %>%
    dplyr::select(x0, x1, log_y0, log_y1, log_area_Simp, A, B, C) %>%
    dplyr::arrange(x0)
  
  return(ret)
}

#------------------------------------------------
#' @title Plot output of adaptive quadrature
#'
#' @description Takes data.frame output from the function
#'   \code{adaptive_quadrature()}, normalises through the function
#'   \code{normalise_quadrature()}, and produces a plot comparing the
#'   appximation with the true function evaluated over a grid.
#'
#' @param f1 the function used in adaptive quadrature. This is needed to compare
#'   against the adaptive quadrature result by running on a grid.
#' @param df_quad a data.frame output from the function
#'   \code{adaptive_quadrature()}.
#'  
# do not document to keep package help simple
# @export
#' @noRd

plot_quadrature <- function(f1, df_quad) {
  
  # normalise data.frame
  df_norm <- normalise_quadrature(df_quad)
  
  # work out total interval range
  left <- min(df_norm$x0)
  right <- max(df_norm$x1)
  
  # create a grid over this interval and evaluate the function by brute force
  x <- seq(left, right, l = 501)
  lx <- rep(NA, length(x))
  for (i in seq_along(x)) {
    lx[i] <- f1(x[i])
  }
  
  # normalise brute force solution
  mx <- max(df_quad$log_area_Simp)
  log_area_sum <- mx + log(sum(exp(df_quad$log_area_Simp - mx)))
  fx <- exp(lx - log_area_sum)
  
  # calculate interpolated curve from Simpson's rule
  w <- findInterval(x, vec = df_norm$x0)
  fx2 <- df_norm$A[w]*x^2 + df_norm$B[w]*x + df_norm$C[w]
  
  # plot curves
  plot(x, fx, type = 'l', lwd = 2)
  lines(x, fx2, col = 4, lwd = 2, lty = 2)
  points(df_norm$x0, exp(df_norm$log_y0), pch = 20, cex = 0.75, col = 4)
  
}

#------------------------------------------------
#' @title Integrate adaptive quadrature object up to p
#'
#' @description Returns the area under the curve of an adaptive quadrature
#'   object up to x = p (similar to e.g. \code{pnorm()}).
#'
#' @param df_norm a normalised quadrature data.frame as produced by
#'   \code{normalise_quadrature()}.
#' @param p the target area.
#'  
# do not document to keep package help simple
# @export
#' @noRd

pquad <- function(df_norm, p = 0.5) {
  
  # check inputs
  assert_dataframe(df_norm)
  assert_in(c("x0", "x1", "log_y0", "log_y1", "log_area_Simp", "A", "B", "C"), names(df_norm))
  assert_single_numeric(p)
  
  # get cumulative area
  area <- exp(df_norm$log_area_Simp)
  area_cs <- cumsum(area)
  
  # find the first interval that goes over p
  w <- which(df_norm$x1 > p)[1]
  
  # integrate this interval up to p
  area_int <- get_Simp_integral(df_norm$x0[w], p, df_norm$A[w], df_norm$B[w], df_norm$C[w])
  
  # get the total area
  ret <- area_cs[w] - area[w] + area_int
  
  return(ret)
}

#------------------------------------------------
#' @title Get quantile from adaptive quadrature object
#'
#' @description Returns the x-value for an adaptive quadrature object at which
#'   the integral up to this value equals the given quantile (similar to e.g.
#'   \code{qnorm()}).
#'
#' @param df_norm a normalised quadrature data.frame as produced by
#'   \code{normalise_quadrature()}.
#' @param q the target quantile.
#'  
# do not document to keep package help simple
# @export
#' @noRd

qquad <- function(df_norm, q = 0.5) {
  
  # check inputs
  assert_dataframe(df_norm)
  assert_in(c("x0", "x1", "log_y0", "log_y1", "log_area_Simp", "A", "B", "C"), names(df_norm))
  assert_single_bounded(q)
  
  # special case if q=0 or q=1
  if (q == 0) {
    return(min(df_norm$x0))
  } else if (q == 1) {
    return(max(df_norm$x1))
  }
  
  # get cumulative area
  area <- exp(df_norm$log_area_Simp)
  area_cs <- cumsum(area)
  
  # find the first interval that pushes over the quantile limit
  w <- which(area_cs >= q)[1]
  
  # get the total area remaining that is unaccounted for at the start of this
  # interval
  area_remaining <- q - (area_cs[w] - area[w])
  
  # solve integral for this remaining area
  ret <- solve_Simp_area(x0 = df_norm$x0[w],
                         x1 = df_norm$x1[w],
                         A = df_norm$A[w],
                         B = df_norm$B[w],
                         C = df_norm$C[w],
                         target_area = area_remaining)
  
  return(ret)
}

#------------------------------------------------
#' @title Get High Density Interval from adaptive quadrature object
#'
#' @description The HDI for a given alpha is defined as the narrowest interval
#'   in which the area under curve equals 1 - alpha. This can be obtained by
#'   searching down in the y-axis until the integrated area equals 1 - alpha.
#'   NB, although the HDI can be defined for multi-modal distributions,
#'   sometimes resulting in non-contiguous intervals, this function only returns
#'   a single interval and will throw an error if multiple intervals are
#'   detected. See details for the method used.
#'
#' @details The method used to calculate the HDI is a simple grid-based
#'   approach. The entire range implied by the quadrature object is divided up
#'   into a grid of equally sizes slices. Simpson's rule is evaluated at all of
#'   these points using the apropriate coefficients, and the size of each slice
#'   is calculated. Slices are ordered in terms of decreasing y-value, and the
#'   cumulative area of slices is calculated until we reach the desired area (1
#'   - alpha). The x-range implied by these slices is calculated, and a check is
#'   performed to ensure this results in a single contiguous interval.
#'
#' @param df_norm a normalised quadrature data.frame as produced by
#'   \code{normalise_quadrature()}.
#' @param alpha the significance level of the interval. Equivalently, one minus
#'   the area under the curve of the interval.
#' @param n_grid the number of intervals the final curve is broken into when
#'   computing the HDI on a grid.
#'  
# do not document to keep package help simple
# @export
#' @noRd

get_HDI <- function(df_norm, alpha = 0.05, n_grid = 1e3) {
  
  # avoid "no visible bindings" message
  x0 <- x1 <- x_mid <- w_interval <- A <- B <- C <- y_mid <- y_area <- cum_area <- NULL
  
  # check inputs
  assert_dataframe(df_norm)
  assert_in(c("x0", "x1", "log_y0", "log_y1", "log_area_Simp", "A", "B", "C"), names(df_norm))
  assert_single_bounded(alpha)
  
  # get full distribution range and define grid over this range
  x_min <- min(df_norm$x0)
  x_max <- max(df_norm$x1)
  x_vec <- seq(x_min, x_max, l = n_grid)
  
  # make data.frame holding the y-value and area of each slice in grid
  df_grid <- data.frame(x0 = x_vec[-length(x_vec)],
                        x1 = x_vec[-1]) %>%
    dplyr::mutate(x_mid = (x0 + x1) / 2,
                  w_interval = findInterval(x = x_mid, vec = c(df_norm$x0, x_max), rightmost.closed = TRUE),
                  A = df_norm$A[w_interval],
                  B = df_norm$B[w_interval],
                  C = df_norm$C[w_interval],
                  y_mid = A*x_mid^2 + B*x_mid + C,
                  y_area = 1/3*A*(x1^3 - x0^3) + 1/2*B*(x1^2 - x0^2) + C*(x1 - x0))
  
  # subset to HDI
  df_subset <- df_grid %>%
    dplyr::arrange(dplyr::desc(y_mid)) %>%
    dplyr::mutate(cum_area = cumsum(y_area)) %>%
    dplyr::arrange(x0) %>%
    dplyr::filter(cum_area <= (1 - alpha)) %>%
    dplyr::select(x0, x1)
  
  # assume that interval is contiguous when returning interval
  return(c(lower = min(df_subset$x0),
           upper = max(df_subset$x1)))
}

#------------------------------------------------
#' @title Get maximum of interpolated curve
#'
#' @description Returns the x-value at which y reaches its maximum value under
#'   an adaptive quadrature object.
#'
#' @param df_norm a normalised quadrature data.frame as produced by
#'   \code{normalise_quadrature()}.
#'  
# do not document to keep package help simple
# @export
#' @noRd

get_max_x <- function(df_norm) {
  
  # avoid "no visible binding" note
  A <- B <- C <- max_x <- max_y <- x0 <- x1 <- NULL
  
  ret <- df_norm %>%
    dplyr::mutate(turn_x = -0.5*B / A,
                  turn_x = ifelse(turn_x > x0, turn_x, NA),
                  turn_x = ifelse(turn_x < x1, turn_x, NA),
                  turn_y = ifelse(is.na(turn_x), NA, get_Simp_y(turn_x, A, B, C)),
                  max_y = ifelse(is.na(turn_y),
                                 ifelse(log_y0 > log_y1, exp(log_y0), exp(log_y1)),
                                 ifelse( (turn_y > exp(log_y0)) & (turn_y > exp(log_y1)),
                                         turn_y,
                                         ifelse(log_y0 > log_y1, exp(log_y0), exp(log_y1) ))),
                  max_x = ifelse(is.na(turn_y),
                                 ifelse(log_y0 > log_y1, x0, x1),
                                 ifelse( (turn_y > exp(log_y0)) & (turn_y > exp(log_y1)),
                                         turn_x,
                                         ifelse(log_y0 > log_y1, x0, x1 )))
                  ) %>%
    dplyr::filter(max_y == max(max_y)) %>%
    dplyr::pull(max_x)
  
  return(ret[1])
}
