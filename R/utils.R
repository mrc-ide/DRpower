
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




