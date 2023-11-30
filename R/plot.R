
#------------------------------------------------
##' @name plot_posterior
##' @rdname plot_posterior
##'
#' @title Plot posterior distribution of prevalence and ICC
#'
#' @description These two functions run \code{get_prevalence()} and
#'   \code{get_ICC()} respectively to obtain the full posterior distribution of
#'   the parameter of interest. Then they plot the posterior density along with
#'   some useful visualisations including the 95% CrI.
#'
#' @param prev_thresh the prevalence threshold that we are testing against
#'   (single value only, proportion between 0 and 1).
NULL

##' @rdname plot_posterior
##' 
#' @inheritParams get_prevalence
#' @param prev_range the range of prevalence values explored. Vector of two
#'   values giving lower and upper limits, defined between 0 and 1.
#'   
#' @examples
#' plot_prevalence(n = c(5, 2, 9), N = c(100, 80, 120))
#' 
#' @export

plot_prevalence <- function(n, N, prev_range = c(0, 1), alpha = 0.05, prev_thresh = 0.05,
                            prior_prev_shape1 = 1.0, prior_prev_shape2 = 1.0,
                            prior_ICC_shape1 = 1.0, prior_ICC_shape2 = 9.0,
                            CrI_type = "HDI", n_intervals = 20,
                            use_cpp = TRUE) {
  
  # avoid "no visible bindings" note
  above <- NULL
  
  # check inputs
  assert_vector_pos_int(n)
  assert_vector_pos_int(N)
  if (length(N) == 1) {
    N <- rep(N, length(n))
  }
  assert_same_length(n, N, message = "N must be either a single value (all clusters the same size) or a vector with the same length as n")
  assert_limit(prev_range)
  assert_bounded(prev_range)
  assert_single_bounded(alpha, inclusive_left = FALSE, inclusive_right = FALSE)
  assert_single_bounded(prev_thresh)
  assert_single_bounded(prior_prev_shape1, left = 1e-3, right = 1e3)
  assert_single_bounded(prior_prev_shape2, left = 1e-3, right = 1e3)
  assert_single_bounded(prior_ICC_shape1, left = 1e-3, right = 1e3)
  assert_single_bounded(prior_ICC_shape2, left = 1e-3, right = 1e3)
  assert_single_string(CrI_type)
  assert_in(CrI_type, c("HDI", "ETI"))
  assert_pos_int(n_intervals)
  assert_greq(n_intervals, 5)
  assert_single_logical(use_cpp)
  
  # get posterior
  x <- seq(prev_range[1], prev_range[2], 0.001)
  post <- get_prevalence(n = n, 
                         N = N,
                         alpha = alpha,
                         prev_thresh = prev_thresh,
                         prior_prev_shape1 = prior_prev_shape1,
                         prior_prev_shape2 = prior_prev_shape2,
                         prior_ICC_shape1 = prior_ICC_shape1,
                         prior_ICC_shape2 = prior_ICC_shape2,
                         MAP_on = TRUE,
                         post_mean_on = FALSE,
                         post_median_on = FALSE,
                         post_CrI_on = TRUE,
                         post_thresh_on = TRUE,
                         post_full_on = TRUE,
                         post_full_breaks = x,
                         CrI_type = CrI_type,
                         n_intervals = n_intervals,
                         round_digits = 2,
                         use_cpp = use_cpp)
  
  y <- post$post_full[[1]]
  
  # define plotting labels
  lab1 <- sprintf("%s%% chance below threshold", round(1e2*(1 - post$prob_above_threshold), 1))
  lab2 <- sprintf("%s%% chance above threshold", round(1e2*post$prob_above_threshold, 1))
  
  data.frame(x = x, y = y) %>%
    mutate(above = ifelse(x > 0.05, lab2, lab1),
           above = factor(above, levels = c(lab1, lab2))) %>%
    ggplot() + theme_bw() +
    geom_ribbon(aes(x = 1e2*x, ymin = 0, ymax = y, fill = above)) +
    geom_line(aes(x = 1e2*x, y = y)) +
    geom_segment(aes(x = 5, xend = 5, y = 0, yend = y[x == 0.05])) +
    geom_errorbar(aes(xmin = post$CrI_lower, xmax = post$CrI_upper, y = 1.1*max(y)), width = 0.5) +
    annotate(geom = "text", x = post$MAP, y = 1.2*max(y), label = "95% Credible Interval", hjust = 0) +
    geom_point(aes(x = post$MAP, y = 1.1*max(y))) +
    scale_fill_manual(values = c("grey", "tomato1"), name = NULL) +
    scale_x_continuous(limits = 1e2*prev_range, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.3*max(y)), expand = c(0, 0)) +
    xlab("Prevalence of pfhrp2/3 deletions") + ylab("Posterior probability density") +
    theme(legend.position = "bottom")
}

##' @rdname plot_posterior
##' 
#' @inheritParams get_prevalence
#' @param ICC_range the range of ICC values explored. Vector of two values
#'   giving lower and upper limits, defined between 0 and 1.
#'   
#' @examples
#' plot_ICC(n = c(5, 2, 9), N = c(100, 80, 120))
#' 
#' @export

plot_ICC <- function(n, N, ICC_range = c(0, 1), alpha = 0.05, prev_thresh = 0.05,
                     prior_prev_shape1 = 1.0, prior_prev_shape2 = 1.0,
                     prior_ICC_shape1 = 1.0, prior_ICC_shape2 = 9.0,
                     CrI_type = "HDI", n_intervals = 20,
                     use_cpp = TRUE) {
  
  # check inputs
  assert_vector_pos_int(n)
  assert_vector_pos_int(N)
  if (length(N) == 1) {
    N <- rep(N, length(n))
  }
  assert_same_length(n, N, message = "N must be either a single value (all clusters the same size) or a vector with the same length as n")
  assert_limit(ICC_range)
  assert_bounded(ICC_range)
  assert_single_bounded(alpha, inclusive_left = FALSE, inclusive_right = FALSE)
  assert_single_bounded(prior_prev_shape1, left = 1e-3, right = 1e3)
  assert_single_bounded(prior_prev_shape2, left = 1e-3, right = 1e3)
  assert_single_bounded(prior_ICC_shape1, left = 1e-3, right = 1e3)
  assert_single_bounded(prior_ICC_shape2, left = 1e-3, right = 1e3)
  assert_single_string(CrI_type)
  assert_in(CrI_type, c("HDI", "ETI"))
  assert_pos_int(n_intervals)
  assert_greq(n_intervals, 5)
  assert_single_logical(use_cpp)
  
  # get posterior
  x <- seq(ICC_range[1], ICC_range[2], 0.001)
  post <- get_ICC(n = n, 
                  N = N,
                  alpha = alpha,
                  prior_prev_shape1 = prior_prev_shape1,
                  prior_prev_shape2 = prior_prev_shape2,
                  prior_ICC_shape1 = prior_ICC_shape1,
                  prior_ICC_shape2 = prior_ICC_shape2,
                  MAP_on = TRUE,
                  post_mean_on = FALSE,
                  post_median_on = FALSE,
                  post_CrI_on = TRUE,
                  post_full_on = TRUE,
                  post_full_breaks = x,
                  CrI_type = CrI_type,
                  n_intervals = n_intervals,
                  round_digits = 2,
                  use_cpp = use_cpp)
  
  y <- post$post_full[[1]]
  
  data.frame(x = x, y = y) %>%
    ggplot() + theme_bw() +
    geom_ribbon(aes(x = x, ymin = 0, ymax = y), fill = "cornflowerblue") +
    geom_line(aes(x = x, y = y)) +
    geom_errorbar(aes(xmin = post$CrI_lower, xmax = post$CrI_upper, y = 1.1*max(y)), width = 0.5) +
    annotate(geom = "text", x = max(0.01, post$MAP), y = 1.2*max(y), label = "95% Credible Interval", hjust = 0) +
    geom_point(aes(x = post$MAP, y = 1.1*max(y))) +
    scale_x_continuous(limits = ICC_range, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.3*max(y)), expand = c(0, 0)) +
    xlab("Intra-cluster correlation") + ylab("Posterior probability density")
}

#------------------------------------------------
#' @title Contour plot of joint posterior distribution of prevalence and ICC
#'
#' @description Runs \code{get_joint()} to obtain the joint posterior
#'   distribution of the prevalence and the ICC. Creates a ggplot contour plot
#'   object from this result.
#'
#' @inheritParams get_joint
#' @param n_bins the number of equally spaced breaks in the contour plot. For
#'   example, 5 bins creates 4 contour lines at 20\%, 40\%, 60\% and 80\% of
#'   the maximum value.
#'
#' @examples
#' plot_joint(n = c(5, 2, 9), N = c(100, 80, 120))
#' 
#' @import ggplot2
#' @export

plot_joint <- function(n, N, 
                       prior_prev_shape1 = 1.0, prior_prev_shape2 = 1.0,
                       prior_ICC_shape1 = 1.0, prior_ICC_shape2 = 9.0,
                       prev_breaks = seq(0, 1, 0.01),
                       ICC_breaks = seq(0, 1, 0.01),
                       n_bins = 5) {
  
  # avoid "no visible bindings" note
  above <- prev <- ICC <- NULL
  
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
  assert_single_pos_int(n_bins, zero_allowed = FALSE)
  
  # get joint posterior
  z <- get_joint(n = n, N = N, 
                 prior_prev_shape1 = prior_prev_shape1,
                 prior_prev_shape2 = prior_prev_shape2,
                 prior_ICC_shape1 = prior_ICC_shape1,
                 prior_ICC_shape2 = prior_ICC_shape2,
                 prev_breaks = prev_breaks,
                 ICC_breaks = ICC_breaks)
  
  # get into long format
  df_plot <- data.frame(prev = rep(prev_breaks, each = length(ICC_breaks)),
                        ICC = ICC_breaks,
                        z = as.vector(z))
  
  # produce contour plot
  df_plot %>%
    ggplot(aes(x = 1e2*prev, y = ICC, z = z)) + theme_bw() +
    stat_contour(bins = n_bins, colour = "black") +
    scale_x_continuous(limits = 1e2*range(prev_breaks), expand = c(0, 0)) +
    scale_y_continuous(limits = range(ICC_breaks), expand = c(0, 0)) +
    xlab("Prevalence of pfhrp2/3 deletions (%)") + ylab("Intra-cluster correlation") +
    theme(plot.margin = unit(c(.2,.5,.2,.2),"cm"))
}

#------------------------------------------------
#' @title Plot a power curve using pre-computed values
#'
#' @description Runs \code{get_joint()} to obtain the joint posterior
#'   distribution of the prevalence and the ICC. Creates a ggplot contour plot
#'   object from this result.
#'
#' @param n_clust the number of clusters. Allowed values are anything in
#'   \code{2:20}, including vectors of values.
#' @param prevalence the assumed true prevalence of pfhrp2/3 deletions in the
#'   domain. Allowed values are anything in \code{seq(0, 0.2, 0.01)}, including
#'   vectors of values.
#' @param ICC the assumed intra-cluster correlation. Allowed values are" \{0,
#'   0.01, 0.02, 0.05, 0.1, 0.2\}.
#' @param prev_thresh the prevalence threshold against which we are comparing.
#'   Allowed values are: \{0.05, 0.08, 0.1\}.
#' @param N_min,N_max plotting limits on the x-axis.
#'
#' @examples
#' plot_power()
#' 
#' @import ggplot2
#' @importFrom dplyr filter near
#' @importFrom utils data
#' @export

plot_power <- function(n_clust = 5, prevalence = 0.1, ICC = 0.05,
                       prev_thresh = 0.05, N_min = 1, N_max = 2000) {
  
  # avoid "no visible binding" note
  df_sim <- N <- power <- lower <- upper <- group <- NULL
  
  # check inputs
  assert_vector_pos_int(n_clust, zero_allowed = FALSE)
  assert_vector_bounded(prevalence)
  assert_vector_bounded(ICC)
  assert_vector_bounded(prev_thresh)
  assert_single_pos_int(N_min)
  assert_single_pos_int(N_max)
  
  # load data
  data(df_sim, envir = environment())
  
  # need to rename as filter does not work when column name equals filter
  # variable
  c <- n_clust
  p <- prevalence
  r <- ICC
  p_thresh <- prev_thresh
  
  # filter data
  df_plot <- df_sim %>%
    filter(mapply(function(x) any(near(x, c)), n_clust)) %>%
    filter(mapply(function(x) any(near(x, p)), prevalence)) %>%
    filter(mapply(function(x) any(near(x, r)), ICC)) %>%
    filter(mapply(function(x) any(near(x, p_thresh)), prev_thresh)) %>%
    dplyr::filter(N >= N_min) %>%
    dplyr::filter(N <= N_max)
  
  # make grouping names based on inputs
  group_names <- ""
  first_name <- TRUE
  if (length(n_clust) > 1) {
    group_names <- sprintf("sites = %s", df_plot$n_clust)
    first_name <- FALSE
  }
  if (length(prevalence) > 1) {
    if (!first_name) {
      group_names <- sprintf("%s, ", group_names)
    }
    group_names <- sprintf("%sprevalence = %s", group_names, df_plot$prevalence)
    first_name <- FALSE
  }
  if (length(ICC) > 1) {
    if (!first_name) {
      group_names <- sprintf("%s, ", group_names)
    }
    group_names <- sprintf("%sICC = %s", group_names, df_plot$ICC)
    first_name <- FALSE
  }
  if (length(prev_thresh) > 1) {
    if (!first_name) {
      group_names <- sprintf("%s, ", group_names)
    }
    group_names <- sprintf("%sthreshold = %s", group_names, df_plot$prev_thresh)
    first_name <- FALSE
  }
  
  # produce plot
  ret <- df_plot %>%
    mutate(group = group_names) %>%
    ggplot() + theme_bw()
  if (first_name) {
    ret <-  ret + geom_pointrange(aes(x = N, y = power, ymin = lower, ymax = upper))
  } else {
    ret <-  ret + geom_pointrange(aes(x = N, y = power, ymin = lower, ymax = upper, col = group))
  }
  ret <- ret +
    geom_hline(yintercept = 80, linetype = "dashed") +
    xlim(c(N_min, N_max)) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    xlab("Sample size per cluster") +
    ylab("Power (%)") +
    theme(legend.title = element_blank())
  
  return(ret)
}
