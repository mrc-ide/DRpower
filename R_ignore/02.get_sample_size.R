# 02.get_sample_size.R
#
# Author: Bob Verity
# Date: 2023-11-27
#
# Inputs:
# df_sim.rds, produced by 01.simulate_power.R
#
# Outputs:
# data/df_ss.rda: estimated sample size for all parameter combinations
#
# Purpose:
# Takes raw simulation results, estimates sample size required to reach a given
# target power threshold using a simple interpolation approach. Find the value
# of N that crosses the threshold, and the value of N preceding this that does
# not, and do linear interpolation between them to get the estimated sample size
# at the threshold. Deal with special cases of N always being below the target
# power or always above the target power.
# 
# Some additional wrangling of final results, for example ensuring that N always
# decreases with increasing numbers of clusters (not always the case due to
# random variation).
#
# ------------------------------------------------------------------

# load packages
library(dplyr)
library(tidyverse)

# define target power
target_power <- 80

# get unique parameter groups and add placeholder for sample size
df_ss <- df_sim %>%
  group_by(n_clust, prevalence, ICC, alpha, prior_prev_shape1, prior_prev_shape2,
           prior_ICC_shape1, prior_ICC_shape2, CrI_type, n_intervals, round_digits,
           rejection_threshold, seed, prev_thresh) %>%
  summarise(group = cur_group_id()) %>%
  ungroup() %>%
  mutate(N_opt = NA)

# split data into list of data.frames for each group
tmp <- df_sim %>%
  left_join(df_ss)
df_list <- split(x = tmp, f = tmp$group)

# estimate sample size for all combinations
for (i in seq_along(df_list)) {
  message(sprintf("%s of %s", i, length(df_list)))
  
  # get estimated power and establish whether this crosses the threshold
  p <- df_list[[i]]$power
  cross_line <- (p[-1] >= target_power) & (p[-length(p)] < target_power)
  
  # simple plot of data
  if (FALSE) {
    plot(df_list[[i]]$N, p, ylim = c(0, 1e2))
    abline(h = target_power, lty = 2)
  }
  
  # estimate N_opt from simulation results
  if (any(cross_line)) {
    
    # find point at which power crosses target line
    w <- which(cross_line)[1] + 1
    
    # linear interpolation to find optimal N
    x0 <- df_list[[i]]$N[w - 1]
    x1 <- df_list[[i]]$N[w]
    y0 <- df_list[[i]]$power[w - 1]
    y1 <- df_list[[i]]$power[w]
    m <- (y1 - y0) / (x1 - x0)
    N_opt <- ceiling((target_power - y0 + m*x0) / m)
    
    df_ss$N_opt[i] <- N_opt
    
  } else {
    
    if (all(p[-length(p)] > target_power)) {
      df_ss$N_opt[i] <- 5
    } else {
      df_ss$N_opt[i] <- NA
    }
  }
}
df_ss <- df_ss %>%
  select(-group)

# checks on sample sizes:
# ideally sample size should be decreasing with n_clust. View cases where this
# is not true
df_ss %>%
  group_by_at(vars(-c(n_clust, N_opt))) %>%
  summarise(pass = all(diff(N_opt) <= 0)) %>%
  filter(!pass) %>%
  left_join(df_ss) %>%
  pivot_wider(names_from = n_clust, values_from = N_opt) %>%
  View()

# force decreasing with n_clust
f1 <- function(x) {
  m <- 0  # m is the maximum value seen so far when searching right to left
  for (i in length(x):1) {
    if (is.na(x[i])) {
      break
    }
    if (x[i] < m) {
      x[i] <- m
    }
    m <- max(m, x[i])
  }
  return(x)
}

df_ss <- df_ss %>%
  group_by_at(vars(-c(n_clust, N_opt))) %>%
  summarise(n_clust = n_clust,
            N_opt = f1(N_opt)) %>%
  ungroup()

# quick look at results
df_ss %>%
  filter(ICC == 0.05) %>%
  filter(prev_thresh == 0.05) %>%
  select(n_clust, prevalence, N_opt) %>%
  pivot_wider(names_from = prevalence, values_from = N_opt) %>%
  View()

# write to file
save(df_ss, file = "data/df_ss.rda")
