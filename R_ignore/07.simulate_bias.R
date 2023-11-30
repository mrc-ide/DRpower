# 07.simulate_bias.R
#
# Author: Bob Verity
# Date: 2023-11-30
#
# Inputs: (none)
#
# Outputs: R_ignore/outputs/sim_bias.rds
#
# Purpose:
# Runs simulations for vignette exploring issues of bias. Simulates data and
# uses get_prevalence() over a range of true prevalences, for both small and
# medium sample sizes. Summarizes results with posterior mean (biased) and MAP
# (less biased). The point is to demonstrate that 1) the MAP estimate is not
# badly biased, 2) in all cases, bias goes down for larger sample sizes.
#
# ------------------------------------------------------------------

# simulate, summarise prevalence and take mean over simulations
get_bias <- function(N, prevalence, reps) {
  
  # draw data for all simulations. This allows us to group identical datasets,
  # which can save time in some situations (e.g. low prevalence where 0 counts
  # are common)
  l_n <- list()
  for (i in 1:reps) {
    l_n[[i]] <- data.frame(N = N,
                           n = rbbinom_reparam(n_clust = length(N), N = N,
                                               p = prevalence, rho = 0.05)) |>
      dplyr::arrange(N, n)
  }
  
  # group duplicates
  l_u <- unique(l_n)
  l_w <- tabulate(match(l_n, l_u))
  
  # simulate
  sim_df <- data.frame(MAP = rep(NA, length(l_u)),
                       post_mean = NA,
                       post_median = NA)
  for (i in seq_along(l_u)) {
    
    p_est <- get_prevalence(n = l_u[[i]]$n,
                            N = l_u[[i]]$N,
                            prior_prev_shape1 = 1,
                            prior_prev_shape2 = 1,
                            prior_ICC_shape1 = 1,
                            prior_ICC_shape2 = 9,
                            MAP_on = TRUE,
                            post_mean_on = TRUE,
                            post_median_on = TRUE,
                            post_CrI_on = FALSE,
                            post_thresh_on = FALSE,
                            post_full_on = FALSE,
                            n_intervals = 20,
                            round_digits = 2,
                            use_cpp = TRUE)
    
      sim_df$MAP[i] <- p_est$MAP
      sim_df$post_mean[i] <- p_est$post_mean
      sim_df$post_median[i] <- p_est$post_median
  }
  
  # unroll identical datasets
  sim_df <- sim_df[rep(seq_along(l_w), times = l_w),]
  
  # calculate mean
  ret <- colMeans(sim_df)
  
  return(ret)
}

# runs get_bias over a data.frame of parameters
f1 <- function(df_sim) {
  
  ret <- list()
  for (i in 1:nrow(df_sim)) {
    set.seed(1)
    ret[[i]] <- get_bias(N = rep(df_sim$N[i], 6),
                         prevalence = df_sim$prevalence[i],
                         reps = 1e3)
  }
  return(ret)
}

# ------------------------------------------------------------------

# load packages
library(tidyverse)
library(parallel)

# define data.frame of all parameter combinations to explore
df_sim <- expand_grid(N = c(5, 50, 100),
                      prevalence = seq(0.01, 0.20, 0.01))

dim(df_sim)

# split the parameters data.frame into chunks
df_split <- split(df_sim, f = rep(1:nrow(df_sim), each = 6)[1:nrow(df_sim)])

# the main function can be run either in serial or in parallel
t0 <- Sys.time()
pow_sim <- mapply(f1, df_split, SIMPLIFY = FALSE)
Sys.time() - t0

#cl <- makeCluster(10)

#clusterExport(cl, varlist = c("get_bias", "rbbinom_reparam", "get_prevalence"))

#t0 <- Sys.time()
#pow_sim <- clusterApplyLB(cl = cl, x = df_split, fun = f1)
#Sys.time() - t0

#stopCluster(cl)

# get results into data.frame
df_bias <- pow_sim %>%
  bind_rows()

# merge with df_sim
df_plot <- cbind(df_sim, df_bias)

# save to file, including to inst/extdata folder so available to vignettes
saveRDS(df_plot, file = "R_ignore/outputs/sim_bias.rds")
saveRDS(df_plot, file = "inst/extdata/sim_bias.rds")
