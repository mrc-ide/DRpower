# 01.simulate_power.R
#
# Author: Bob Verity
# Date: 2023-11-23
#
# Inputs: (none)
#
# Outputs:
# data/df_sim.rda: power estimates from bank of simulations over parameter combinations
#
# Purpose:
# Estimates power using get_power() over a large range of parameter
# combinations. Results are saved to the data/ folder.
#
# ------------------------------------------------------------------

# load packages
library(tidyverse)
library(parallel)

# define data.frame of all parameter combinations to explore
N_vec <- c(seq(5, 100, 5),
           seq(110, 500, 10),
           seq(550, 1000, 50),
           seq(1100, 2000, 100))
df_sim <- expand_grid(n_clust = c(2:20),
                      N = N_vec,
                      prevalence = seq(0.01, 0.20, 0.01),
                      ICC = c(0.0, 0.01, 0.02, 0.05, 0.10, 0.20),
                      alpha = 0.05,
                      prior_prev_shape1 = 1,
                      prior_prev_shape2 = 1,
                      prior_ICC_shape1 = 1,
                      prior_ICC_shape2 = 9,
                      CrI_type = "HDI",
                      n_intervals = 20,
                      round_digits = 2,
                      rejection_threshold = 0.95,
                      reps = 1e3,
                      seed = 1)

dim(df_sim)

# we can run multiple prevalence thresholds in a single analysis
prev_thresh <- c(0.05, 0.08, 0.10)

# function that does the actual simulation work
f1 <- function(df_sim) {
  
  ret <- list()
  for (i in 1:nrow(df_sim)) {
    
    v <- df_sim[i,]
    set.seed(v$seed)
    
    ret[[i]] <- DRpower::get_power(N = rep(v$N, v$n_clust),
                                   prevalence = v$prevalence,
                                   ICC = v$ICC,
                                   prev_thresh = prev_thresh,
                                   rejection_threshold = v$rejection_threshold,
                                   prior_prev_shape1 = v$prior_prev_shape1,
                                   prior_prev_shape2 = v$prior_prev_shape2,
                                   prior_ICC_shape1 = v$prior_ICC_shape1,
                                   prior_ICC_shape2 = v$prior_ICC_shape2,
                                   n_intervals = v$n_intervals,
                                   round_digits = v$round_digits,
                                   reps = v$reps,
                                   silent = TRUE)
  }
  return(ret)
}

# split the parameters data.frame into chunks
df_split <- split(df_sim, f = rep(1:nrow(df_sim), each = 100)[1:nrow(df_sim)])

# the main function can be run either in serial or in parallel
t0 <- Sys.time()
pow_sim <- mapply(f1, df_split, SIMPLIFY = FALSE)
Sys.time() - t0

#cl <- makeCluster(10)

#t0 <- Sys.time()
#pow_sim <- clusterApplyLB(cl = cl, x = df_split, fun = f1)
#Sys.time() - t0

#stopCluster(cl)

# get results into data.frame
df_pow <- pow_sim %>%
  bind_rows()

# merge with df_sim
df_sim <- expand_grid(df_sim, prev_thresh) %>%
  select(-prev_thresh) %>%
  cbind(df_pow)

# save to file
save(df_sim, file = "data/df_sim.rda")

# resave in better compression (avoids build warnings)
tools::resaveRdaFiles("data/df_sim.rda")
