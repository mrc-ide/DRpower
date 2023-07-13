pkgname <- "DRpower"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('DRpower')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("df_sim")
### * df_sim

flush(stderr()); flush(stdout())

### Name: df_sim
### Title: Summary of simulations from the threshold analysis
### Aliases: df_sim
### Keywords: datasets

### ** Examples

data(df_sim)




cleanEx()
nameEx("df_ss")
### * df_ss

flush(stderr()); flush(stdout())

### Name: df_ss
### Title: Minimum sample sizes for the threshold analysis
### Aliases: df_ss
### Keywords: datasets

### ** Examples

data(df_ss)




cleanEx()
nameEx("get_joint_grid")
### * get_joint_grid

flush(stderr()); flush(stdout())

### Name: get_joint_grid
### Title: Get posterior distribution of both prevalence and the ICC on a
###   grid
### Aliases: get_joint_grid

### ** Examples

get_joint_grid(n = c(5, 2, 9), N = c(100, 80, 120))




cleanEx()
nameEx("get_posterior")
### * get_posterior

flush(stderr()); flush(stdout())

### Name: get_posterior
### Title: Estimate prevalence and intra-cluster correlation from raw
###   counts
### Aliases: get_posterior get_prevalence get_ICC

### ** Examples

# basic example of estimating prevalence and ICC from observed counts
df_counts <- data.frame(sample_size = c(80, 110, 120),
                        deletions = c(3, 5, 6))
get_prevalence(n = df_counts$deletions, N = df_counts$sample_size)
get_ICC(n = df_counts$deletions, N = df_counts$sample_size)




cleanEx()
nameEx("get_power_presence")
### * get_power_presence

flush(stderr()); flush(stdout())

### Name: get_power_presence
### Title: Calculate power when testing for presence of deletions
### Aliases: get_power_presence

### ** Examples

get_power_presence(N = c(120, 90, 150), prevalence = 0.01, ICC = 0.1)




cleanEx()
nameEx("get_power_threshold")
### * get_power_threshold

flush(stderr()); flush(stdout())

### Name: get_power_threshold
### Title: Estimate power when testing prevalence against a threshold
### Aliases: get_power_threshold

### ** Examples

get_power_threshold(N = c(120, 90, 150), prevalence = 0.15, ICC = 0.1 , reps = 1e2)




cleanEx()
nameEx("get_sample_size_presence")
### * get_sample_size_presence

flush(stderr()); flush(stdout())

### Name: get_sample_size_presence
### Title: Get minimum sample size when testing for presence of deletions
### Aliases: get_sample_size_presence

### ** Examples

get_sample_size_presence(n_clust = 5, prevalence = 0.01, ICC = 0.1)




cleanEx()
nameEx("historical_data")
### * historical_data

flush(stderr()); flush(stdout())

### Name: historical_data
### Title: TODO
### Aliases: historical_data
### Keywords: datasets

### ** Examples

data(historical_data)




cleanEx()
nameEx("studies_inclusion")
### * studies_inclusion

flush(stderr()); flush(stdout())

### Name: studies_inclusion
### Title: TODO
### Aliases: studies_inclusion
### Keywords: datasets

### ** Examples

data(studies_inclusion)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
