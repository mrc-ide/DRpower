# 08.homepage_figure.R
#
# Author: Bob Verity
# Date: 2023-11-30
#
# Inputs: (none)
#
# Outputs: R_ignore/outputs/homepage_figure.png
#
# Purpose:
# Make a pretty plot to put on the landing page of the DRpower website.
#
# ------------------------------------------------------------------

library(DRpower)

# define parameters
n_clust <- 4
n <- rep(6, n_clust)
N <- rep(100, n_clust)

# produce plot
plot_prevalence(n, N, prev_range = c(0, 0.3))

# save to file
ggsave(filename = "R_ignore/outputs/homepage_figure.png", width = 7, height = 4.5)
