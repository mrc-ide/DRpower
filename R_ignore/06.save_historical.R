# 06.save_historical.R
#
# Author: Bob Verity
# Date: 2023-11-29
#
# Inputs:
# R_ignore/outputs/02.dat_map.rds
#
# Outputs:
# data/historical_data.rda
#
# Purpose:
# Reads in filtered historical data, saves back to base package data/ folder as
# .rda format.
#
# ------------------------------------------------------------------

# read in data
historical_data <- readRDS("R_ignore/outputs/02.dat_map.rds")

# write to data folder
save(historical_data, file = "data/historical_data.rda")
