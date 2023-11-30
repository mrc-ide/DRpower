# 04.map_historical.R
#
# Author: Bob Verity
# Date: 2023-11-29
#
# Inputs:
# R_ignore/outputs/01.dat_premap.rds
# R_ignore/data/shp_combined.rds
# R_ignore/data/additional_data.csv
#
# Outputs:
# R_ignore/outputs/02.dat_map.rds
#
# Purpose:
# Starts with a pre-filtered dataset, and a shapefile of the relevant countries
# obtained from GADM version 4.1.0. All points are mapped against their ADMIN1
# unit, and additional data points are added for which we do not have the
# lat/lon coordinates but we do know these represent distinct spatial sites.
# Data is then filtered based on the number of sites in a given ADMIN1 unit. The
# final filtered dataset is saved to file for separate plotting.
#
# ------------------------------------------------------------------

library(sf)

# read in data that has been cleaned and pre-filtered but not mapped
dat <- readRDS("R_ignore/outputs/01.dat_premap.rds")

# read in sf object containins all the relevant countries
shp_combined <- readRDS("R_ignore/data/shp_combined.rds")

# make an sf version of the site coordinates
dat_sf <- data.frame(
  x = dat$LONGITUDE,
  y = dat$LATITUDE) |>
  sf::st_as_sf(coords = c("x", "y"), crs = 4326)

# return properties at the sample sites
intersections <- shp_combined[unlist(st_intersects(dat_sf, shp_combined)),]

# merge the ADMIN1 unit back with the original data
dat <- dat %>%
  mutate(COUNTRY = intersections$COUNTRY,
         ADMIN1_NAME = intersections$NAME_1) |>
  select(CONTINENT_NAME, COUNTRY_NAME, ADMIN1_NAME,SITE_NAME, LONGITUDE, LATITUDE, 
         YEAR_START, YEAR_END, HRP2_TESTED, HRP2_NUM_DELETION, CITATION_URL)

# load additional data and merge
dat_add <- read.csv("R_ignore/data/additional_data.csv")

dat <- rbind(dat, dat_add)

# filter to ADMIN1 with 3 or more sites from the same study in the same year
dat_n_site <- dat %>%
  group_by(COUNTRY_NAME, ADMIN1_NAME, YEAR_START, CITATION_URL) %>%
  summarise(n_site = n()) %>%
  ungroup()

dat <- dat |>
  left_join(dat_n_site) |>
  filter(n_site >= 3) |>
  select(-n_site)

# write to file in current state (cleaned and mapped)
saveRDS(dat, file = "R_ignore/outputs/02.dat_map.rds")
