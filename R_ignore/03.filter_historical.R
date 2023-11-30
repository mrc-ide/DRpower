# 03.filter_historical.R
#
# Author: Bob Verity
# Date: 2023-11-29
#
# Inputs:
# R_ignore/data/MTM_PFHRP23_GENE_DELETIONS_20231127_edited.xlsx
# R_ignore/data/continent_key.csv
#
# Outputs:
# R_ignore/outputs/01.dat_premap.rds
#
# Purpose:
# Starts with data downloaded from WHO malaria threats map on 27 Nov 2023. Note
# that extra columns were added to this raw dataset to specify certain rows that
# contain issues or data entry mistakes so these can be filtered out. Carries
# out a further initial filtering step based on e.g. patient type and study
# type, but does not yet filter based on number of sites as this requires
# mapping all sites to their corresponding ADMIN1 unit, which occurs in the next
# script. The intermediate object is saved to file.
#
# ------------------------------------------------------------------

library(tidyverse)
library(openxlsx)

# read in raw data
dat_raw <- openxlsx::read.xlsx("R_ignore/data/MTM_PFHRP23_GENE_DELETIONS_20231127_edited.xlsx", sheet = 2)

# drop rows that are pre-filtered due to issues with raw data
dat <- dat_raw %>%
  filter(!discard)

# filter by continent, symptomatic patients, and survey type
continent_key <- read.csv("R_ignore/data/continent_key.csv")
dat <- dat %>%
  left_join(continent_key) %>%
  filter(CONTINENT_NAME %in% c("Africa", "Asia", "South America")) %>%
  filter(PATIENT_TYPE == "Symptomatic") %>%
  filter(SURVEY_TYPE %in% c("Convenience survey", "Cross-sectional prospective survey"))

# clean up columns
dat <- dat %>%
  mutate(LATITUDE = as.numeric(LATITUDE),
         LONGITUDE = as.numeric(LONGITUDE),
         HRP2_TESTED = as.numeric(HRP2_TESTED),
         HRP2_PROPORTION_DELETION = as.numeric(HRP2_PROPORTION_DELETION),
         HRP2_NUM_DELETION = round(HRP2_TESTED * HRP2_PROPORTION_DELETION)) %>%
  select(CONTINENT_NAME, COUNTRY_NAME, SITE_NAME, LONGITUDE, LATITUDE, YEAR_START,
         YEAR_END, HRP2_TESTED, HRP2_NUM_DELETION, CITATION_URL)

# combine data collected in exact same location (lat/lon) in same year
dat <- dat %>%
  group_by(CONTINENT_NAME, COUNTRY_NAME, SITE_NAME, LONGITUDE, LATITUDE, YEAR_START,
           YEAR_END, CITATION_URL) %>%
  summarise(HRP2_TESTED = sum(HRP2_TESTED),
            HRP2_NUM_DELETION = sum(HRP2_NUM_DELETION)) %>%
  ungroup()

# must have at least 10 samples per site
dat <- dat %>%
  filter(HRP2_TESTED >= 10)

# must have at least three sites per country, collected in the same year. NB,
# this will be superceded by a filter at a lower administrative level, but
# discarding at this stage reduces the number of rows for which we need to map
# to this lower level
dat_country_year <- dat %>%
  group_by(COUNTRY_NAME, YEAR_START) %>%
  summarise(NSITES_COUNTRY = n()) %>%
  ungroup()

dat <- dat %>%
  left_join(dat_country_year) %>%
  filter(NSITES_COUNTRY >= 3) %>%
  select(-NSITES_COUNTRY)

# write to file in current state (cleaned but not mapped)
saveRDS(dat, file = "R_ignore/outputs/01.dat_premap.rds")

