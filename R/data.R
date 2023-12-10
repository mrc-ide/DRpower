#------------------------------------------------
#' Summary of simulations from the threshold analysis
#' 
#' @description This object was produced by running the function
#'   \code{get_power_threshold()} over a wide range of parameter combinations.
#'   This data.frame contains the results of these simulations attached to the
#'   parameter values used in simulation. The most obvious use of this object is
#'   in constructing power curves over a range of sample sizes (see
#'   \code{?power_curve()}).
#'
#' @docType data
#'
#' @usage data(df_sim)
#'
#' @format A data.frame of 547200 rows and 19 columns. The first 13 columns give
#'   parameter combinations that were used in simulating and analysing data. The
#'   "reps" column gives the number of times simulation was repeated, and "seed"
#'   gives the value of the seed that was used at the start of this loop (to
#'   ensure reproducibility). "prev_thresh" gives the prevalence threshold used
#'   in hypothesis testing. The final three columns give the estimates power
#'   over simulations along with lower and upper 95\% CIs calculated using the
#'   method of Clopper and Pearson (1934).
#'
#' @keywords datasets
#'
#' @examples
#' data(df_sim)
#' 
"df_sim"

#------------------------------------------------
#' Minimum sample sizes for the threshold analysis
#' 
#' @description This object was produced by finding the point at which
#'   \code{df_sim} crossed the target power threshold of 80\% (see details).
#'
#' @details Minimum sample sizes were calculated as follows:
#'   \enumerate{
#'   \item Find the value of N that crosses the threshold, and the value of N
#'   preceding it that does not.
#'   \item Do linear interpolation between these two values to get the estimated
#'   sample size at the threshold. \item Deal with special cases of N always
#'   being below the target power or always above the target power.
#'   \item Some additional manual wrangling of final results. Ensure that N
#'   always decreases with increasing numbers of clusters (this is not always
#'   the case due to random variation).
#'   }
#'
#' @docType data
#'
#' @usage data(df_ss)
#'
#' @format A data.frame of 6840 rows and 15 columns. The first 14 columns give
#'   parameter combinations that were used in simulating and analysing data. The
#'   final "N_opt" column gives the optimal sample size to achieve a power of
#'   80\%.
#'
#' @keywords datasets
#'
#' @examples
#' data(df_ss)
#' 
"df_ss"

#------------------------------------------------
#' Data from historical pfhrp2 studies that passed filters for inclusion into an
#' ICC analysis.
#' 
#' @description A data.frame of sites that were used to estimate the ICC based
#'   on previously published data. These sites passed strict inclusion criteria
#'   to ensure they are maximally informative (see details).
#'
#' @details The raw dataset of historical pfhrp2/3 studies was downloaded from
#'   the WHO malaria threats map on 27 Nov 2023. This spreadsheet can be found
#'   in this package in the R_ignore/data folder (see the
#'   \href{https://github.com/mrc-ide/DRpower/}{Github repos}) under the name
#'   "MTM_PFHRP23_GENE_DELETIONS_20231127_edited.xlss". Note that this
#'   spreadsheet has the term "_edited" added to the name because two extra
#'   columns were added to the original data: "discard" and "discard_reason".
#'   These columns specify certain rows that should be discarded in the original
#'   data due to data entry mistakes. The following steps were then taken. All
#'   scripts to perform these steps can be found in the same R_ignore folder:
#'   \enumerate{
#'     \item Rows were dropped that were identified to discard based on problems
#'     in the original data.
#'     \item Filtered to Africa, Asia or South America.
#'     \item Filtered to Symptomatic patients.
#'     \item Filtered to convenience surveys or cross-sectional prospective
#'     surveys only.
#'     \item Combined counts (tested and positive) of studies conducted in the
#'     same exact location (based on lat/lon) in the same year and from the same
#'     source publication. These are considered a single site.
#'     \item Filtered to have 10 or more samples per site. Avoids very small
#'     sample sizes which would have very little information from the data and
#'     therefore would be driven by our prior assumptions.
#'     \item All sites were mapped to ADMIN1 level by comparing the lat/lon
#'     coordinates against a shapefile from \href{https://gadm.org/}{GADM}
#'     version 4.1.0, first administrative unit.
#'     \item Results were combined with studies that contain additional
#'     information not reflected in the WHO malaria threats map data. For
#'     example, some studies have site-level resolution despite not being
#'     apparent in the original data download. These additional studies can be
#'     found in the \href{https://github.com/mrc-ide/DRpower/}{R_ignore/data
#'     folder} under the name "additional_data.csv".
#'     \item Filter to ADMIN1 regions that contain at least 3 distinct sites
#'     within the same year and from the same source publication.
#'   }
#'   This final filtered dataset is what is available here.
#'
#' @docType data
#'
#' @usage data(historical_data)
#'
#' @format A data.frame of 30 rows and 11 columns. Each row gives a different
#'   site that made it through filtering steps in the ICC analysis from
#'   historical data. Coluns give geographic properties, sampling times, the
#'   number of samples tested and positive for pfhrp2 deletions, and the
#'   citation from which the data originates.
#'
#' @keywords datasets
#'
#' @examples
#' data(historical_data)
#' 
"historical_data"
