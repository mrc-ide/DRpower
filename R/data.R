#------------------------------------------------
#' Summary of simulations from the threshold analysis
#' 
#' TODO
#'
#' @docType data
#'
#' @usage data(df_sim)
#'
#' @format A data frame of 547200 rows and 32 columns. The first 13 columns give
#'   parameter combinations that were used in simulating and analysing data. The
#'   "reps" column gives the number of times simulation was repeated, and "seed"
#'   gives the value of the seed that was used at the start of this loop (to
#'   ensure reproducibility). "prev_thresh" gives the prevalence threshold used
#'   in hypothesis testing. The remaining 16 columns give summary results over
#'   simulations. The MAP, posterior mean, posterior median, lower and upper
#'   credible intervals, and the probability of being above the target threshold
#'   are all summarised in terms of their mean and variance. "n_reject" gives
#'   the number of times the hypothesis of being below the threshold was
#'   rejected. This is used to estimate empirical power as \code{n_reject /
#'   reps}, and lower and upper CIs are calculated around this using the method
#'   of Clopper and Pearson.
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
#' TODO
#'
#' @docType data
#'
#' @usage data(df_ss)
#'
#' @format A data frame of 6840 rows and 15 columns. The first 14 columns give
#'   parameter combinations that were used in simulating and analysing data. The
#'   final "N_opt" column gives the optimal sample size to achieve a power of
#'   80%.
#'
#' @keywords datasets
#'
#' @examples
#' data(df_ss)
#' 
"df_ss"

#------------------------------------------------
#' TODO
#' 
#' TODO
#'
#' @docType data
#'
#' @usage data(historical_data)
#'
#' @format TODO
#'
#' @keywords datasets
#'
#' @examples
#' data(historical_data)
#' 
"historical_data"

#------------------------------------------------
#' TODO
#' 
#' TODO
#'
#' @docType data
#'
#' @usage data(studies_inclusion)
#'
#' @format TODO
#'
#' @keywords datasets
#'
#' @examples
#' data(studies_inclusion)
#' 
"studies_inclusion"
