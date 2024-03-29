% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_sample_size_table}
\alias{get_sample_size_table}
\title{Get pre-computed sample size tables}
\usage{
get_sample_size_table(
  prevalence = seq(0, 0.2, 0.01),
  ICC = 0.05,
  prev_thresh = 0.05
)
}
\arguments{
\item{prevalence}{the assumed true prevalence of pfhrp2/3 deletions in the
domain. Allowed values are anything in \code{seq(0, 0.2, 0.01)}, including
vectors of values.}

\item{ICC}{the assumed intra-cluster correlation. Allowed values are" \{0,
0.01, 0.02, 0.05, 0.1, 0.2\}.}

\item{prev_thresh}{the prevalence threshold against which we are comparing.
Allowed values are: \{0.05, 0.08, 0.1\}.}
}
\description{
Produce a sample size table giving the minimum sample size per
  cluster for given values of the ICC and the prevalence threshold against
  which we are comparing.
}
\details{
The function \code{get_power_threshold()} was run over a large range
  of parameter combinations and results were stored within the \code{df_sim}
  object (see \code{?df_sim}). These simulations were then used to produce
  minimum sample size estimates by linear interpolation that were stored
  within the \code{df_ss} object (see \code{?df_ss}). This function provides
  a simple way of querying the \code{df_ss} object for given parameter
  values.
}
\examples{
get_sample_size_table()

}
