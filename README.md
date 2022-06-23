
[![master checks](https://github.com/mrc-ide/DRpower/workflows/checks_master/badge.svg)](https://github.com/mrc-ide/DRpower/actions)
[![develop checks](https://github.com/mrc-ide/DRpower/workflows/checks_develop/badge.svg)](https://github.com/mrc-ide/DRpower/actions)

# DRpower

This package is designed to simplify the process of designing and running Plasmodium hrp2/3
gene deletion studies. It can be used to performi sample size and power calculations
before collecting any samples, and to analyse data once it arrives.

The current version of the package is very simple and contains only these functions:
  * `estimate_prevalence()` takes raw counts of pfhrp2 deletions in each of
  several clusters (e.g. clinics), along with total sample size per cluster, and
  returns an estimate of prevalence of pfhrp2 deletions and 95% confidence
  intervals. This analysis takes account of intra-cluster correlation (ICC)
  when estimating prevalence.
  * `estimate_power()` uses a simulation-and-analysis approach to determine the
  power (probability of correctly concluding that prevalence is above/below the
  5% threshold) given various parameters such as sample size and number of
  clusters.
  * `get_credible_ICC()` and `get_credible_prevalence()` are works in progress.
  They allow estimation of the ICC and prevalence using a Baysian statistical
  model.

The package can be installed using the following command:
```r
devtools::install_github("mrc-ide/DRpower", ref = "v1.0.0")
```

# Version History

The current active version is v1.0.0, which is the first release of this software.
