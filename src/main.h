
#include <Rcpp.h>

#include <vector>
#include "misc_v16.h"

//------------------------------------------------
double dbbinom_reparam(const std::vector<int> &n,
                       const std::vector<int> &N,
                       double p, double rho);

//------------------------------------------------
double loglike_joint(const std::vector<int> &n,
                     const std::vector<int> &N,
                     double p, double rho,
                     double prior_p_shape1, double prior_p_shape2,
                     double prior_rho_shape1, double prior_rho_shape2);

//------------------------------------------------
double loglike_marginal_p(const std::vector<int> &n,
                          const std::vector<int> &N,
                          double p, int n_intervals,
                          double prior_p_shape1, double prior_p_shape2,
                          double prior_rho_shape1, double prior_rho_shape2,
                          std::vector<double> &x0, std::vector<double> &xm,
                          std::vector<double> &x1, std::vector<double> &log_y0,
                          std::vector<double> &log_ym, std::vector<double> &log_y1,
                          std::vector<double> &log_area_trap,
                          std::vector<double> &log_area_Simp,
                          std::vector<double> &log_area_diff);

//------------------------------------------------
double loglike_marginal_rho(const std::vector<int> &n,
                            const std::vector<int> &N,
                            double rho, int n_intervals,
                            double prior_p_shape1, double prior_p_shape2,
                            double prior_rho_shape1, double prior_rho_shape2,
                            std::vector<double> &x0, std::vector<double> &xm,
                            std::vector<double> &x1, std::vector<double> &log_y0,
                            std::vector<double> &log_ym, std::vector<double> &log_y1,
                            std::vector<double> &log_area_trap,
                            std::vector<double> &log_area_Simp,
                            std::vector<double> &log_area_diff);

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List get_credible_prevalence_cpp(Rcpp::List args_params);

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List get_credible_ICC_cpp(Rcpp::List args_params);
