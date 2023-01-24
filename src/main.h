
#include <Rcpp.h>

#include <vector>
#include "misc_v16.h"

//------------------------------------------------
// TODO
// [[Rcpp::export]]
Rcpp::NumericVector get_credible_prevalence_cpp(Rcpp::List args_params, Rcpp::List args_functions);

//------------------------------------------------
// TODO
double get_marginal_p_max(const std::vector<int> &k,
                          const std::vector<int> &m,
                          int GQ_intervals, int GQ_nodes,
                          double prior_p_shape1, double prior_p_shape2,
                          double prior_rho_shape1, double prior_rho_shape2,
                          const std::vector<double> &GQ_node_pos,
                          const std::vector<double> &GQ_node_weights,
                          std::vector<double> &dummy_vec,
                          double &final_loglike);

//------------------------------------------------
// TODO
double get_marginal_p_bound(const std::vector<int> &k, const std::vector<int> &m,
                            int GQ_intervals, int GQ_nodes,
                            double prior_p_shape1, double prior_p_shape2,
                            double prior_rho_shape1, double prior_rho_shape2,
                            const std::vector<double> &GQ_node_pos,
                            const std::vector<double> &GQ_node_weights,
                            std::vector<double> &dummy_vec,
                            double p_lower, double p_upper,
                            double target_ll);

//------------------------------------------------
// TODO
double loglike_marginal_p(const std::vector<int> &k,
                          const std::vector<int> &m,
                          double p, int n_intervals, int n_nodes,
                          double prior_p_shape1, double prior_p_shape2,
                          double prior_rho_shape1, double prior_rho_shape2,
                          const std::vector<double> &GQ_node_pos,
                          const std::vector<double> &GQ_node_weights,
                          std::vector<double> &dummy_vec);

//------------------------------------------------
// TODO
double get_rho_max(const std::vector<int> &k, const std::vector<int> &m, double p,
                   double prior_p_shape1, double prior_p_shape2,
                   double prior_rho_shape1, double prior_rho_shape2,
                   double &final_loglike);

//------------------------------------------------
// TODO
double get_rho_bound(const std::vector<int> &k, const std::vector<int> &m, double p,
                     double prior_p_shape1, double prior_p_shape2,
                     double prior_rho_shape1, double prior_rho_shape2,
                     double rho_lower, double rho_upper,
                     double target_ll);

//------------------------------------------------
// TODO
double loglike_joint(const std::vector<int> &k,
                     const std::vector<int> &m,
                     double p, double rho,
                     double prior_p_shape1, double prior_p_shape2,
                     double prior_rho_shape1, double prior_rho_shape2);
