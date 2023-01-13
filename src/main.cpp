
#include "main.h"

using namespace std;

//------------------------------------------------
// TODO
Rcpp::NumericVector get_credible_prevalence_cpp(Rcpp::List args_params, Rcpp::List args_functions) {
  
  // extract inputs
  vector<int> n = rcpp_to_vector_int(args_params["n"]);
  vector<int> N = rcpp_to_vector_int(args_params["N"]);
  double alpha = rcpp_to_double(args_params["alpha"]);
  double prior_prev_shape1 = rcpp_to_double(args_params["prior_prev_shape1"]);
  double prior_prev_shape2 = rcpp_to_double(args_params["prior_prev_shape2"]);
  double prior_ICC_shape1 = rcpp_to_double(args_params["prior_ICC_shape1"]);
  double prior_ICC_shape2 = rcpp_to_double(args_params["prior_ICC_shape2"]);
  int GQ_intervals = rcpp_to_int(args_params["GQ_intervals"]);
  int GQ_nodes = rcpp_to_int(args_params["GQ_nodes"]);
  vector<double> GQ_node_pos = rcpp_to_vector_double(args_params["GQ_node_pos"]);
  vector<double> GQ_node_weights = rcpp_to_vector_double(args_params["GQ_node_weights"]);
  
  // extract R functions
  Rcpp::Function solve_Simpsons_area = args_functions["solve_Simpsons_area"];
  
  // define some hard-coded parameters and dummy variables
  int n_intervals = 40;
  double precision_limit = 6*log(10);
  vector<double> dummy_vec(GQ_intervals * GQ_nodes);
  
  // get maximum value of p
  double p_ml_ll = 0.0;
  double p_ml = get_marginal_p_max(n, N, GQ_intervals, GQ_nodes, precision_limit,
                                   prior_prev_shape1, prior_prev_shape2,
                                   prior_ICC_shape1, prior_ICC_shape2,
                                   GQ_node_pos, GQ_node_weights, dummy_vec,
                                   p_ml_ll);
  
  // get lower and upper bounds at which distribution is a factor
  // exp(precision_limit) smaller than the maximum we just found. We will only
  // explore this interval as it contains nearly all the probability mass
  double p_lower = get_marginal_p_bound(n, N, GQ_intervals, GQ_nodes,
                                        prior_prev_shape1, prior_prev_shape2,
                                        prior_ICC_shape1, prior_ICC_shape2,
                                        GQ_node_pos, GQ_node_weights, dummy_vec,
                                        0.0, p_ml, p_ml_ll, precision_limit);
  double p_upper = get_marginal_p_bound(n, N, GQ_intervals, GQ_nodes,
                                        prior_prev_shape1, prior_prev_shape2,
                                        prior_ICC_shape1, prior_ICC_shape2,
                                        GQ_node_pos, GQ_node_weights, dummy_vec,
                                        p_ml, 1.0, p_ml_ll, precision_limit);
  
  // split domain into intervals and use quadratic smoothing over each interval
  // to calculate area under curve. Store interval areas, and also y-values at
  // interval ends and midpoints
  vector<double> interval_area(n_intervals);
  vector<double> interval_y_ends(n_intervals + 1);
  vector<double> interval_y_mids(n_intervals);
  
  double interval_width = (p_upper - p_lower) / double(n_intervals);
  interval_y_ends[0] = exp(loglike_marginal_p(n, N, p_lower, GQ_intervals, GQ_nodes, precision_limit,
                                              prior_prev_shape1, prior_prev_shape2,
                                              prior_ICC_shape1, prior_ICC_shape2,
                                              GQ_node_pos, GQ_node_weights, dummy_vec) - p_ml_ll);
  for (int i = 0; i < n_intervals; ++i) {
    double x_mid = p_lower + (i + 0.5)*interval_width;
    interval_y_mids[i] = exp(loglike_marginal_p(n, N, x_mid, GQ_intervals, GQ_nodes, precision_limit,
                                                prior_prev_shape1, prior_prev_shape2,
                                                prior_ICC_shape1, prior_ICC_shape2,
                                                GQ_node_pos, GQ_node_weights, dummy_vec) - p_ml_ll);
    
    double x_right = p_lower + (i + 1)*interval_width;
    interval_y_ends[i+1] = exp(loglike_marginal_p(n, N, x_right, GQ_intervals, GQ_nodes, precision_limit,
                                                  prior_prev_shape1, prior_prev_shape2,
                                                  prior_ICC_shape1, prior_ICC_shape2,
                                                  GQ_node_pos, GQ_node_weights, dummy_vec) - p_ml_ll);
    
    interval_area[i] = interval_width / 6.0 * (interval_y_ends[i] + 4*interval_y_mids[i] + interval_y_ends[i+1]);
  }
  
  // get target area of integral for lower and upper CrI
  double target_area_lower = 0.5* alpha * sum(interval_area);
  double target_area_upper = (1.0 - 0.5*alpha) * sum(interval_area);
  
  // find which intervals these target areas fall within
  int w_lower = 0;
  int w_upper = 0;
  double area_remaining_lower = 0.0;
  double area_remaining_upper = 0.0;
  double area_cumsum = 0.0;
  for (int i = 0; i < n_intervals; ++i) {
    area_cumsum += interval_area[i];
    if (area_cumsum > target_area_lower) {
      w_lower = i;
      area_remaining_lower = target_area_lower - (area_cumsum - interval_area[i]);
      break;
    }
  }
  for (int i = (w_lower + 1); i < n_intervals; ++i) {
    area_cumsum += interval_area[i];
    if (area_cumsum > target_area_upper) {
      w_upper = i;
      area_remaining_upper = target_area_upper - (area_cumsum - interval_area[i]);
      break;
    }
  }
  
  // solve for lower and upper CrIs
  double CrI_lower = rcpp_to_double(solve_Simpsons_area(p_lower + w_lower*interval_width,
                                                        p_lower + (w_lower + 1)*interval_width,
                                                        interval_y_ends[w_lower],
                                                        interval_y_mids[w_lower],
                                                        interval_y_ends[w_lower + 1],
                                                        area_remaining_lower));
  
  double CrI_upper = rcpp_to_double(solve_Simpsons_area(p_lower + w_upper*interval_width,
                                                        p_lower + (w_upper + 1)*interval_width,
                                                        interval_y_ends[w_upper],
                                                        interval_y_mids[w_upper],
                                                        interval_y_ends[w_upper + 1],
                                                        area_remaining_upper));
  
  // return outputs as vector
  Rcpp::NumericVector ret = {CrI_lower, CrI_upper};
  return ret;
}

//------------------------------------------------
// TODO
double get_marginal_p_max(const vector<int> &k, const vector<int> &m,
                          int GQ_intervals, int GQ_nodes,
                          double precision_limit,
                          double prior_p_shape1, double prior_p_shape2,
                          double prior_rho_shape1, double prior_rho_shape2,
                          const vector<double> &GQ_node_pos,
                          const vector<double> &GQ_node_weights,
                          vector<double> &dummy_vec,
                          double &final_loglike) {
  
  // set parameters of optimisation
  double delta = 1e-4;
  int iterations = 10;
  
  // perform binary search (i.e. subdivide interval) based on gradient
  double p_lower = 0.0;
  double p_upper = 1.0;
  double p = 0.5;
  for (int i = 0; i < iterations; ++i) {
    double l1 = loglike_marginal_p(k, m, p, GQ_intervals, GQ_nodes, precision_limit,
                                   prior_p_shape1, prior_p_shape2,
                                   prior_rho_shape1, prior_rho_shape2,
                                   GQ_node_pos, GQ_node_weights, dummy_vec);
    double l2 = loglike_marginal_p(k, m, p + delta, GQ_intervals, GQ_nodes, precision_limit,
                                   prior_p_shape1, prior_p_shape2,
                                   prior_rho_shape1, prior_rho_shape2,
                                   GQ_node_pos, GQ_node_weights, dummy_vec);
    double grad = l2 - l1;
    if (grad > 0.0) {
      p_lower = p;
    } else {
      p_upper = p;
    }
    p = 0.5*(p_lower + p_upper);
  }
  
  // calculate loglikelihood at final value (returned by reference)
  final_loglike = loglike_marginal_p(k, m, p, GQ_intervals, GQ_nodes, precision_limit,
                                     prior_p_shape1, prior_p_shape2,
                                     prior_rho_shape1, prior_rho_shape2,
                                     GQ_node_pos, GQ_node_weights, dummy_vec);
  
  return p;
}

//------------------------------------------------
// TODO
double get_marginal_p_bound(const vector<int> &k, const vector<int> &m,
                            int GQ_intervals, int GQ_nodes,
                            double prior_p_shape1, double prior_p_shape2,
                            double prior_rho_shape1, double prior_rho_shape2,
                            const vector<double> &GQ_node_pos,
                            const vector<double> &GQ_node_weights,
                            vector<double> &dummy_vec,
                            double p_lower, double p_upper,
                            double target_ll, double precision_limit) {
  
  // set parameters of optimisation
  double delta = 1e-4;
  int iterations = 10;
  
  // perform binary search (i.e. subdivide interval) based on gradient
  double p = 0.5*(p_lower + p_upper);
  for (int i = 0; i < iterations; ++i) {
    double l1 = loglike_marginal_p(k, m, p, GQ_intervals, GQ_nodes, precision_limit,
                                   prior_p_shape1, prior_p_shape2,
                                   prior_rho_shape1, prior_rho_shape2,
                                   GQ_node_pos, GQ_node_weights, dummy_vec);
    l1 = -abs(l1 - target_ll + precision_limit);
    double l2 = loglike_marginal_p(k, m, p + delta, GQ_intervals, GQ_nodes, precision_limit,
                                   prior_p_shape1, prior_p_shape2,
                                   prior_rho_shape1, prior_rho_shape2,
                                   GQ_node_pos, GQ_node_weights, dummy_vec);
    l2 = -abs(l2 - target_ll + precision_limit);
    double grad = l2 - l1;
    if (grad > 0.0) {
      p_lower = p;
    } else {
      p_upper = p;
    }
    p = 0.5*(p_lower + p_upper);
  }
  
  return p;
}

//------------------------------------------------
// TODO
double loglike_marginal_p(const vector<int> &k,
                          const vector<int> &m,
                          double p, int n_intervals, int n_nodes,
                          double precision_limit,
                          double prior_p_shape1, double prior_p_shape2,
                          double prior_rho_shape1, double prior_rho_shape2,
                          const vector<double> &GQ_node_pos,
                          const vector<double> &GQ_node_weights,
                          vector<double> &dummy_vec) {
  
  // get maximum likelihood estimate of rho and log-likelihood at this max value
  double rho_ml_ll = 0.0;
  double rho_ml = get_rho_max(k, m, p, prior_p_shape1, prior_p_shape2,
                              prior_rho_shape1, prior_rho_shape2,
                              rho_ml_ll);
  
  // get lower and upper bounds at which distribution is a factor
  // exp(precision_limit) smaller than the maximum we just found. We will only
  // integrate over this interval as it contains nearly all the probability mass
  double rho_lower = get_rho_bound(k, m, p, prior_p_shape1, prior_p_shape2,
                                   prior_rho_shape1, prior_rho_shape2,
                                   0.0,  rho_ml, rho_ml_ll, precision_limit);
  double rho_upper = get_rho_bound(k, m, p, prior_p_shape1, prior_p_shape2,
                                   prior_rho_shape1, prior_rho_shape2,
                                   rho_ml, 1.0, rho_ml_ll, precision_limit);
  
  // calculate log-likelihood at selected points
  double interval_width = (rho_upper - rho_lower) / double(n_intervals);
  int i2 = 0;
  for (int i = 0; i < n_intervals; ++i) {
    for (int j = 0; j < n_nodes; ++j) {
      double rho = rho_lower + (i + GQ_node_pos[j])*interval_width;
      dummy_vec[i2] = loglike_joint(k, m, p, rho, prior_p_shape1, prior_p_shape2,
                                    prior_rho_shape1, prior_rho_shape2);
      i2++;
    }
  }
  
  // perform quadrature in underflow-safe way
  double ll_max = max(dummy_vec);
  double ll_sum = 0.0;
  i2 = 0;
  for (int i = 0; i < n_intervals; ++i) {
    for (int j = 0; j < n_nodes; ++j) {
      ll_sum += 0.5 * interval_width * GQ_node_weights[j] * exp(dummy_vec[i2] - ll_max);
      i2++;
    }
  }
  double ret = ll_max + log(ll_sum);
  
  return ret;
}

//------------------------------------------------
// TODO
double get_rho_max(const vector<int> &k, const vector<int> &m, double p,
                   double prior_p_shape1, double prior_p_shape2,
                   double prior_rho_shape1, double prior_rho_shape2,
                   double &final_loglike) {
  
  // set parameters of optimisation
  double delta = 1e-4;
  int iterations = 10;
  
  // perform binary search (i.e. subdivide interval) based on gradient
  double rho_lower = 0.0;
  double rho_upper = 1.0;
  double rho = 0.5;
  for (int i = 0; i < iterations; ++i) {
    double l1 = loglike_joint(k, m, p, rho, prior_p_shape1, prior_p_shape2,
                              prior_rho_shape1, prior_rho_shape2);
    double l2 = loglike_joint(k, m, p, rho + delta, prior_p_shape1, prior_p_shape2,
                              prior_rho_shape1, prior_rho_shape2);
    double grad = l2 - l1;
    if (grad > 0.0) {
      rho_lower = rho;
    } else {
      rho_upper = rho;
    }
    rho = 0.5*(rho_lower + rho_upper);
  }
  
  // calculate loglikelihood at final value (returned by reference)
  final_loglike = loglike_joint(k, m, p, rho, prior_p_shape1, prior_p_shape2,
                                prior_rho_shape1, prior_rho_shape2);
  
  return rho;
}

//------------------------------------------------
// TODO
double get_rho_bound(const vector<int> &k, const vector<int> &m, double p,
                     double prior_p_shape1, double prior_p_shape2,
                     double prior_rho_shape1, double prior_rho_shape2,
                     double rho_lower, double rho_upper,
                     double target_ll, double precision_limit) {
  
  // set parameters of optimisation
  double delta = 1e-4;
  int iterations = 10;
  
  // perform binary search (i.e. subdivide interval) based on gradient
  double rho = 0.5*(rho_lower + rho_upper);
  for (int i = 0; i < iterations; ++i) {
    double l1 = loglike_joint(k, m, p, rho, prior_p_shape1, prior_p_shape2,
                              prior_rho_shape1, prior_rho_shape2);
    l1 = -abs(l1 - target_ll + precision_limit);
    double l2 = loglike_joint(k, m, p, rho + delta, prior_p_shape1, prior_p_shape2,
                              prior_rho_shape1, prior_rho_shape2);
    l2 = -abs(l2 - target_ll + precision_limit);
    double grad = l2 - l1;
    if (grad > 0.0) {
      rho_lower = rho;
    } else {
      rho_upper = rho;
    }
    rho = 0.5*(rho_lower + rho_upper);
  }
  
  return rho;
}

//------------------------------------------------
// TODO
double loglike_joint(const vector<int> &k,
                     const vector<int> &m,
                     double p, double rho,
                     double prior_p_shape1, double prior_p_shape2,
                     double prior_rho_shape1, double prior_rho_shape2) {
  
  // calculate likelihood over all clusters
  int n_clust = k.size();
  double alpha = p*(1.0 / rho - 1.0);
  double beta = (1.0 - p)*(1.0 / rho - 1.0);
  double ret = n_clust * (lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta));
  for (int i = 0; i < n_clust; ++ i) {
    ret += lgamma(m[i] + 1) - lgamma(k[i] + 1) - lgamma(m[i] - k[i] + 1);
    ret += lgamma(k[i] + alpha) + lgamma(m[i] - k[i] + beta) - lgamma(m[i] + alpha + beta);
  }
  
  // apply priors
  if ((prior_p_shape1 != 1.0) || (prior_p_shape2 != 1.0)) {
    ret += lgamma(prior_p_shape1 + prior_p_shape2) - lgamma(prior_p_shape1) - lgamma(prior_p_shape2) +
      (prior_p_shape1 - 1)*log(p) + (prior_p_shape2 - 1)*log(1.0 - p);
  }
  if ((prior_rho_shape1 != 1.0) || (prior_rho_shape2 != 1.0)) {
    ret += lgamma(prior_rho_shape1 + prior_rho_shape2) - lgamma(prior_rho_shape1) - lgamma(prior_rho_shape2) +
      (prior_rho_shape1 - 1)*log(p) + (prior_rho_shape2 - 1)*log(1.0 - p);
  }
  
  return ret;
}
