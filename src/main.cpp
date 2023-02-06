
#include "main.h"
#include "adaptive_quadrature.h"

using namespace std;

//------------------------------------------------
// reparameterisation of the beta-binomial distribution in terms of a mean (p)
// and an intra-cluster correlation coefficient (rho). Deals with special cases
// that simplify to the binomial or bernoulli distributions

double dbbinom_reparam(const vector<int> &n,
                       const vector<int> &N,
                       double p, double rho) {
  
  // deal with special cases
  int n_clust = n.size();
  if (rho == 0.0) {
    // simplifies to binomial distribution
    double ret = 0.0;
    for (int i = 0; i < n_clust; ++ i) {
      ret += lgamma(N[i] + 1) - lgamma(n[i] + 1) - lgamma(N[i] - n[i] + 1) +
        n[i]*log(p) + (N[i] - n[i])*log(1.0 - p);
    }
    return ret;
  } else if (rho == 1.0) {
    // likelihood still positive whenever k == 0 or k == m
    double ret = 0.0;
    for (int i = 0; i < n_clust; ++ i) {
      if (n[i] == 0) {
        ret += log(1.0 - p);
      } else if (n[i] == N[i]) {
        ret += log(p);
      } else {
        ret = -INFINITY;
        break;
      }
    }
    return ret;
  }
  if (p == 0.0) {
    // likelihood 1 whenever k == 0
    for (int i = 0; i < n_clust; ++ i) {
      if (n[i] != 0) {
        return -INFINITY;
      }
    }
    return 0.0;
  } else if (p == 1.0) {
    // likelihood 1 whenever k == m
    for (int i = 0; i < n_clust; ++ i) {
      if (n[i] != N[i]) {
        return -INFINITY;
      }
    }
    return 0.0;
  }
  
  // calculate Beta-binomial likelihood over all clusters
  double alpha = p*(1.0 / rho - 1.0);
  double beta = (1.0 - p)*(1.0 / rho - 1.0);
  double ret = n_clust * (lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta));
  for (int i = 0; i < n_clust; ++ i) {
    ret += lgamma(N[i] + 1) - lgamma(n[i] + 1) - lgamma(N[i] - n[i] + 1) +
      lgamma(n[i] + alpha) + lgamma(N[i] - n[i] + beta) - lgamma(N[i] + alpha + beta);
  }
  
  return ret;
}

//------------------------------------------------
// joint probability of data multiplied by priors on p and rho

double loglike_joint(const vector<int> &n,
                     const vector<int> &N,
                     double p, double rho,
                     double prior_p_shape1, double prior_p_shape2,
                     double prior_rho_shape1, double prior_rho_shape2) {
  
  // calculate log-likelihood
  double ret = dbbinom_reparam(n, N, p, rho);
  
  // apply priors (skip if uniform)
  if ((prior_p_shape1 != 1.0) || (prior_p_shape2 != 1.0)) {
    ret += lgamma(prior_p_shape1 + prior_p_shape2) - lgamma(prior_p_shape1) - lgamma(prior_p_shape2);
    if (p > 0) {
      ret += (prior_p_shape1 - 1)*log(p);
    }
    if (p < 1) {
      ret += (prior_p_shape2 - 1)*log(1.0 - p);
    }
  }
  if ((prior_rho_shape1 != 1.0) || (prior_rho_shape2 != 1.0)) {
    ret += lgamma(prior_rho_shape1 + prior_rho_shape2) - lgamma(prior_rho_shape1) - lgamma(prior_rho_shape2);
    if (rho > 0) {
      ret += (prior_rho_shape1 - 1)*log(rho);
    }
    if (rho < 1) {
      ret += (prior_rho_shape2 - 1)*log(1.0 - rho);
    }
  }
  
  return ret;
}

//------------------------------------------------
// prior*likelihood integrated over rho

double loglike_marginal_p(const vector<int> &n,
                          const vector<int> &N,
                          double p, int n_intervals,
                          double prior_p_shape1, double prior_p_shape2,
                          double prior_rho_shape1, double prior_rho_shape2,
                          vector<double> &x0, vector<double> &xm,
                          vector<double> &x1, vector<double> &log_y0,
                          vector<double> &log_ym, vector<double> &log_y1,
                          vector<double> &log_area_trap,
                          vector<double> &log_area_Simp,
                          vector<double> &log_area_diff) {
  
  // use lambda method to define a version of loglike with p fixed and with rho
  // as the only free parameter
  auto loglike_interms_rho = [&cref_n = n, &cref_N = N, &cref_p = p, &cref_prior_p_shape1 = prior_p_shape1,
                              &cref_prior_p_shape2 = prior_p_shape2, &cref_prior_rho_shape1 = prior_rho_shape1,
                              &cref_prior_rho_shape2 = prior_rho_shape2](auto rho) {
                                return loglike_joint(cref_n, cref_N, cref_p, rho, cref_prior_p_shape1, cref_prior_p_shape2, cref_prior_rho_shape1, cref_prior_rho_shape2);
                              };
  
  // integrate over rho via adaptive quadrature
  adaptive_quadrature(loglike_interms_rho, n_intervals, 0.0, 1.0, x0, xm, x1,
                      log_y0, log_ym, log_y1, log_area_trap, log_area_Simp, log_area_diff);
  
  // get total area
  double log_area_sum = 0.0;
  auto it = max_element(log_area_Simp.begin(), log_area_Simp.end());
  double mx = (*it);
  if (!isfinite(mx)) {
    log_area_sum = -INFINITY;
  } else {
    double tmp = 0.0;
    for (int i = 0; i < n_intervals; ++i) {
      tmp += exp(log_area_Simp[i] - mx);
    }
    log_area_sum = mx + log(tmp);
  }
  
  return log_area_sum;
}

//------------------------------------------------
// prior*likelihood integrated over p

double loglike_marginal_rho(const vector<int> &n,
                            const vector<int> &N,
                            double rho, int n_intervals,
                            double prior_p_shape1, double prior_p_shape2,
                            double prior_rho_shape1, double prior_rho_shape2,
                            vector<double> &x0, vector<double> &xm,
                            vector<double> &x1, vector<double> &log_y0,
                            vector<double> &log_ym, vector<double> &log_y1,
                            vector<double> &log_area_trap,
                            vector<double> &log_area_Simp,
                            vector<double> &log_area_diff) {
  
  // use lambda method to define a version of loglike with rho fixed and with p
  // as the only free parameter
  auto loglike_interms_p = [&cref_n = n, &cref_N = N, &cref_rho = rho, &cref_prior_p_shape1 = prior_p_shape1,
                            &cref_prior_p_shape2 = prior_p_shape2, &cref_prior_rho_shape1 = prior_rho_shape1,
                            &cref_prior_rho_shape2 = prior_rho_shape2](auto p) {
                              return loglike_joint(cref_n, cref_N, p, cref_rho, cref_prior_p_shape1, cref_prior_p_shape2, cref_prior_rho_shape1, cref_prior_rho_shape2);
                            };
  
  // integrate over rho via adaptive quadrature
  adaptive_quadrature(loglike_interms_p, n_intervals, 0.0, 1.0, x0, xm, x1,
                      log_y0, log_ym, log_y1, log_area_trap, log_area_Simp, log_area_diff);
  
  // get total area
  double log_area_sum = 0.0;
  auto it = max_element(log_area_Simp.begin(), log_area_Simp.end());
  double mx = (*it);
  if (!isfinite(mx)) {
    log_area_sum = -INFINITY;
  } else {
    double tmp = 0.0;
    for (int i = 0; i < n_intervals; ++i) {
      tmp += exp(log_area_Simp[i] - mx);
    }
    log_area_sum = mx + log(tmp);
  }
  
  return log_area_sum;
}

//------------------------------------------------
// TODO

Rcpp::List get_credible_prevalence_cpp(Rcpp::List args_params) {
  
  // extract inputs
  vector<int> n = rcpp_to_vector_int(args_params["n"]);
  vector<int> N = rcpp_to_vector_int(args_params["N"]);
  double alpha = rcpp_to_double(args_params["alpha"]);
  double prior_p_shape1 = rcpp_to_double(args_params["prior_prev_shape1"]);
  double prior_p_shape2 = rcpp_to_double(args_params["prior_prev_shape2"]);
  double prior_rho_shape1 = rcpp_to_double(args_params["prior_ICC_shape1"]);
  double prior_rho_shape2 = rcpp_to_double(args_params["prior_ICC_shape2"]);
  int n_intervals = rcpp_to_int(args_params["n_intervals"]);
  
  // define objects to be used repeatedly within inner adaptive quadrature
  vector<double> Qx0(n_intervals);
  vector<double> Qxm(n_intervals);
  vector<double> Qx1(n_intervals);
  vector<double> Qlog_y0(n_intervals);
  vector<double> Qlog_ym(n_intervals);
  vector<double> Qlog_y1(n_intervals);
  vector<double> Qlog_area_trap(n_intervals);
  vector<double> Qlog_area_Simp(n_intervals);
  vector<double> Qlog_area_diff(n_intervals);
  
  // use lambda method to define a version of loglike with p as the only free
  // parameter, marginalised over rho
  auto loglike_interms_p = [&cref_n = n, &cref_N = N, &cref_n_intervals = n_intervals, &cref_prior_p_shape1 = prior_p_shape1,
                            &cref_prior_p_shape2 = prior_p_shape2, &cref_prior_rho_shape1 = prior_rho_shape1,
                            &cref_prior_rho_shape2 = prior_rho_shape2, &cref_x0 = Qx0, &cref_xm = Qxm, &cref_x1 = Qx1,
                            &cref_log_y0 = Qlog_y0, &cref_log_ym = Qlog_ym, &cref_log_y1 = Qlog_y1,
                            &cref_log_area_trap = Qlog_area_trap, &cref_log_area_Simp = Qlog_area_Simp,
                            &cref_log_area_diff = Qlog_area_diff](auto p) {
                              return loglike_marginal_p(cref_n, cref_N, p, cref_n_intervals, cref_prior_p_shape1, cref_prior_p_shape2,
                                                        cref_prior_rho_shape1, cref_prior_rho_shape2,
                                                        cref_x0, cref_xm, cref_x1, cref_log_y0, cref_log_ym, cref_log_y1,
                                                        cref_log_area_trap, cref_log_area_Simp, cref_log_area_diff);
                            };
  
  // define objects for storing outer quadrature results
  vector<double> x0(n_intervals);
  vector<double> xm(n_intervals);
  vector<double> x1(n_intervals);
  vector<double> log_y0(n_intervals);
  vector<double> log_ym(n_intervals);
  vector<double> log_y1(n_intervals);
  vector<double> log_area_trap(n_intervals);
  vector<double> log_area_Simp(n_intervals);
  vector<double> log_area_diff(n_intervals);
  
  // integrate over p via adaptive quadrature
  adaptive_quadrature(loglike_interms_p, n_intervals, 0.0, 1.0, x0, xm, x1,
                      log_y0, log_ym, log_y1, log_area_trap, log_area_Simp, log_area_diff);
  
  // return outputs as vector
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("x0") = x0,
                                      Rcpp::Named("xm") = xm,
                                      Rcpp::Named("x1") = x1,
                                      Rcpp::Named("log_y0") = log_y0,
                                      Rcpp::Named("log_ym") = log_ym,
                                      Rcpp::Named("log_y1") = log_y1,
                                      Rcpp::Named("log_area_trap") = log_area_trap,
                                      Rcpp::Named("log_area_Simp") = log_area_Simp,
                                      Rcpp::Named("log_area_diff") = log_area_diff);
  return ret;
}

//------------------------------------------------
// TODO

Rcpp::List get_credible_ICC_cpp(Rcpp::List args_params) {
  
  // extract inputs
  vector<int> n = rcpp_to_vector_int(args_params["n"]);
  vector<int> N = rcpp_to_vector_int(args_params["N"]);
  double alpha = rcpp_to_double(args_params["alpha"]);
  double prior_p_shape1 = rcpp_to_double(args_params["prior_prev_shape1"]);
  double prior_p_shape2 = rcpp_to_double(args_params["prior_prev_shape2"]);
  double prior_rho_shape1 = rcpp_to_double(args_params["prior_ICC_shape1"]);
  double prior_rho_shape2 = rcpp_to_double(args_params["prior_ICC_shape2"]);
  int n_intervals = rcpp_to_int(args_params["n_intervals"]);
  
  // define objects to be used repeatedly within inner adaptive quadrature
  vector<double> Qx0(n_intervals);
  vector<double> Qxm(n_intervals);
  vector<double> Qx1(n_intervals);
  vector<double> Qlog_y0(n_intervals);
  vector<double> Qlog_ym(n_intervals);
  vector<double> Qlog_y1(n_intervals);
  vector<double> Qlog_area_trap(n_intervals);
  vector<double> Qlog_area_Simp(n_intervals);
  vector<double> Qlog_area_diff(n_intervals);
  
  // use lambda method to define a version of loglike with rho as the only free
  // parameter, marginalised over p
  auto loglike_interms_rho = [&cref_n = n, &cref_N = N, &cref_n_intervals = n_intervals, &cref_prior_p_shape1 = prior_p_shape1,
                              &cref_prior_p_shape2 = prior_p_shape2, &cref_prior_rho_shape1 = prior_rho_shape1,
                              &cref_prior_rho_shape2 = prior_rho_shape2, &cref_x0 = Qx0, &cref_xm = Qxm, &cref_x1 = Qx1,
                              &cref_log_y0 = Qlog_y0, &cref_log_ym = Qlog_ym, &cref_log_y1 = Qlog_y1,
                              &cref_log_area_trap = Qlog_area_trap, &cref_log_area_Simp = Qlog_area_Simp,
                              &cref_log_area_diff = Qlog_area_diff](auto rho) {
                                return loglike_marginal_rho(cref_n, cref_N, rho, cref_n_intervals, cref_prior_p_shape1, cref_prior_p_shape2,
                                                            cref_prior_rho_shape1, cref_prior_rho_shape2,
                                                            cref_x0, cref_xm, cref_x1, cref_log_y0, cref_log_ym, cref_log_y1,
                                                            cref_log_area_trap, cref_log_area_Simp, cref_log_area_diff);
                              };
  
  // define objects for storing outer quadrature results
  vector<double> x0(n_intervals);
  vector<double> xm(n_intervals);
  vector<double> x1(n_intervals);
  vector<double> log_y0(n_intervals);
  vector<double> log_ym(n_intervals);
  vector<double> log_y1(n_intervals);
  vector<double> log_area_trap(n_intervals);
  vector<double> log_area_Simp(n_intervals);
  vector<double> log_area_diff(n_intervals);
  
  // integrate over p via adaptive quadrature
  adaptive_quadrature(loglike_interms_rho, n_intervals, 0.0, 1.0, x0, xm, x1,
                      log_y0, log_ym, log_y1, log_area_trap, log_area_Simp, log_area_diff);
  
  // return outputs as vector
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("x0") = x0,
                                      Rcpp::Named("xm") = xm,
                                      Rcpp::Named("x1") = x1,
                                      Rcpp::Named("log_y0") = log_y0,
                                      Rcpp::Named("log_ym") = log_ym,
                                      Rcpp::Named("log_y1") = log_y1,
                                      Rcpp::Named("log_area_trap") = log_area_trap,
                                      Rcpp::Named("log_area_Simp") = log_area_Simp,
                                      Rcpp::Named("log_area_diff") = log_area_diff);
  return ret;
}
