
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
    // special case if p=0 or p=1
    if (p == 0) {
      for (int i = 0; i < n_clust; ++ i) {
        if (n[i] != 0) {
          return -INFINITY;
        }
      }
      return 0.0;
    } else if (p == 1) {
      for (int i = 0; i < n_clust; ++ i) {
        if (n[i] != N[i]) {
          return -INFINITY;
        }
      }
      return 0.0;
    }
    // simplifies to binomial distribution
    double ret = 0.0;
    for (int i = 0; i < n_clust; ++ i) {
      ret += lgamma(N[i] + 1) - lgamma(n[i] + 1) - lgamma(N[i] - n[i] + 1) +
        n[i]*log(p) + (N[i] - n[i])*log(1.0 - p);
    }
    return ret;
  } else if (rho == 1.0) {
    // likelihood still positive whenever n == 0 or n == N
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
    // likelihood 1 if all k == 0, otherwise likelihood 0
    for (int i = 0; i < n_clust; ++ i) {
      if (n[i] != 0) {
        return -INFINITY;
      }
    }
    return 0.0;
  } else if (p == 1.0) {
    // likelihood 1 if all k == m, otherwise likelihood 0
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
  
  // deal first with edge cases under the prior
  double ret = 0.0;
  if (p == 0) {
    if (prior_p_shape1 > 1.0) {
      return -INFINITY;
    } else {
      ret += lgamma(prior_p_shape1 + prior_p_shape2) - lgamma(prior_p_shape1) - lgamma(prior_p_shape2);
    }
  } else if (p == 1) {
    if (prior_p_shape2 > 1.0) {
      return -INFINITY;
    } else {
      ret += lgamma(prior_p_shape1 + prior_p_shape2) - lgamma(prior_p_shape1) - lgamma(prior_p_shape2);
    }
  } else {
    ret += lgamma(prior_p_shape1 + prior_p_shape2) - lgamma(prior_p_shape1) - lgamma(prior_p_shape2) +
      (prior_p_shape1 - 1)*log(p) + (prior_p_shape2 - 1)*log(1.0 - p);
  }
  if (rho == 0) {
    if (prior_rho_shape1 > 1.0) {
      return -INFINITY;
    } else {
      ret += lgamma(prior_rho_shape1 + prior_rho_shape2) - lgamma(prior_rho_shape1) - lgamma(prior_rho_shape2);
    }
  } else if (rho == 1) {
    if (prior_rho_shape2 > 1.0) {
      return -INFINITY;
    } else {
      ret += lgamma(prior_rho_shape1 + prior_rho_shape2) - lgamma(prior_rho_shape1) - lgamma(prior_rho_shape2);
    }
  } else {
    ret += lgamma(prior_rho_shape1 + prior_rho_shape2) - lgamma(prior_rho_shape1) - lgamma(prior_rho_shape2) +
      (prior_rho_shape1 - 1)*log(rho) + (prior_rho_shape2 - 1)*log(1.0 - rho);
  }
  
  // combine with log-likelihood
  ret += dbbinom_reparam(n, N, p, rho);
  
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
                          vector<double> &total_diff) {
  
  // use lambda method to define a version of loglike with p fixed and with rho
  // as the only free parameter
  auto loglike_interms_rho = [&cref_n = n, &cref_N = N, &cref_p = p, &cref_prior_p_shape1 = prior_p_shape1,
                              &cref_prior_p_shape2 = prior_p_shape2, &cref_prior_rho_shape1 = prior_rho_shape1,
                              &cref_prior_rho_shape2 = prior_rho_shape2](auto rho) {
                                return loglike_joint(cref_n, cref_N, cref_p, rho, cref_prior_p_shape1, cref_prior_p_shape2, cref_prior_rho_shape1, cref_prior_rho_shape2);
                              };
  
  // integrate over rho via adaptive quadrature
  adaptive_quadrature(loglike_interms_rho, n_intervals, 0.0, 1.0, 1e-3, x0, xm, x1,
                      log_y0, log_ym, log_y1, log_area_trap, log_area_Simp, total_diff);
  
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
                            vector<double> &total_diff) {
  
  // use lambda method to define a version of loglike with rho fixed and with p
  // as the only free parameter
  auto loglike_interms_p = [&cref_n = n, &cref_N = N, &cref_rho = rho, &cref_prior_p_shape1 = prior_p_shape1,
                            &cref_prior_p_shape2 = prior_p_shape2, &cref_prior_rho_shape1 = prior_rho_shape1,
                            &cref_prior_rho_shape2 = prior_rho_shape2](auto p) {
                              return loglike_joint(cref_n, cref_N, p, cref_rho, cref_prior_p_shape1, cref_prior_p_shape2, cref_prior_rho_shape1, cref_prior_rho_shape2);
                            };
  
  // integrate over rho via adaptive quadrature
  adaptive_quadrature(loglike_interms_p, n_intervals, 0.0, 1.0, 1e-3, x0, xm, x1,
                      log_y0, log_ym, log_y1, log_area_trap, log_area_Simp, total_diff);
  
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

Rcpp::List get_prevalence_cpp(Rcpp::List args_params) {
  
  // extract inputs
  vector<int> n = rcpp_to_vector_int(args_params["n"]);
  vector<int> N = rcpp_to_vector_int(args_params["N"]);
  double prior_p_shape1 = rcpp_to_double(args_params["prior_prev_shape1"]);
  double prior_p_shape2 = rcpp_to_double(args_params["prior_prev_shape2"]);
  double prior_rho_shape1 = rcpp_to_double(args_params["prior_ICC_shape1"]);
  double prior_rho_shape2 = rcpp_to_double(args_params["prior_ICC_shape2"]);
  int n_intervals = rcpp_to_int(args_params["n_intervals"]);
  
  // define objects to be used repeatedly within inner adaptive quadrature
  vector<double> x0_inner(n_intervals);
  vector<double> xm_inner(n_intervals);
  vector<double> x1_inner(n_intervals);
  vector<double> log_y0_inner(n_intervals);
  vector<double> log_ym_inner(n_intervals);
  vector<double> log_y1_inner(n_intervals);
  vector<double> log_area_trap_inner(n_intervals);
  vector<double> log_area_Simp_inner(n_intervals);
  vector<double> total_diff_inner(n_intervals);
  
  // use lambda method to define a version of loglike with p as the only free
  // parameter, marginalised over rho
  auto loglike_interms_p = [&cref_n = n,
                            &cref_N = N,
                            &cref_n_intervals = n_intervals,
                            &cref_prior_p_shape1 = prior_p_shape1,
                            &cref_prior_p_shape2 = prior_p_shape2,
                            &cref_prior_rho_shape1 = prior_rho_shape1,
                            &cref_prior_rho_shape2 = prior_rho_shape2,
                            &cref_x0 = x0_inner,
                            &cref_xm = xm_inner,
                            &cref_x1 = x1_inner,
                            &cref_log_y0 = log_y0_inner,
                            &cref_log_ym = log_ym_inner,
                            &cref_log_y1 = log_y1_inner,
                            &cref_log_area_trap = log_area_trap_inner,
                            &cref_log_area_Simp = log_area_Simp_inner,
                            &cref_total_diff = total_diff_inner](auto p) {
                              return loglike_marginal_p(cref_n, cref_N, p, cref_n_intervals, cref_prior_p_shape1, cref_prior_p_shape2,
                                                        cref_prior_rho_shape1, cref_prior_rho_shape2, cref_x0, cref_xm, cref_x1,
                                                        cref_log_y0, cref_log_ym, cref_log_y1,
                                                        cref_log_area_trap, cref_log_area_Simp, cref_total_diff);
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
  vector<double> total_diff(n_intervals);
  
  // integrate over p via adaptive quadrature
  adaptive_quadrature(loglike_interms_p, n_intervals, 0.0, 1.0, 1e-3, x0, xm, x1,
                      log_y0, log_ym, log_y1, log_area_trap, log_area_Simp, total_diff);
  
  // return outputs as vector
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("x0") = x0,
                                      Rcpp::Named("xm") = xm,
                                      Rcpp::Named("x1") = x1,
                                      Rcpp::Named("log_y0") = log_y0,
                                      Rcpp::Named("log_ym") = log_ym,
                                      Rcpp::Named("log_y1") = log_y1,
                                      Rcpp::Named("log_area_trap") = log_area_trap,
                                      Rcpp::Named("log_area_Simp") = log_area_Simp,
                                      Rcpp::Named("total_diff") = total_diff);
  return ret;
}

//------------------------------------------------
// TODO

Rcpp::List get_ICC_cpp(Rcpp::List args_params) {
  
  // extract inputs
  vector<int> n = rcpp_to_vector_int(args_params["n"]);
  vector<int> N = rcpp_to_vector_int(args_params["N"]);
  double prior_p_shape1 = rcpp_to_double(args_params["prior_prev_shape1"]);
  double prior_p_shape2 = rcpp_to_double(args_params["prior_prev_shape2"]);
  double prior_rho_shape1 = rcpp_to_double(args_params["prior_ICC_shape1"]);
  double prior_rho_shape2 = rcpp_to_double(args_params["prior_ICC_shape2"]);
  int n_intervals = rcpp_to_int(args_params["n_intervals"]);
  
  // define objects to be used repeatedly within inner adaptive quadrature
  vector<double> x0_inner(n_intervals);
  vector<double> xm_inner(n_intervals);
  vector<double> x1_inner(n_intervals);
  vector<double> log_y0_inner(n_intervals);
  vector<double> log_ym_inner(n_intervals);
  vector<double> log_y1_inner(n_intervals);
  vector<double> log_area_trap_inner(n_intervals);
  vector<double> log_area_Simp_inner(n_intervals);
  vector<double> total_diff_inner(n_intervals);
  
  // use lambda method to define a version of loglike with rho as the only free
  // parameter, marginalised over p
  auto loglike_interms_rho = [&cref_n = n,
                              &cref_N = N,
                              &cref_n_intervals = n_intervals,
                              &cref_prior_p_shape1 = prior_p_shape1,
                              &cref_prior_p_shape2 = prior_p_shape2,
                              &cref_prior_rho_shape1 = prior_rho_shape1,
                              &cref_prior_rho_shape2 = prior_rho_shape2,
                              &cref_x0 = x0_inner,
                              &cref_xm = xm_inner, 
                              &cref_x1 = x1_inner,
                              &cref_log_y0 = log_y0_inner,
                              &cref_log_ym = log_ym_inner,
                              &cref_log_y1 = log_y1_inner,
                              &cref_log_area_trap = log_area_trap_inner,
                              &cref_log_area_Simp = log_area_Simp_inner,
                              &cref_total_diff = total_diff_inner](auto rho) {
                                return loglike_marginal_rho(cref_n, cref_N, rho, cref_n_intervals, cref_prior_p_shape1, cref_prior_p_shape2,
                                                            cref_prior_rho_shape1, cref_prior_rho_shape2, cref_x0, cref_xm, cref_x1,
                                                            cref_log_y0, cref_log_ym, cref_log_y1,
                                                            cref_log_area_trap, cref_log_area_Simp, cref_total_diff);
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
  vector<double> total_diff(n_intervals);
  
  // integrate over p via adaptive quadrature
  adaptive_quadrature(loglike_interms_rho, n_intervals, 0.0, 1.0, 1e-3, x0, xm, x1,
                      log_y0, log_ym, log_y1, log_area_trap, log_area_Simp, total_diff);
  
  // return outputs as vector
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("x0") = x0,
                                      Rcpp::Named("xm") = xm,
                                      Rcpp::Named("x1") = x1,
                                      Rcpp::Named("log_y0") = log_y0,
                                      Rcpp::Named("log_ym") = log_ym,
                                      Rcpp::Named("log_y1") = log_y1,
                                      Rcpp::Named("log_area_trap") = log_area_trap,
                                      Rcpp::Named("log_area_Simp") = log_area_Simp,
                                      Rcpp::Named("total_diff") = total_diff);
  return ret;
}
