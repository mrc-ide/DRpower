
#include "adaptive_quadrature.h"

using namespace std;

//------------------------------------------------
// TODO

double get_logarea_trap(double x0, double x1, double log_y0, double log_y1) {
  double ret = 0.0;
  if (!isfinite(log_y0) && !isfinite(log_y1)) {
    ret = -INFINITY;
  } else if (log_y0 > log_y1) {
    ret = log_y0 + log(1 + exp(log_y1 - log_y0)) + log(x1 - x0) - log(2);
  } else {
    ret = log_y1 + log(1 + exp(log_y0 - log_y1)) + log(x1 - x0) - log(2);
  }
  return ret;
}

//------------------------------------------------
// TODO

double get_logarea_doubletrap(double x0, double x1, double log_y0, double log_ym,
                              double log_y1) {
  double a1 = get_logarea_trap(x0, 0.5*(x0 + x1), log_y0, log_ym);
  double a2 = get_logarea_trap(0.5*(x0 + x1), x1, log_ym, log_y1);
  return log_sum(a1, a2);
}

//------------------------------------------------
// TODO

double get_logarea_Simp(double x0, double x1, double log_y0, double log_ym,
                        double log_y1) {
  double ret = 0.0;
  if (!isfinite(log_y0) && !isfinite(log_ym) && !isfinite(log_y1)) {
    ret = -INFINITY;
  } else {
    double mx = max(max(log_y0, log_ym), log_y1);
    ret = mx + log(exp(log_y0 - mx) + 4*exp(log_ym - mx) + exp(log_y1 - mx)) + log(x1 - x0) - log(6.0);
  }
  return ret;
}

//------------------------------------------------
double compare_grad(double x0, double x1, double log_y0, double log_y1,
                    double log_ym, double log_ydelta, double delta) {
  double ret = 0.0;
  if (!isfinite(log_y0) && !isfinite(log_ym) && !isfinite(log_y1)) {
    ret = 1;
  } else {
    double grad_mid = (exp(log_ydelta) - exp(log_ym)) / delta;
    double grad_trap = (exp(log_y1) - exp(log_y0)) / (x1 - x0);
    ret = abs(grad_mid / grad_trap - 1.0);
  }
  return ret;
}

//------------------------------------------------
// TODO

void adaptive_quadrature(function<double(double)> loglike, int n_intervals,
                         double left, double right,
                         vector<double> &x0, vector<double> &xm,
                         vector<double> &x1, vector<double> &log_y0,
                         vector<double> &log_ym, vector<double> &log_y1,
                         vector<double> &log_area_trap,
                         vector<double> &log_area_Simp,
                         vector<double> &log_area_diff) {
  
  double delta = 1e-4;
  
  // create first interval
  x0[0] = left;
  x1[0] = right;
  xm[0] = 0.5*(left + right);
  log_y0[0] = loglike(x0[0]);
  log_y1[0] = loglike(x1[0]);
  log_ym[0] = loglike(xm[0]);
  
  // loop through remaining intervals
  for (int i = 1; i < n_intervals; ++i) {
    
    // find which element contains largest discrepancy in area between the two
    // methods
    auto it = max_element(log_area_diff.begin(), log_area_diff.begin() + i);
    int w = distance(log_area_diff.begin(), it);
    
    // create new entry and modify existing
    x0[i] = xm[w];
    x1[i] = x1[w];
    xm[i] = 0.5*(x0[i] + x1[i]);
    log_y0[i] = log_ym[w];
    log_y1[i] = log_y1[w];
    log_ym[i] = loglike(xm[i]);
    log_area_trap[i] = get_logarea_doubletrap(x0[i], x1[i], log_y0[i], log_ym[i], log_y1[i]);
    log_area_Simp[i] = get_logarea_Simp(x0[i], x1[i], log_y0[i], log_ym[i], log_y1[i]);
    double log_ydelta = loglike(xm[i] + delta);
    double grad_diff = compare_grad(x0[i], x1[i], log_y0[i], log_y1[i], log_ym[i], log_ydelta, delta);
    log_area_diff[i] = grad_diff * abs(exp(log_area_trap[i]) - exp(log_area_Simp[i]));
    
    x1[w] = xm[w];
    xm[w] = 0.5*(x0[w] + x1[w]);
    log_y1[w] = log_ym[w];
    log_ym[w] = loglike(xm[w]);
    log_area_trap[w] = get_logarea_doubletrap(x0[w], x1[w], log_y0[w], log_ym[w], log_y1[w]);
    log_area_Simp[w] = get_logarea_Simp(x0[w], x1[w], log_y0[w], log_ym[w], log_y1[w]);
    log_ydelta = loglike(xm[w] + delta);
    grad_diff = compare_grad(x0[w], x1[w], log_y0[w], log_y1[w], log_ym[w], log_ydelta, delta);
    log_area_diff[w] = grad_diff * abs(exp(log_area_trap[w]) - exp(log_area_Simp[w]));
    
  }
  
}
