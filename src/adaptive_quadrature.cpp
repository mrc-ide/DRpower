
#include <math.h>

#include "adaptive_quadrature.h"

using namespace std;

//------------------------------------------------
// TODO

double get_trap_area(double x0, double x1, double log_y0, double log_y1) {
  double ret = log_sum(log_y0, log_y1) + log(x1 - x0) - log(2);
  return ret;
}

//------------------------------------------------
// TODO

double get_trap_area_double(double x0, double xm, double x1, double log_y0,
                            double log_ym, double log_y1) {
  double a1 = get_trap_area(x0, xm, log_y0, log_ym);
  double a2 = get_trap_area(xm, x1, log_ym, log_y1);
  return log_sum(a1, a2);
}

//------------------------------------------------
// TODO

double get_Simp_area(double x0, double x1, double log_y0, double log_ym,
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
double get_log_area_diff(double log_area1, double log_area2) {
  double ret = 0.0;
  if (!isfinite(log_area1) && !isfinite(log_area2)) {
    ret = -INFINITY;
  } else if (log_area1 > log_area2) {
    ret = log_area1 + log(1 - exp(log_area2 - log_area1));
  } else {
    ret = log_area2 + log(1 - exp(log_area1 - log_area2));
  }
  return ret;
}

//------------------------------------------------
double get_grad_diff(double x0, double x1, double log_y0, double log_y1,
                     double log_ym, double log_yd, double delta) {
  
  // get two angles between [-0.5, 0.5] and subtract to get relative difference,
  // then normalise to [0,1] scale where 0 represents no differece and 1
  // represents maximal difference
  double theta1 = atan((exp(log_y1) - exp(log_y0)) / (x1 - x0)) / M_PI;
  double theta2 = atan((exp(log_yd) - exp(log_ym)) / (delta*(x1 - x0)*0.5)) / M_PI;
  return abs(theta1 - theta2);
}

//------------------------------------------------
double get_total_diff(double log_area_diff, double grad_diff) {
  return exp(log_area_diff) * grad_diff;
}

//------------------------------------------------
// TODO

void adaptive_quadrature(function<double(double)> loglike, int n_intervals,
                         double left, double right, double delta,
                         vector<double> &x0, vector<double> &xm,
                         vector<double> &x1, vector<double> &log_y0,
                         vector<double> &log_ym, vector<double> &log_y1,
                         vector<double> &log_area_trap,
                         vector<double> &log_area_Simp,
                         vector<double> &total_diff) {
  
  // create first interval over entire domain
  x0[0] = left;
  x1[0] = right;
  xm[0] = 0.5*(left + right);
  
  log_y0[0] = loglike(x0[0]);
  log_y1[0] = loglike(x1[0]);
  log_ym[0] = loglike(xm[0]);
  
  total_diff[0] = 1;
  
  // loop through remaining intervals
  for (int i = 1; i < n_intervals; ++i) {
    
    // find which interval contains largest diff
    auto it = max_element(total_diff.begin(), total_diff.begin() + i);
    int w = distance(total_diff.begin(), it);
    
    // create new entry from xm to x1
    x0[i] = xm[w];
    x1[i] = x1[w];
    xm[i] = 0.5*(x0[i] + x1[i]);
    double xd = xm[i] + (x1[i] - xm[i])*delta;
    
    log_y0[i] = log_ym[w];
    log_y1[i] = log_y1[w];
    log_ym[i] = loglike(xm[i]);
    double log_yd = loglike(xd);
    
    log_area_trap[i] = get_trap_area_double(x0[i], xm[i], x1[i], log_y0[i], log_ym[i], log_y1[i]);
    log_area_Simp[i] = get_Simp_area(x0[i], x1[i], log_y0[i], log_ym[i], log_y1[i]);
    double log_area_diff = get_log_area_diff(log_area_trap[i], log_area_Simp[i]);
    double grad_diff = get_grad_diff(x0[i], x1[i], log_y0[i], log_y1[i], log_ym[i], log_yd, delta);
    total_diff[i] = get_total_diff(log_area_diff, grad_diff);
    
    // modify existing entry to go from x0 to xm
    x1[w] = xm[w];
    xm[w] = 0.5*(x0[w] + x1[w]);
    xd = xm[w] + (x1[w] - xm[w])*delta;
    
    log_y1[w] = log_ym[w];
    log_ym[w] = loglike(xm[w]);
    log_yd = loglike(xd);
    
    log_area_trap[w] = get_trap_area_double(x0[w], xm[w], x1[w], log_y0[w], log_ym[w], log_y1[w]);
    log_area_Simp[w] = get_Simp_area(x0[w], x1[w], log_y0[w], log_ym[w], log_y1[w]);
    log_area_diff = get_log_area_diff(log_area_trap[w], log_area_Simp[w]);
    grad_diff = get_grad_diff(x0[w], x1[w], log_y0[w], log_y1[w], log_ym[w], log_yd, delta);
    total_diff[w] = get_total_diff(log_area_diff, grad_diff);
  }
  
}
