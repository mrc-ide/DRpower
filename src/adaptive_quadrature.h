
#include <Rcpp.h>

#include <vector>
#include "misc_v16.h"

//------------------------------------------------
double get_trap_area(double x0, double x1, double log_y0, double log_y1);

//------------------------------------------------
double get_trap_area_double(double x0, double xm, double x1, double log_y0,
                            double log_ym, double log_y1);

//------------------------------------------------
double get_Simp_area(double x0, double x1, double log_y0, double log_ym,
                     double log_y1);

//------------------------------------------------
double get_log_area_diff(double log_area1, double log_area2);

//------------------------------------------------
double get_grad_diff(double x0, double x1, double log_y0, double log_y1,
                     double log_ym, double log_yd, double delta);

//------------------------------------------------
double get_total_diff(double log_area_diff, double grad_diff);

//------------------------------------------------
void adaptive_quadrature(std::function<double(double)> loglike, int n_intervals,
                         double left, double right, double delta,
                         std::vector<double> &x0, std::vector<double> &xm,
                         std::vector<double> &x1, std::vector<double> &log_y0,
                         std::vector<double> &log_ym, std::vector<double> &log_y1,
                         std::vector<double> &log_area_trap,
                         std::vector<double> &log_area_Simp,
                         std::vector<double> &total_diff);

