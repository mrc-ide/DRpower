
#include <Rcpp.h>

#include <vector>
#include "misc_v16.h"

//------------------------------------------------
double get_logarea_trap(double x0, double x1, double log_y0, double log_y1);

//------------------------------------------------
double get_logarea_doubletrap(double x0, double x1, double log_y0, double log_ym,
                              double log_y1);

//------------------------------------------------
double get_logarea_Simp(double x0, double x1, double log_y0, double log_ym,
                        double log_y1);

//------------------------------------------------
double compare_grad(double x0, double x1, double log_y0, double log_y1,
                    double log_ym, double log_ydelta, double delta);

//------------------------------------------------
void adaptive_quadrature(std::function<double(double)> loglike, int n_intervals,
                         double left, double right,
                         std::vector<double> &x0, std::vector<double> &xm,
                         std::vector<double> &x1, std::vector<double> &log_y0,
                         std::vector<double> &log_ym, std::vector<double> &log_y1,
                         std::vector<double> &log_area_trap,
                         std::vector<double> &log_area_Simp,
                         std::vector<double> &log_area_diff);

