#ifndef RPROBIT_DISTRIB_H
#define RPROBIT_DISTRIB_H

#include <RcppEigen.h>

extern double tol;
extern double bound_max;

// Normal distribution helpers with bounded evaluation points.
double std_normal_cdf (double x);
double std_normal_pdf (double x);
double normal_cdf (double x, double v);
double normal_pdf (double x, double v);

// Bivariate standard normal distribution and derivatives.
double biv_normal_cdf (double x0, double x1, double r12);
double biv_normal_pdf (double x0, double x1, double r12);
Eigen::VectorXd grad_cdf (double x0, double x1, double r12);
Eigen::MatrixXd Hess_cdf (double x0, double x1, double r12);

#endif
