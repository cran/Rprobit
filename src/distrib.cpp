#include <cmath>

#include "distrib.h"
#include "toms462.h"

double tol = 0.000000001;
double bound_max = 6.0;

//' Standard normal cumulative distribution function
//'
//' @description
//' Computes the standard normal CDF after bounding the evaluation point.
//'
//' @param x \[`numeric(1)`\]\cr
//' Evaluation point.
//'
//' @return
//' CDF value.
//'
//' @keywords internal

double std_normal_cdf (double x) {
  double m = 0;
  double s = 1;
  if (x > bound_max) {
    x = bound_max;
  }
  if (x < -bound_max) {
    x = -bound_max;
  }
  return 0.5 * (1 + std::erf ((x - m) / (s * std::sqrt (2.))));
}

//' Standard normal probability density function
//'
//' @description
//' Computes the standard normal PDF after bounding the standardized point.
//'
//' @param x \[`numeric(1)`\]\cr
//' Evaluation point.
//'
//' @return
//' PDF value.
//'
//' @keywords internal

double std_normal_pdf (double x) {
  double m = 0;
  double s = 1;
  static const double inv_sqrt_2pi = 0.3989422804014327;
  double a = (x - m) / s;
  if (a > bound_max) {
    a = bound_max;
  }
  if (a < -bound_max) {
    a = -bound_max;
  }
  return inv_sqrt_2pi / s * std::exp (-0.5 * a * a);
}

//' Normal cumulative distribution function
//'
//' @description
//' Computes the CDF of a zero-mean normal distribution.
//'
//' @param x \[`numeric(1)`\]\cr
//' Evaluation point.
//'
//' @param v \[`numeric(1)`\]\cr
//' Variance.
//'
//' @return
//' CDF value.
//'
//' @keywords internal

double normal_cdf (double x, double v) {
  double s = std::sqrt (v);
  double m = x / s;
  if (m > bound_max) {
    m = bound_max;
  }
  if (m < -bound_max) {
    m = -bound_max;
  }
  return 0.5 * (1 + std::erf (m / std::sqrt (2.)));
}

//' Normal probability density function
//'
//' @description
//' Computes the PDF of a zero-mean normal distribution.
//'
//' @param x \[`numeric(1)`\]\cr
//' Evaluation point.
//'
//' @param v \[`numeric(1)`\]\cr
//' Variance.
//'
//' @return
//' PDF value.
//'
//' @keywords internal

double normal_pdf (double x, double v) {
  double m = 0;
  double s = std::sqrt (v);
  static const double inv_sqrt_2pi = 0.3989422804014327;
  double a = (x - m) / s;
  if (a > bound_max) {
    a = bound_max;
  }
  if (a < -bound_max) {
    a = -bound_max;
  }
  return inv_sqrt_2pi / s * std::exp (-0.5 * a * a);
}


//' Bivariate standard normal cumulative distribution function
//'
//' @description
//' Computes the bivariate standard normal CDF.
//'
//' @param x0 \[`numeric(1)`\]\cr
//' First evaluation point.
//'
//' @param x1 \[`numeric(1)`\]\cr
//' Second evaluation point.
//'
//' @param r12 \[`numeric(1)`\]\cr
//' Correlation.
//'
//' @return
//' CDF value, bounded below by `tol`.
//'
//' @keywords internal

double biv_normal_cdf (double x0, double x1, double r12) {
  double cd = 0;
  if (r12 > 1 - tol) {
    r12 = 1 - tol;
  }
  if (r12 < -tol - 1) {
    r12 = tol - 1;
  }
  if (x0 > bound_max) {
    x0 = bound_max;
  }
  if (x1 > bound_max) {
    x1 = bound_max;
  }
  if (x0 < -bound_max) {
    x0 = -bound_max;
  }
  if (x1 < -bound_max) {
    x1 = -bound_max;
  }
  cd = bivnor (-x0, -x1, r12);
  if (cd < tol) {
    cd = tol;
  }
  return cd;
}

//' Bivariate standard normal probability density function
//'
//' @description
//' Computes the bivariate standard normal PDF.
//'
//' @param x0 \[`numeric(1)`\]\cr
//' First evaluation point.
//'
//' @param x1 \[`numeric(1)`\]\cr
//' Second evaluation point.
//'
//' @param r12 \[`numeric(1)`\]\cr
//' Correlation.
//'
//' @return
//' PDF value.
//'
//' @keywords internal

double biv_normal_pdf (double x0, double x1, double r12) {
  static const double inv_2pi = 0.1591549; /* 1/(2pi) */
  double pd = 0;
  if (r12 > 1 - tol) {
    r12 = 1 - tol;
  }
  if (r12 < -tol - 1) {
    r12 = tol - 1;
  }
  if (x0 > bound_max) {
    x0 = bound_max;
  }
  if (x1 > bound_max) {
    x1 = bound_max;
  }
  if (x0 < -bound_max) {
    x0 = -bound_max;
  }
  if (x1 < -bound_max) {
    x1 = -bound_max;
  }
  pd = inv_2pi / std::sqrt (1 - r12 * r12)
    * std::exp (
        -0.5 * (x0 * x0 + x1 * x1 - 2 * x0 * x1 * r12) / (1 - r12 * r12));
        return pd;
}

//' Gradient of the bivariate standard normal CDF
//'
//' @description
//' Computes the gradient of the bivariate standard normal CDF. The
//' derivatives are ordered with respect to `x0`, `x1`, and `r12`.
//'
//' @param x0 \[`numeric(1)`\]\cr
//' First evaluation point.
//'
//' @param x1 \[`numeric(1)`\]\cr
//' Second evaluation point.
//'
//' @param r12 \[`numeric(1)`\]\cr
//' Correlation.
//'
//' @return
//' A numeric vector of length 3.
//'
//' @keywords internal

Eigen::VectorXd grad_cdf (double x0, double x1, double r12) {
  Eigen::VectorXd grad (3);
  grad.setZero ();
  // initialize
  if (r12 > 1 - tol) {
    r12 = 1 - tol;
  }
  if (r12 < -tol - 1) {
    r12 = tol - 1;
  }
  if (x0 > bound_max) {
    x0 = bound_max;
  }
  if (x1 > bound_max) {
    x1 = bound_max;
  }
  if (x0 < -bound_max) {
    x0 = -bound_max;
  }
  if (x1 < -bound_max) {
    x1 = -bound_max;
  }
  // derivatives
  grad(0) = std_normal_pdf (x0) * normal_cdf (x1 - r12 * x0, 1 - r12 * r12);
  grad(1) = std_normal_pdf (x1) * normal_cdf (x0 - r12 * x1, 1 - r12 * r12);
  grad(2) = biv_normal_pdf (x0, x1, r12);
  return grad;
}

//' Hessian of the bivariate standard normal CDF
//'
//' @description
//' Computes the Hessian of the bivariate standard normal CDF. Rows and
//' columns follow the derivative order `x0`, `x1`, `r12`.
//'
//' @param x0 \[`numeric(1)`\]\cr
//' First evaluation point.
//'
//' @param x1 \[`numeric(1)`\]\cr
//' Second evaluation point.
//'
//' @param r12 \[`numeric(1)`\]\cr
//' Correlation.
//'
//' @return
//' A 3x3 numeric Hessian matrix.
//'
//' @keywords internal

Eigen::MatrixXd Hess_cdf (double x0, double x1, double r12) {
  Eigen::MatrixXd Hess (3, 3);
  Hess.setZero ();
  // initialize
  if (r12 > 1 - tol) {
    r12 = 1 - tol;
  }
  if (r12 < -tol - 1) {
    r12 = tol - 1;
  }
  if (x0 > bound_max) {
    x0 = bound_max;
  }
  if (x1 > bound_max) {
    x1 = bound_max;
  }
  if (x0 < -bound_max) {
    x0 = -bound_max;
  }
  if (x1 < -bound_max) {
    x1 = -bound_max;
  }
  // gradient needed?
  Eigen::VectorXd gr (3);
  gr.setZero ();
  gr = grad_cdf (x0, x1, r12);
  // Hessian
  Hess (0, 0) = gr (0) * (-x0) + biv_normal_pdf (x0, x1, r12) * (-r12);
  Hess (0, 1) = std_normal_pdf (x0)
    * normal_pdf (x1 - r12 * x0, 1 - r12 * r12);
  Hess (0, 2) = gr (2) * (-x0 + x1 * r12) / (1 - r12 * r12);
  Hess (1, 0) = Hess (0, 1);
  Hess (1, 1) = gr (1) * (-x1) + biv_normal_pdf (x0, x1, r12) * (-r12);
  Hess (1, 2) = gr (2) * (-x1 + x0 * r12) / (1 - r12 * r12);
  Hess (2, 0) = Hess (0, 2);
  Hess (2, 1) = Hess (1, 2);
  // last with respect to two times the corr.
  double onemr2 = 1 - r12 * r12;
  double he = x0 * x1 * onemr2 - (x0 * x0 + x1 * x1 - 2 * r12 * x0 * x1) * r12;
  Hess (2, 2) = biv_normal_pdf (x0, x1, r12)
    * (he / (onemr2 * onemr2) + r12 / onemr2);
  return Hess;
}
