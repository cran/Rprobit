#ifndef RPROBIT_SJ_H
#define RPROBIT_SJ_H

#include <RcppEigen.h>

double SJ (Eigen::VectorXd x, Eigen::MatrixXd r);
Eigen::VectorXd SJ_grad (Eigen::VectorXd x, Eigen::MatrixXd r);
Eigen::MatrixXd SJ_hess (Eigen::VectorXd x, Eigen::MatrixXd r);

#endif
