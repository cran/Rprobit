#ifndef RPROBIT_LIN_ALG_H
#define RPROBIT_LIN_ALG_H

#include <RcppEigen.h>

Eigen::VectorXd trilCols (int M);
Eigen::VectorXd trilRows (int M);

Eigen::MatrixXd duplmat (int n);
Eigen::MatrixXd vechor (int k);
Eigen::VectorXd vech (Eigen::MatrixXd Q);
Eigen::VectorXd fdiag (int k);
Eigen::MatrixXd vechor_diag (int k);
Eigen::MatrixXd commmat (int r, int c);
Eigen::MatrixXd elimmat (int n);

void removeRow (Eigen::MatrixXd& matrix, unsigned int rowToRemove);
void removeColumn (Eigen::MatrixXd& matrix, unsigned int colToRemove);
Rcpp::IntegerVector NA_ind (Eigen::MatrixXd Xc);

#endif
