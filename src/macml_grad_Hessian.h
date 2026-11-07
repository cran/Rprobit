#ifndef RPROBIT_MACML_GRAD_HESSIAN_H
#define RPROBIT_MACML_GRAD_HESSIAN_H

#include <RcppEigen.h>

#include "SJ.h"
#include "lin_alg.h"

Rcpp::List syst_from_mod (Eigen::VectorXd theta, Rcpp::List mod);
Eigen::MatrixXd J_diff (int mm, int kk);
Eigen::MatrixXd J_theta (
  Eigen::MatrixXd X,
  Eigen::VectorXd b,
  Eigen::MatrixXd sigma
);
Eigen::MatrixXd H_theta (
  Eigen::MatrixXd X,
  Eigen::VectorXd b,
  Eigen::MatrixXd sigma
);
Eigen::MatrixXd J_chol1 (
  Eigen::VectorXd x,
  Eigen::MatrixXd L,
  Eigen::MatrixXd K,
  Eigen::MatrixXd D
);
Eigen::MatrixXd J_bvOvSig_bLOLSig (
  Eigen::VectorXd tho,
  Eigen::VectorXd thl,
  int M,
  Eigen::MatrixXd Lo,
  Eigen::MatrixXd Ko,
  Eigen::MatrixXd Do,
  Eigen::MatrixXd Ll,
  Eigen::MatrixXd Kl,
  Eigen::MatrixXd Dl
);
Eigen::MatrixXd H_bvOvSig_bLOLSig (
  Eigen::MatrixXd J,
  int M,
  int O,
  int L
);
Eigen::MatrixXd J_Omega1 (
  Eigen::MatrixXd X,
  Eigen::MatrixXd L,
  Eigen::MatrixXd D
);
Eigen::MatrixXd makeM (int m, int y);
Eigen::MatrixXd J_TpdM_red (int alt, Eigen::VectorXd y);
Eigen::MatrixXd J_M_bOmSig_red (
  int alt,
  Rcpp::IntegerVector ind_BI,
  Eigen::VectorXd y,
  int lthbb,
  Eigen::MatrixXd X,
  Eigen::MatrixXd L,
  Eigen::MatrixXd D
);
Eigen::MatrixXd cal_pairs_matrix (int Tp_n, int cml_pair_type);
Rcpp::List cal_grad_Hess (
  Eigen::VectorXd x_norm,
  Eigen::MatrixXd Lambda_cor,
  Eigen::MatrixXd Lambda,
  Rcpp::List mod,
  Eigen::MatrixXd X_red,
  int cur,
  int alt,
  Rcpp::IntegerVector ind,
  Eigen::VectorXd y_n,
  Eigen::MatrixXd D_JOm,
  Rcpp::List syst,
  Eigen::MatrixXd L_JcholO,
  Eigen::MatrixXd K_JcholO,
  Eigen::MatrixXd D_JcholO,
  Eigen::MatrixXd L_JcholL,
  Eigen::MatrixXd K_JcholL,
  Eigen::MatrixXd D_JcholL,
  Eigen::MatrixXd HJ2,
  int hess
);

#endif
