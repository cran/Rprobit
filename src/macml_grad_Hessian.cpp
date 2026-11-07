#include "macml_grad_Hessian.h"

#include <cmath>

// This file contains helper functions that are needed to
// calculate the gradient and the Hessian of the probit and MacML likelihoods.
//
// J_XY denotes functions for Jacobians needed for gradients and Hessians.
// H_XY denotes functions of second derivatives.

// Returns the vectors and covariance matrices that are implied by theta and the
// linear restrictions stored in mod.

Rcpp::List syst_from_mod (Eigen::VectorXd theta, Rcpp::List mod) {

  // initialize output
  Rcpp::List out;

  //............................................................//
  // UNPACK MODEL SPEC
  //............................................................//

  int length_theta_b = mod["length_theta_b"];
  int length_theta_O = mod["length_theta_O"];
  int length_theta_L = mod["length_theta_L"];
  Rcpp::IntegerVector b_count = mod["b_count"];
  int length_RE = b_count[1];
  int alt = mod["alt"];
  Eigen::MatrixXd H_b = mod["H_b"];
  Eigen::VectorXd f_b = mod["f_b"];
  Eigen::VectorXd f_L = mod["f_L"];
  Eigen::MatrixXd H_L = mod["H_L"];
  Eigen::VectorXd f_O = mod["f_O"];
  Eigen::MatrixXd H_O = mod["H_O"];

  //............................................................//
  // RECOVER PARAMETERS
  //............................................................//

  // Recover the vector b
  Eigen::VectorXd b = ((H_b * theta.head (length_theta_b)) + f_b);

  out["b"] = b;

  // Recover the Sigma matrix
  Eigen::MatrixXd ELt = elimmat (alt).transpose ();
  Eigen::MatrixXd Sigma_chol_vech;

  // Without free parameters only f_L is used
  if (length_theta_L == 0) {
    Sigma_chol_vech = f_L;
  } else {
    // Sigma_chol_vech = ((H_L*l) + f_L);
    Eigen::VectorXd theta_L = theta.segment (
      length_theta_O + length_theta_b,
      length_theta_L
    );
    Sigma_chol_vech = H_L * theta_L + f_L;
  }

  out["Sigma_chol_vech"] = Sigma_chol_vech;

  Eigen::MatrixXd Sigma_chol = ELt * Sigma_chol_vech;
  Sigma_chol.resize (alt, alt);
  Eigen::MatrixXd Sigma_chol_t = Sigma_chol.transpose ();
  Eigen::MatrixXd Sigma = Sigma_chol * Sigma_chol_t;

  out["Sigma"] = Sigma;
  out["Sigma_chol"] = Sigma_chol;
  // Recover the matrix Omega
  Eigen::MatrixXd EOt = elimmat (length_RE).transpose ();
  Eigen::MatrixXd Omega_chol_vech, Omega_chol, Omega_chol_t, Omega;

  // check if we have mixed parameters
  if (length_RE != 0) {
    // Without free parameters only f_O is used
    if (length_theta_O == 0) {
      Omega_chol_vech = f_O;
    } else {
      // Omega_chol_vech = ((H_O*o) + f_O);
      Eigen::VectorXd theta_O = theta.segment (
        length_theta_b,
        length_theta_O
      );
      Omega_chol_vech = H_O * theta_O + f_O;
    }
    Omega_chol = EOt * Omega_chol_vech;
    Omega_chol.resize (length_RE, length_RE);
    Omega_chol_t = Omega_chol.transpose ();
    Omega = Omega_chol * Omega_chol_t;
  } else {
    Omega_chol_vech = f_O;
    Omega_chol = f_O;
    Omega_chol_t = f_O;
    Omega = f_O;
  }

  out["Omega"] = Omega;
  out["Omega_chol_vech"] = Omega_chol_vech;
  out["Omega_chol"] = Omega_chol;

  return out;
}


Eigen::MatrixXd J_diff (int mm, int kk) {
  // Computes the first derivative of vech(M_k Sig M_k') w.r.t. vech(Sig),
  // where Sig is an mm-by-mm matrix. The Hessian is zero. The function uses
  // k = kk - 1 to comply with C++ indexing.

  int m = (((mm * mm) - mm) / 2 + mm);
  int nn = mm - 1;
  int n = (((nn * nn) - nn) / 2 + nn);

  Eigen::MatrixXd J (n, m);
  J.setZero ();

  // compute entries of Sig

  Eigen::MatrixXd indSig = vechor_diag (mm);

  int k = kk - 1;

  // compute entries of M * Sig * t(M)

  Eigen::VectorXi row_selected = (indSig.col (0).array () != k).cast<int> ();
  Eigen::VectorXi col_selected = (indSig.col (1).array () != k).cast<int> ();
  Eigen::VectorXi selected = (row_selected.array () * col_selected.array ());

  Eigen::MatrixXd indMSigM (selected.sum (), 2);
  int rownew = 0;
  for (int i = 0; i < indSig.rows (); ++i) {
    if (selected[i]) {
      indMSigM.row (rownew) = indSig.row (i);
      rownew++;
    }
  }

  // j is the rows of J
  // i is the cols of J
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {

      // Every entry gets diffed by the diag entry.
      if ((indSig (i, 0) == k) && (indSig (i, 1) == k)) {
        J (j, i) = 1;
      }

      // Derivative for each entry with itself.
      if (indSig.row (i) == indMSigM.row (j)) {
        J (j, i) = 1;
      }

      // Off-diagonal entry in the diffed matrix I.
      if ((indSig (i, 1) == k) && (indSig (i, 0) == indMSigM (j, 0))) {
        J (j, i) = -1;
      }

      // Off-diagonal entry in the diffed matrix II.
      if ((indSig (i, 1) == k) && (indSig (i, 0) == indMSigM (j, 1))) {
        J (j, i) = -1;
      }

      // Off-diagonal entry in the diffed matrix III.
      if ((indSig (i, 0) == k) && (indSig (i, 1) == indMSigM (j, 1))) {
        J (j, i) = -1;
      }

      // Off-diagonal entry in the diffed matrix IV.
      if ((indSig (i, 0) == k) && (indSig (i, 1) == indMSigM (j, 0))) {
        J (j, i) = -1;
      }

      // Diagonal entry in the diffed matrix I.
      if ((indSig (i, 0) == k) && (indSig (i, 1) == indMSigM (j, 0))
      && (indMSigM (j, 0) == indMSigM (j, 1))) {
        J (j, i) = -2;
      }

      // Diagonal entry in the diffed matrix II.
      if ((indSig (i, 1) == k) && (indSig (i, 0) == indMSigM (j, 0))
      && (indMSigM (j, 0) == indMSigM (j, 1))) {
        J (j, i) = -2;
      }
    }

  }

  return J;
}

Eigen::MatrixXd J_theta (
  Eigen::MatrixXd X, Eigen::VectorXd b,Eigen::MatrixXd sigma
) {

  // x: covariate matrix
  // b: linear coefficients (k)
  // sigma: covariance matrix XOOX + LL (k x k)
  //
  // Returns the Jacobian w.r.t. b and vech(M), M = XOOX + LL (k x k). Hence,
  // this moves from theta = [b / diag(sqrt(s_11)), ..., s*_12] to b, vech(M).

  // Picking the length of b from X avoids problems with t(b).
  int n = X.cols ();

  Eigen::VectorXd sig = sigma.diagonal ();
  int k = sigma.rows ();

  Eigen::VectorXd kk = (-1 * (X * b).array ());
  Eigen::VectorXd kk1 = kk.array () / sig.array ().sqrt ();
  Eigen::MatrixXd J (k + ((k * k) - k) / 2, n + ((k * k) - k) / 2 + k);
  J.setZero ();

  // Derivative of the lin-utilities w.r.t. b.
  Eigen::MatrixXd XX = X;
  XX.array ().colwise () /= sig.array ().sqrt ();
  J.block (0, 0, X.rows (), X.cols ()) = -1 * XX;

  int m = n; // C++ indexing (-1)-R.

  Eigen::VectorXd id (k);
  id.setZero ();

  // Derivative of the lin-utilities w.r.t. diag(sqrt(Sig)).

  Eigen::VectorXd ssig = 1 / sig.array ().pow (1.5);
  Eigen::VectorXd sqsig = sig.array ().sqrt ();

  for (int i = 0; i < k; ++i) {
    id (i) = m;

    J (i, m) = -0.5 * kk[i] * ssig[i];
    int mm = m + (k - i);
    m = mm;
  }

  Eigen::MatrixXd ind = vechor (k);

  // Derivative of the off-diagonals w.r.t. diag(sqrt(Sig)).
  // Could be done in one if-clause but it works and hence will not be changed.
  for (int i = 0; i < k; ++i) {
    for (int j = k; j < (k + ((k * k) - k) / 2); ++j) {

      if ((ind ((j - k), 0) == i) && (ind ((j - k), 1) != i)) {

        J (j, (int) id[i]) =
          (sigma ((int) ind ((j - k), 0), (int) ind ((j - k), 1))
          / sqsig[ind ((j - k), 1)]) * -0.5 * ssig[i];
      }

      if ((ind ((j - k), 1) == i) && (ind ((j - k), 0) != i)) {
        J (j, (int) id[i]) =
          (sigma ((int) ind ((j - k), 0), (int) ind ((j - k), 1))
          / sqsig[(int) ind ((j - k), 0)]) * (-0.5 * ssig[i]);
      }

    }
  }

  Eigen::MatrixXd ind1 = vechor_diag (k);

  // Derivative of the off-diagonals w.r.t. the off-diagonals.

  for (int i = k; i < J.rows (); ++i) {
    for (int j = n; j < J.cols (); ++j) {

      if ((ind (i - k, 0) == ind1 (j - n, 0))
          && (ind (i - k, 1) == ind1 (j - n, 1))) {
        J (i, j) = 1
          / (sqsig[ind ((i - k), 0)] * sqsig[ind ((i - k), 1)]);
      }

    }
  }

  return J;
}

Eigen::MatrixXd H_theta (
  Eigen::MatrixXd X,
  Eigen::VectorXd b,
  Eigen::MatrixXd sigma
) {

  // x: covariate matrix
  // b: linear coefficients (k)
  // sigma: covariance matrix XOOX + LL (k x k)
  //
  // Returns the Hessian w.r.t. b and vech(sigma), where
  // sigma = XOOX + LL (k x k). Hence, this moves from
  // theta = [b / diag(sqrt(s_11)), ..., s*_12] to b, vech(sigma).

  // Picking the length of b from X avoids problems with t(b).
  int n = X.cols ();
  Eigen::VectorXd sig = sigma.diagonal ();
  int k = sigma.rows ();

  // Set up the Hessian.
  int h1 = k + ((k * k) - k) / 2;
  int h2 = n + ((k * k) - k) / 2 + k;
  Eigen::MatrixXd H (h2 * h1, h2);
  H.setZero ();

  // Identify indices and diagonals.
  Eigen::MatrixXd ind = vechor (k);
  Eigen::MatrixXd ind1 = vechor_diag (k);
  Eigen::VectorXd dd = fdiag (k);

  // Precompute relevant values.
  Eigen::VectorXd kk = (-1 * (X * b).array ());
  Eigen::VectorXd ssig15 = 1 / sig.array ().pow (1.5);
  Eigen::VectorXd ssig25 = 1 / sig.array ().pow (2.5);
  Eigen::MatrixXd ssigma15 = 1 / sigma.array ().pow (1.5);
  Eigen::MatrixXd ssigma25 = 1 / sigma.array ().pow (2.5);
  Eigen::MatrixXd ssigma = 1 / sigma.array ().sqrt ();

  int bb = 0;
  for (int i = 0; i < h1; ++i) {

    // Hessian for d b / d diag and d diag / d diag.
    // dbdb = 0
    // dbdOD = 0

    if (i < (k)) {

      H.block (bb, dd[i] + n, n, 1) =
        ((-1) * X.row (i) * (-0.5 * ssig15[i])).transpose ();
      H.block (dd[i] + n + bb, 0, 1, n) =
        ((-1) * X.row (i) * (-0.5 * ssig15[i]));
      H ((int) (dd[i] + n + bb), (int) (dd[i] + n)) = (kk[i])
        * ((0.75) * ssig25[i]);

    } else {

      // Off-diagonal entry shifted from the original index by the diagonals.
      int od = (i - k) + 1 + ind ((i - k), 1);
      int d1 = dd[ind ((i - k), 1)];
      int d2 = dd[ind ((i - k), 0)];

      // d od/d d1 (d1 and d2 are diagonals).
      H (bb + n + od, n + d1) = ssigma ((int) ind1 (d2, 0),
        (int) ind1 (d2, 0)) * (-0.5)
        * ssigma15 ((int) ind1 (d1, 0), (int) ind1 (d1, 0));
      H (bb + n + d1, n + od) = H (bb + n + od, n + d1);

      // d od/d d2.
      H (bb + n + od, n + d2) = ssigma ((int) ind1 (d1, 0),
        (int) ind1 (d1, 0)) * (-0.5)
        * ssigma15 ((int) ind1 (d2, 0), (int) ind1 (d2, 0));
      H (bb + n + d2, n + od) = H (bb + n + od, n + d2);

      // d d2/d d1.

      H (bb + n + d1, n + d2) = sigma ((int) ind1 (od, 0),
        (int) ind1 (od, 1)) * (-0.5)
        * ssigma15 ((int) ind1 (d2, 0), (int) ind1 (d2, 0)) * (-0.5)
        * ssigma15 ((int) ind1 (d1, 0), (int) ind1 (d1, 0));
      H (bb + n + d2, n + d1) = H (bb + n + d1, n + d2);

      // d d1/d d1.
      H (bb + n + d1, n + d1) = sigma ((int) ind1 (od, 0),
        (int) ind1 (od, 1))
        * ssigma ((int) ind1 (d2, 0), (int) ind1 (d2, 0)) * (0.75)
        * ssigma25 ((int) ind1 (d1, 0), (int) ind1 (d1, 0));

      // d d2/d d2.
      H (bb + n + d2, n + d2) = sigma ((int) ind1 (od, 0),
        (int) ind1 (od, 1))
        * ssigma ((int) ind1 (d1, 0), (int) ind1 (d1, 0)) * (0.75)
        * ssigma25 ((int) ind1 (d2, 0), (int) ind1 (d2, 0));

    }
    bb += h2;
  }

  return H;
}

Eigen::MatrixXd J_chol1 (
  Eigen::VectorXd x,
  Eigen::MatrixXd L,
  Eigen::MatrixXd K,
  Eigen::MatrixXd D
) {

  // Takes a vector of Cholesky factors (lower triangular, column-major) and
  // returns the column-major dvech(LL') / dvech(L) Jacobian.
  //
  // Uses precomputed elimination, commutation, and duplication matrices.

  int p = x.rows ();
  int m = -0.5 + std::sqrt (0.25 + 2 * p);

  Eigen::MatrixXd C (m, m);
  C.setZero ();

  int ii = 0;
  int jj = 0;

  for (int i = 0; i < p; ++i) {
    C (ii, jj) = x[i];
    ii++;
    if (ii == m) {
      ii = jj + 1;
      jj = ii;
    }
  }

  Eigen::MatrixXd Lt = L.transpose ();

  Eigen::MatrixXd I (m, m);
  Eigen::MatrixXd II (m * m, m * m);
  I.setIdentity ();
  II.setIdentity ();

  Eigen::MatrixXd J =
    L * (II + K) * Eigen::kroneckerProduct (C, I) * Lt;

  return J;
}

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
) {

  // Takes two vectors of lower-triangular, column-major Cholesky factors for O
  // and L and the number of linear parameters M.
  //
  // partial (b, vech(Omega), vech(Sigma)) /
  // partial (b, vech(L_Omega), vech(L_Sigma))

  int O = tho.rows ();
  int L = thl.rows ();

  Eigen::MatrixXd J1 = J_chol1 (tho, Lo, Ko, Do);
  Eigen::MatrixXd J2 = J_chol1 (thl, Ll, Kl, Dl);

  Eigen::MatrixXd JJ (M + O + L, M + O + L), Ind (M, M);
  JJ.setZero ();
  Ind.setIdentity ();

  JJ.block (0, 0, M, M) = Ind;
  JJ.block (M, M, O, O) = J1;
  JJ.block (M + O, M + O, L, L) = J2;

  return JJ;
}

Eigen::MatrixXd H_bvOvSig_bLOLSig (
  Eigen::MatrixXd J,
  int M,
  int O,
  int L
) {

  // Takes the Jacobian matrix and computes the Hessian in stacked form.
  // Entries are 1 except for diagonal entries, which are 2.
  int h1 = J.rows ();
  int h2 = J.cols ();

  Eigen::MatrixXd H (h1 * h2, h2);
  H.setZero ();

  Eigen::MatrixXd JJ (h1, h2);
  JJ.setZero ();
  JJ = (J.array () != 0).select (1, JJ);

  int LL = -0.5 + std::sqrt (0.25 + 2 * L);
  int OO = -0.5 + std::sqrt (0.25 + 2 * O);

  int bb = 0; // Number of parameters.

  int kl = (M + O);
  int ko = M;

  int jjl = 0;
  int jjo = 0;

  for (int i = 0; i < h1; ++i) {

    // Skip all the linear parameters because their Hessian is 0.
    if ((i > (M - 1)) && (i < (M + O))) {
      // Move through the O parameters.

      int e1 = -1;
      // Identify the diagonal entries in the final matrix. Those are squares
      // whenever J is positive.
      if (ko == i) {
        for (int j = M; j < h2; ++j) {
          if (JJ (i, j) != 0) {
            H (bb + j, j) = 2;
          }
        }
        ko = ko + (OO - jjo);
        jjo++;
      } else {
        for (int j = M; j < h2; ++j) {

          if (JJ (i, j) != 0) {
            if (e1 == -1) {
              e1 = j;
            } else {
              H (bb + j, (e1)) = 1;
              H (bb + (e1), (j)) = 1;
              e1 = -1;
            }
          }

        }
      }
    }

    if (i > (M + O - 1)) {

      // Move through the L parameters.

      int e1 = -1;
      // Identify the diagonal entries in the final matrix. Those are squares
      // whenever J is positive.
      if (kl == i) {
        for (int j = (M + O); j < h2; ++j) {
          if (JJ (i, j) != 0) {
            H (bb + j, j) = 2;
          }
        }
        kl = kl + (LL - jjl);
        jjl++;
      } else {
        for (int j = (M + O); j < h2; ++j) {

          if (JJ (i, j) != 0) {
            if (e1 == -1) {
              e1 = j;
            } else {
              H (bb + j, (e1)) = 1;
              H (bb + (e1), (j)) = 1;
              e1 = -1;
            }
          }
        }
      }
    }

    bb += h2;

  }

  return H;
}

Eigen::MatrixXd J_Omega1 (
  Eigen::MatrixXd X,
  Eigen::MatrixXd L,
  Eigen::MatrixXd D
) {
  // Precomputed duplication and elimination matrices improve performance.

  // Computes the Jacobian for vech(Omega).
  // X: coefficient matrix or, better, submatrix relevant to RE.

  // Basically just Theorem 2, page 35 in Magnus/Neudecker, plus the
  // elimination matrix for vec() to vech().
  // Uses precomputed elimination and duplication matrices.
  Eigen::MatrixXd J = L * Eigen::kroneckerProduct (X, X) * D;

  return J;

}

Eigen::MatrixXd makeM (int m, int y) {
  // m: size of the unrestricted/undiffed matrix.
  // y: row to remove, column to set to -1. Usually this will be choice - 1 due
  // to C++ indexing.
  Eigen::MatrixXd M (m, m);
  M.setIdentity ();
  M.col (y).setConstant (-1);
  removeRow (M, y);
  return M;
}

Eigen::MatrixXd J_TpdM_red (int alt, Eigen::VectorXd y) {
  // Computes the first derivative of M_Sig w.r.t. Sigma,
  // M_Sig = diag(M_k1 Sigma M_k1', M_k2 Sigma M_k2', ...).
  //
  // Sig is an LL matrix that appears as a block diagonal matrix with different
  // differencing. alt is the dimension of LL, equal to the number of choice
  // alternatives at the time point.
  int Tp = y.rows ();
  int dim_w = Tp * (alt - 1);

  Eigen::MatrixXd J (dim_w * (dim_w + 1) / 2, (alt + 1) * alt / 2);
  J.setZero ();

  int cur = 0;
  // int cur_col = 0;

  for (int jj = 0; jj < Tp; jj++) {
    Eigen::MatrixXd J1 = J_diff (alt, y[jj]); // choice.
    int curJ = 0;
    for (int a = 0; a < alt - 1; a++) {
      // a is column index
      for (int b = 0; b < alt - 1 - a; b++) {
        // b is row index
        if (cur < J.rows ()) {
          J.row (cur) = J1.row (curJ);
          curJ++;
          cur++;
        }
      }
      cur += dim_w - (jj + 1) * (alt - 1);
    }
  }

  return J;
}

Eigen::MatrixXd J_M_bOmSig_red (
  int alt,
  Rcpp::IntegerVector ind_BI,
  Eigen::VectorXd y,
  int lthbb,
  Eigen::MatrixXd X,
  Eigen::MatrixXd L,
  Eigen::MatrixXd D
) {
  // partial (b, vech(M)) / partial (b, vech(Omega), vech(Sigma)).
  // Precomputed duplication and elimination matrices improve performance.

  // Computes the first derivative for the first step.
  // b is identity.
  // O is D(XOOX) / D(OO).
  // L is D(block(diff(LL, k1), diff(LL, k2))) / D(LL).
  // At the moment only situations where the person has the same number of
  // choice alternatives are supported.

  // Inputs:

  // mm1, mm2 is dimension of LL, equal to the number of choice alternatives at
  // the time point.
  // k1, k2 choices at the time points.
  // M number of linear parameters.
  // X covariate submatrix of the REs.
  // L elimination matrix for J_Omega1.
  // D duplication matrix for J_Omega1.

  int cur = X.rows ();

  Eigen::MatrixXd J = J_TpdM_red (alt, y);
  int Tp = y.rows ();
  int rows_full = Tp * (alt - 1);
  Eigen::VectorXi ind (rows_full);
  ind.setZero ();

  for (int jj = 0; jj < cur; jj++) {
    ind[ind_BI[jj]] = 1;
  }

  // reduce due to missing values
  int cur_red = 0;
  int cur_full = 0;

  for (int a = 0; a < rows_full; a++) {
    if (ind[a] == 1) {
      // column is not reduced
      for (int b = a; b < rows_full; b++) {
        if (ind[b] == 1) {
          // row is not reduced
          J.row (cur_red) = J.row (cur_full);
          cur_red++;
        }
        cur_full++;
      }
    }
    if (ind[a] == 0) {
      // column is reduced, increase counter accordingly.
      cur_full = cur_full + rows_full - (a);
    }
  }

  Eigen::MatrixXd J_red = J.block (0, 0, cur_red, J.cols ());

  // Derivative of vech(X * Omega * X') w.r.t. entries in the lower triangle of
  // Omega, including the diagonal.
  Eigen::MatrixXd JO = J_Omega1 (X, L, D);

  // initialize output matrix.
  Eigen::MatrixXd JJ (J_red.rows () + lthbb,
    J_red.cols () + JO.cols () + lthbb), BdB (lthbb, lthbb);
  JJ.setZero ();
  BdB.setIdentity ();

  JJ.block (0, 0, lthbb, lthbb) = BdB;

  JJ.block (lthbb, lthbb, JO.rows (), JO.cols ()) = JO;
  JJ.block (lthbb, lthbb + JO.cols (), J_red.rows (), J_red.cols ()) = J_red;

  return JJ;
}

// Calculate the pairs matrix for Tp choices in an ordered probit model.

Eigen::MatrixXd cal_pairs_matrix (int Tp_n, int cml_pair_type) {

  // Define a matrix that allows to specify the pairs. This can be extended to
  // higher order CMLs.
  Eigen::MatrixXd pairs_n;

  Eigen::MatrixXd pairs_temp (Tp_n * (Tp_n - 1) / 2, 3);
  pairs_temp.setOnes ();
  // if the pairs are not predefined, build them
  int i_pair = 0;
  for (int i_pair_0 = 0; i_pair_0 < (Tp_n - 1); i_pair_0++) {
    for (int i_pair_1 = (i_pair_0 + 1); i_pair_1 < Tp_n; i_pair_1++) {

      if (cml_pair_type == 0) {
        // all pairwise weights
        // first choice occasion of the pair
        pairs_temp (i_pair, 0) = i_pair_0;
        // second choice occasion of the pair
        pairs_temp (i_pair, 1) = i_pair_1;
        pairs_temp (i_pair, 2) = 1;
        i_pair++;
      }
      // define the according weights to exclude undesired pairs
      if (cml_pair_type == 1) {
        // for adjacent pairs looped:
        if (((i_pair_0 + 1) == i_pair_1)
            || ((i_pair_0 == 0) && (i_pair_1 == (Tp_n - 1)))) {
          // first choice occasion of the pair
          pairs_temp (i_pair, 0) = i_pair_0;
          // second choice occasion of the pair
          pairs_temp (i_pair, 1) = i_pair_1;
          pairs_temp (i_pair, 2) = 1;
          i_pair++;
        }
      }
      if (cml_pair_type == 2) {
        // for adjacent pairs chained
        if ((i_pair_0 + 1) == i_pair_1) {
          // first choice occasion of the pair
          pairs_temp (i_pair, 0) = i_pair_0;
          // second choice occasion of the pair
          pairs_temp (i_pair, 1) = i_pair_1;
          pairs_temp (i_pair, 2) = 1;
          i_pair++;
        }
      }
    }
  }
  pairs_n = pairs_temp.block (0, 0, i_pair, 3);

  return pairs_n;
}

// Calculate the gradient and Hessian for different MNP model versions.

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
) {

  // initialize output
  Rcpp::List out (3);

  int length_theta_b = mod["length_theta_b"];
  int length_theta_O = mod["length_theta_O"];
  int length_theta_L = mod["length_theta_L"];
  int length_theta = length_theta_b + length_theta_O + length_theta_L;
  Rcpp::IntegerVector b_count = mod["b_count"];
  int length_RE = b_count[1];
  Eigen::MatrixXd H_b = mod["H_b"];
  Eigen::VectorXd f_b = mod["f_b"];
  Eigen::VectorXd f_L = mod["f_L"];
  Eigen::MatrixXd H_L = mod["H_L"];
  Eigen::VectorXd f_O = mod["f_O"];
  Eigen::MatrixXd H_O = mod["H_O"];


  Eigen::VectorXd b = syst["b"];
  int lthbb = b.rows ();
  int freeth = H_b.rows () + H_L.rows () + H_O.rows (); // bsum, Osum, Lsum;

  Eigen::MatrixXd Hess_n (freeth, freeth);
  Hess_n.setZero ();

  Eigen::VectorXd grad_n (freeth);
  grad_n.setZero ();

  // Recover the Sigma matrix
  Eigen::MatrixXd Sigma_chol_vech = syst["Sigma_chol_vech"];
  Eigen::MatrixXd Sigma_chol = syst["Sigma_chol"];
  Eigen::MatrixXd Sigma = syst["Sigma"];
  int lth_Sigma_chol_vech = Sigma_chol_vech.rows ();

  // Recover the matrix Omega
  Eigen::MatrixXd Omega_chol_vech = syst["Omega_chol_vech"];
  Eigen::MatrixXd Omega_chol = syst["Omega_chol"];
  Eigen::MatrixXd Omega = syst["Omega"];
  int lth_Omega_chol_vech = Omega_chol_vech.rows ();

  //............................................................//
  // Likelihood and Gradient
  //............................................................//

  Eigen::VectorXd grad_output;
  grad_output.setZero ();
  int lth_output;

  grad_output = SJ_grad (x_norm, Lambda_cor);

  lth_output = grad_output.size ();
  out["ll"] = grad_output[lth_output - 1];

  //
  //............................................................//
  // Gradient
  //............................................................//
  //

  Eigen::VectorXd grad_MVNCDF_approx;
  grad_MVNCDF_approx = grad_output.head (lth_output - 1);

  // initialize gradient contribution of the current pair
  Eigen::VectorXd grad_n_pair (length_theta);
  grad_n_pair.setZero ();
  //
  //............................................................//
  // START: Manipulations to go from MVNCDF approximation gradient
  //        to gradient of the original problem
  //............................................................//
  //
  // Calculate the Jacobian matrix of theta with respect to b and vech(Lambda).
  // Also calculate the Jacobian matrix of w, vech(Lambda_cor) with respect to
  // b and vech(Lambda).
  Eigen::MatrixXd J_wLambda_cor_bvechLambda = J_theta (
    X_red.block (0, 0, cur, X_red.cols()), b, Lambda);


  // partial log prob / partial (b, vech(Lambda)).
  Eigen::VectorXd g_bvechLambda = grad_MVNCDF_approx.transpose ()
    * J_wLambda_cor_bvechLambda;

  // partial (b, vech(Lambda), vech(X Omega X')) /
  // partial (b, vech(Omega), vech(Sigma))
  Eigen::MatrixXd L = elimmat (cur);

  Eigen::MatrixXd X_RE = X_red.block (0, X_red.cols() - length_RE , cur, length_RE);

  Eigen::MatrixXd J_dM_dbOmSig = J_M_bOmSig_red (
    alt, ind, y_n, lthbb, X_RE.block (0, 0, cur, X_RE.cols ()), L,
    D_JOm);

  // partial (b, vech(Omega), vech(Sigma)) /
  // partial (b, vech(L_Omega), vech(L_Sigma))
  Eigen::MatrixXd J_chol = J_bvOvSig_bLOLSig (Omega_chol_vech,
    Sigma_chol_vech, lthbb, L_JcholO, K_JcholO, D_JcholO, L_JcholL,
    K_JcholL, D_JcholL);

  Eigen::VectorXd g_bOmSig = g_bvechLambda.transpose () * J_dM_dbOmSig;

  Eigen::VectorXd g_bLOLSig = g_bOmSig.transpose () * J_chol;

  Eigen::MatrixXd gradbb = g_bLOLSig.head (lthbb);
  // H_b ... partial b / partial th_b.
  Eigen::MatrixXd gradb = H_b.transpose () * gradbb;

  // H_L ... partial vech(L) / partial th_Sig.
  Eigen::VectorXd gradL1 = g_bLOLSig.segment ((lth_Omega_chol_vech + lthbb),
    lth_Sigma_chol_vech);
  Eigen::MatrixXd gradL = H_L.transpose () * gradL1;

  // H_O ... partial vech(L_Omega) / partial th_Om.
  Eigen::VectorXd gradO1 = g_bLOLSig.segment ((lthbb), lth_Omega_chol_vech);
  Eigen::MatrixXd gradO = H_O.transpose () * gradO1;

  // fill in the gradients, depending on which is present.
  if (length_theta_L == 0) {
    if (length_theta_O == 0) {
      grad_n_pair << gradb;
    } else {
      grad_n_pair << gradb, gradO;
    }
  } else {
    if (length_theta_O == 0) {
      grad_n_pair << gradb, gradL;
    } else {
      grad_n_pair << gradb, gradO, gradL;
    }
  }

  //............................................................//
  // END of gradient manipulations
  //............................................................//

  // add gradient of the current pair to the gradient of the decision maker
  grad_n = (-1) * grad_n_pair;

  //............................................................//
  // Hessian
  //............................................................//

  if (hess == 1) {
    //............................................................//
    // exact hessian
    //............................................................//

    // First step
    Eigen::MatrixXd hes;
    hes = SJ_hess (x_norm, Lambda_cor);

    Eigen::MatrixXd H_th = H_theta (X_red.block (0, 0, cur, X_red.cols()), b,
      Lambda);

    Eigen::MatrixXd Ide (H_th.cols (), H_th.cols ());
    Ide.setIdentity ();
    Eigen::MatrixXd H_bvechLambda = J_wLambda_cor_bvechLambda.transpose ()
      * hes * J_wLambda_cor_bvechLambda
      + Eigen::kroneckerProduct (grad_MVNCDF_approx.transpose (), Ide)
      * H_th;

    // Second step
    Eigen::MatrixXd H_bOmSig = J_dM_dbOmSig.transpose () * H_bvechLambda
      * J_dM_dbOmSig;

    // Third step
    // b identity
    // Sigma_chol_vech & oL Cholesky
    Eigen::MatrixXd H_bLOLSig = H_bvOvSig_bLOLSig (HJ2, lthbb,
      lth_Omega_chol_vech, lth_Sigma_chol_vech);

    Eigen::SparseMatrix<double> sH_bLOLSig = H_bLOLSig.sparseView ();
    Eigen::SparseMatrix<double> Ide1 (H_bLOLSig.cols (),
      H_bLOLSig.cols ());
    Ide1.setIdentity ();

    Eigen::SparseMatrix<double> temp = Eigen::kroneckerProduct (
      g_bOmSig.transpose (), Ide1) * sH_bLOLSig;

    Eigen::MatrixXd H_t = Eigen::MatrixXd (temp);

    Eigen::MatrixXd Hess_n_pair =
      -J_chol.transpose () * H_bOmSig * J_chol - H_t;

    Eigen::MatrixXd J_H = Eigen::MatrixXd (length_theta, freeth);
    J_H.setZero ();

    J_H.block (0, 0, length_theta_b, lthbb) = H_b.transpose ();
    J_H.block (
      length_theta_b,
      lthbb,
      length_theta_O,
      lth_Omega_chol_vech
    ) = H_O.transpose ();
    J_H.block (
      length_theta_b + length_theta_O,
      lthbb + lth_Omega_chol_vech,
      length_theta_L,
      lth_Sigma_chol_vech
    ) = H_L.transpose ();

    // It is already negative in the pairs loop.
    Hess_n = J_H * Hess_n_pair * J_H.transpose ();

  }

  out["gradient"] = -(grad_n.array ());
  out["hessian"] = (Hess_n.array ());

  return out;
}
