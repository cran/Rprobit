#include "lin_alg.h"

//' Lower triangular column indices
//'
//' @description
//' Returns column indices for the lower triangular part of a square
//' matrix, excluding the diagonal.
//'
//' @param M \[`integer(1)`\]\cr
//' Matrix dimension.
//'
//' @return
//' A numeric vector of column indices.
//'
//' @keywords internal

Eigen::VectorXd trilCols (int M) {
  Eigen::VectorXd col ((M * (M - 1)) / 2);
  col.setZero ();
  int k;
  k = (M - 1);
  for (int j = 1; j < (M - 1); j++) {
    for (int i = 0; i < (M - j - 1); i++) {
      col (k) = j;
      k++;
    }
  }
  return col;
}

//' Lower triangular row indices
//'
//' @description
//' Returns row indices for the lower triangular part of a square
//' matrix, excluding the diagonal.
//'
//' @param M \[`integer(1)`\]\cr
//' Matrix dimension.
//'
//' @return
//' A numeric vector of row indices.
//'
//' @keywords internal

Eigen::VectorXd trilRows (int M) {
  Eigen::VectorXd row ((M * (M - 1)) / 2);
  int k;
  k = 0;
  for (int i = 1; i < M; i++) {
    for (int j = i; j < M; j++) {
      row (k) = j;
      k++;
    }
  }
  return row;
}

//' Duplication matrix
//'
//' @description
//' Computes the duplication matrix converting `vech(A)` to `vec(A)` for
//' a symmetric matrix `A`.
//'
//' @param n \[`integer(1)`\]\cr
//' Matrix dimension.
//'
//' @return
//' A numeric duplication matrix with `n * n` rows and `n * (n + 1) / 2`
//' columns.
//'
//' @keywords internal

Eigen::MatrixXd duplmat (int n) {
  int p = n * (n + 1) / 2;
  int nsq = n * n;
  Eigen::MatrixXd Dt (p, nsq);
  Dt.setZero ();
  Eigen::MatrixXd II (n, n);
  II.setIdentity ();
  Eigen::MatrixXd I (p, p);
  I.setIdentity ();
  Eigen::MatrixXd k (n, n);
  k.setZero ();
  int ii = 0;
  int jj = 0;
  for (int i = 0; i < p; ++i) {
    k (ii, jj) = i;
    ii++;
    if (ii == n) {
      ii = jj + 1;
      jj = ii;
    }
  }
  for (int j = 0; j < n; ++j) {
    for (int i = j; i < n; ++i) {
      Eigen::MatrixXd Iij = I.col (k (i, j));
      Eigen::VectorXd cc = II.row (j);
      Eigen::VectorXd rr = II.row (i);
      Eigen::MatrixXd Eij = (rr * cc.transpose ());
      Eigen::MatrixXd Eji = (cc * rr.transpose ());
      Eigen::MatrixXd Tij (Eij.rows (), Eij.cols ());
      if (i == j) {
        Tij = Eij;
      } else {
        Tij = Eij + Eji;
      }
      Tij.resize (1, n * n);
      Dt += (Iij * Tij);
    }
  }
  Eigen::MatrixXd DD = Dt.transpose ();
  return DD;
}

//' Half-vectorization indices
//'
//' @description
//' Returns row-column indices of the off-diagonal entries in `vech(X)`.
//'
//' @param k \[`integer(1)`\]\cr
//' Matrix dimension.
//'
//' @return
//' A numeric matrix with two columns.
//'
//' @keywords internal

Eigen::MatrixXd vechor (int k) {
  int m = ((k * k) - k) / 2 + k; //number of paramters in vech(X)
  Eigen::MatrixXd ind ((m - k), 2);
  ind.setZero ();
  int gg = 1; //indices for cpp
  int hh = 0; //indices for cpp
  for (int i = 0; i < (m - k); ++i) {
    ind (i, 0) = gg;
    ind (i, 1) = hh;
    gg++;
    if (gg == (k)) {
      hh = hh + 1;
      gg = hh + 1;
    }
  }
  return ind;
}

//' Half-vectorize a matrix
//'
//' @description
//' Returns the lower triangular part of a square matrix, including the
//' diagonal, stacked columnwise.
//'
//' @param Q \[`matrix`\]\cr
//' Square matrix.
//'
//' @return
//' A numeric vector.
//'
//' @keywords internal

Eigen::VectorXd vech (Eigen::MatrixXd Q) {
  int n = Q.rows ();
  Eigen::VectorXd vQ (n * (n + 1) / 2);
  vQ.setZero ();
  int cur = 0;
  for (int jj = 0; jj < n; jj++) {
    vQ.segment (cur, n - jj) = Q.block (jj, jj, n - jj, 1);
    cur = cur + n - jj;
  }
  return vQ;
}

//' Diagonal indices of a vectorized square matrix
//'
//' @description
//' Returns the positions of diagonal elements in the columnwise vectorization
//' of a `k` by `k` matrix.
//'
//' @param k \[`integer(1)`\]\cr
//' Matrix dimension.
//'
//' @return
//' A numeric vector of length `k`.
//'
//' @keywords internal

Eigen::VectorXd fdiag (int k) {
  Eigen::VectorXd ind (k);
  ind.setZero ();
  int kk = k;
  ind[0] = 0; //indices for cpp
  for (int i = 1; i < (k); i++) {
    ind[i] = ind[i - 1] + kk;
    kk = (kk - 1);
  }
  return ind;
}

//' Half-vectorization indices including the diagonal
//'
//' @description
//' Returns row-column indices of entries in `vech(X)`, including diagonal
//' entries.
//'
//' @param k \[`integer(1)`\]\cr
//' Matrix dimension.
//'
//' @return
//' A numeric matrix with two columns.
//'
//' @keywords internal

Eigen::MatrixXd vechor_diag (int k) {
  int m = ((k * k) - k) / 2 + k; //number of paramters in vech(X)
  Eigen::MatrixXd ind (m, 2);
  ind.setZero ();
  int gg = 0; //indices for cpp
  int hh = 0; //indices for cpp
  for (int i = 0; i < (m); ++i) {
    ind (i, 0) = gg;
    ind (i, 1) = hh;
    gg++;
    if (gg == (k)) {
      hh = hh + 1;
      gg = hh;
    }
  }
  return ind;
}

//' Commutation matrix
//'
//' @description
//' Computes the commutation matrix `K` such that `K vec(X) = vec(t(X))`.
//'
//' @param r \[`integer(1)`\]\cr
//' Number of rows in `X`.
//'
//' @param c \[`integer(1)`\]\cr
//' Number of columns in `X`.
//'
//' @return
//' A numeric commutation matrix.
//'
//' @keywords internal

Eigen::MatrixXd commmat (int r, int c) {
  int p = r * c;
  Eigen::MatrixXd K (p, p);
  K.setZero ();
  Eigen::MatrixXd Ir (r, r);
  Eigen::MatrixXd Ic (c, c);
  Ir.setIdentity ();
  Ic.setIdentity ();
  for (int i = 0; i < r; ++i) {
    for (int j = 0; j < c; ++j) {
      Eigen::VectorXd cc = Ic.row (j);
      Eigen::VectorXd rr = Ir.row (i);
      Eigen::MatrixXd Hij = rr * cc.transpose ();
      Eigen::MatrixXd tHij = Hij.transpose ();
      K += kroneckerProduct (Hij, tHij);
    }
  }
  return K;
}

//' Elimination matrix
//'
//' @description
//' Computes the elimination matrix `L` such that `L vec(X) = vech(X)`.
//'
//' @param n \[`integer(1)`\]\cr
//' Matrix dimension.
//'
//' @return
//' A numeric elimination matrix.
//'
//' @keywords internal

Eigen::MatrixXd elimmat (int n) {
  int p = n * (n + 1) / 2;
  int nsq = n * n;
  Eigen::MatrixXd L (p, nsq);
  L.setZero ();
  Eigen::MatrixXd I (p, p);
  I.setIdentity ();
  Eigen::MatrixXd II (n, n);
  II.setIdentity ();
  Eigen::MatrixXd k (n, n);
  k.setZero ();
  int ii = 0;
  int jj = 0;
  for (int i = 0; i < p; ++i) {
    k (ii, jj) = i;
    ii++;
    if (ii == n) {
      ii = jj + 1;
      jj = ii;
    }
  }
  for (int j = 0; j < n; ++j) {
    for (int i = j; i < n; ++i) {
      Eigen::MatrixXd Iij = I.col (k (i, j));
      Eigen::VectorXd cc = II.row (j);
      Eigen::VectorXd rr = II.row (i);
      Eigen::MatrixXd IIij = (rr * cc.transpose ());
      IIij.resize (1, n * n);
      L += (Iij * IIij);
    }
  }
  return L;
}

//' Remove a matrix row
//'
//' @description
//' Removes one row from an Eigen matrix in place.
//'
//' @param matrix \[`matrix`\]\cr
//' Matrix to modify.
//'
//' @param rowToRemove \[`integer(1)`\]\cr
//' Zero-based row index to remove.
//'
//' @return
//' No return value. The matrix is modified in place.
//'
//' @keywords internal

void removeRow (Eigen::MatrixXd &matrix, unsigned int rowToRemove) {
  unsigned int numRows = matrix.rows () - 1;
  unsigned int numCols = matrix.cols ();
  if (rowToRemove < numRows) {
    matrix.block (rowToRemove, 0, numRows - rowToRemove, numCols) =
        matrix.bottomRows (numRows - rowToRemove);
  }
  matrix.conservativeResize (numRows, numCols);
}

//' Remove a matrix column
//'
//' @description
//' Removes one column from an Eigen matrix in place.
//'
//' @param matrix \[`matrix`\]\cr
//' Matrix to modify.
//'
//' @param colToRemove \[`integer(1)`\]\cr
//' Zero-based column index to remove.
//'
//' @return
//' No return value. The matrix is modified in place.
//'
//' @keywords internal

void removeColumn (Eigen::MatrixXd &matrix, unsigned int colToRemove) {
  unsigned int numRows = matrix.rows ();
  unsigned int numCols = matrix.cols () - 1;
  if (colToRemove < numCols) {
    matrix.block (0, colToRemove, numRows, numCols - colToRemove) =
        matrix.rightCols (numCols - colToRemove);
  }
  matrix.conservativeResize (numRows, numCols);
}

//' Non-missing row indices
//'
//' @description
//' Returns the indices of rows that do not contain any `NA` entries.
//'
//' @param Xc \[`matrix`\]\cr
//' Matrix to scan.
//'
//' @return
//' An integer vector of zero-based row indices.
//'
//' @keywords internal

Rcpp::IntegerVector NA_ind (Eigen::MatrixXd Xc) {
  Rcpp::IntegerVector vecNA (Xc.rows ());

  int j = 0;
  for (int i = 0; i < vecNA.size (); i++) {
    bool complete = true;
    for (int jc = 0; jc < Xc.cols (); jc++) {
      if (R_IsNA (Xc (i, jc))) {
        complete = false;
        break;
      }
    }
    if (complete) {
      vecNA[j] = i;
      j++;
    }
  }

  if (j == 0) {
    return Rcpp::IntegerVector::create (-1);
  }
  return vecNA[Rcpp::Range (0, j - 1)];
}
