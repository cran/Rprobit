#include "Rcpp.h"

#include <cmath>

#include "distrib.h"
#include "lin_alg.h"
#include "toms462.h"

// Solow-Joe covariance approximation.
//
// Computes the covariance matrix used in the Solow-Joe regression
// approximation of conditional probabilities.
//
// x: Evaluation point.
// r: Correlation matrix.
//
// Returns a numeric covariance approximation matrix.

Eigen::MatrixXd Omega (Eigen::VectorXd x, Eigen::MatrixXd r) {
  // save the number of rows of x in m
  int m = x.rows ();

  // Define Vector P such that
  // P(k) = Pr(X_k < x_k)

  // reserve memory for the vector P of length m with double entries
  Eigen::VectorXd P (m);
  for (int k = 0; k < m; k++) {
    // go through all entries k of P
    // compute the standard normal cumulative distribution function at x=x(k)
    // and save it in P(k)
    P (k) = std_normal_cdf (x (k));
  }

  // Define Matrix z such that
  // z(k,l) = Pr(X_k < x_k, X_l < x_l)
  // z(k,k) = Pr(X_k < x_k)

  // reserve memory for an m-by-m matrix of doubles and call it z
  Eigen::MatrixXd z (m, m);
  // set all entries of z to zero
  z.setZero ();
  for (int k = 0; k < (m); k++) {
    // go trough the rows k of z
    for (int l = (k + 1); l < (m); l++) {
      // go trough the elements l of the row k which are past the diagonal

      // save element r(k,l) in rrkl
      double rrkl = r (k, l);
      if (r (k, l) > 1 - tol) {
        // if r(k,l)>1
        // set rrkl = 0.99
        rrkl = 1 - tol;
      }
      if (r (k, l) < tol - 1) {
        // if r(k,l)<-1
        // set rrkl=-0.99
        rrkl = tol - 1;
      }

      // computes the probability for two normal variates X and Y whose
      // correlation is rrkl, that -x(k) <= X and -x(l) <= Y.
      z (k, l) = bivnor (-x (k), -x (l), rrkl);
      // make the matrix symmetric
      z (l, k) = z (k, l);
    }
  }

  // reserve memory for an m-by-m-Matrix of doubles called tt
  Eigen::MatrixXd tt (m, m);
  // make tt the m-by-m-Zero-Matrix
  tt.setZero ();
  // set the diagonal entries of tt to the entries of P (cdf-values of x(k))
  tt.diagonal () = P;
  // add to z (the above defined matrix of bivariate probabilities) the matrix
  // tt of diagonal univariate probabilities
  z = z + tt;

  // Define Matrix z1 such that
  // z1(k,l) = Pr(X_k < x(k)) * Pr(X_l < x(l))

  // reserve memory for an m-by-m-Matrix of doubles called z1
  Eigen::MatrixXd z1 (m, m);
  // define z1(k,l) = P(k)*P(l)
  z1 = P * P.transpose ();

  // calculate Omega
  Eigen::MatrixXd Omega (m, m);

  Omega = z - z1;

  // return the matrix Omega
  return Omega;
}

//' Solow-Joe approximation
//'
//' @description
//' Computes the log CDF, gradient, and Hessian for a multivariate Gaussian
//' distribution using the Solow-Joe approximation.
//'
//' @param x \[`numeric()`\]\cr
//' Evaluation point.
//'
//' @param r \[`matrix`\]\cr
//' Correlation matrix.
//'
//' @return
//' \describe{
//'   \item{`SJ(x, r)`}{Log-CDF value.}
//'   \item{`SJ_grad(x, r)`}{
//'   Gradient vector wrt evaluation point and correlations.
//'   Last entry is the log-CDF value.
//'   }
//'   \item{`SJ_hess(x, r)`}{Hessian matrix.}
//' }
//'
//' @examples
//' x <- c(0, 0)
//' r <- diag(2)
//' exp(SJ(x, r))
//' SJ_grad(x, r) # gradient + log-CDF value
//' SJ_hess(x, r) # Hessian matrix
//'
//' @export
//'
// [[Rcpp::export]]

double SJ (Eigen::VectorXd x, Eigen::MatrixXd r) {
  double lcond = 0;

  // save the number of rows of x in m
  int m = x.rows ();

  Eigen::VectorXd tol_which (m);
  tol_which.setZero ();

  // reserve memory for a vector of doubles of length m called zeroes
  Eigen::VectorXd zeros (m);
  // set all entries of zeroes to 0
  zeros.setZero ();

  // reserve memory for an m-by-m-Matrix of doubles called Mones
  Eigen::MatrixXd Mones (m, m);
  // set all entries of the m-by-m-Matrix Mones to 1
  Mones.setOnes ();

  // reserve memory for an m-by-m-Matrix of doubles called Mones
  Eigen::MatrixXd Mzeros (m, m);
  // set all entries of the m-by-m-Matrix Mones to 1
  Mzeros.setZero ();

  // reserve memory for vector of length m called ones with double entries
  Eigen::VectorXd ones (m);
  // set all entries of ones to 1
  ones.setOnes ();

  // Define Vector P such that P(k) = Pr(X_k < x_k)
  // reserve memory for the vector P of length m with double entries
  Eigen::VectorXd P (m);
  P.setZero ();

  for (int k = 0; k < m; k++) {
    // go through all entries k of P
    // compute the sandard normal cumulative distribution function at x=x(k)
    // and save it in P(k)
    P (k) = std_normal_cdf (x (k));
  }

  // calculate PDF for later usage.
  // reserve memory for the vector P of length m with double entries
  Eigen::VectorXd p (m);
  p.setZero ();

  for (int k = 0; k < m; k++) {
    // go through all entries k of P
    // compute the sandard normal cumulative distribution function at x=x(k)
    // and save it in P(k)
    p (k) = std_normal_pdf (x (k));
  }

  if (m > 1) {
    // get Omega and its derivatives
    Eigen::MatrixXd Om (m, m);
    Om.setZero ();
    Om = Omega (x, r);

    // start calc with first bivariate CDF
    // reserve memory for an m-by-m matrix of doubles and call it z
    Eigen::MatrixXd z (m, m);
    // set all entries of z to zero
    z.setZero ();
    for (int k = 0; k < (m); k++) {
      // go trough the rows k of z
      for (int l = (k + 1); l < (m); l++) {
        // go trough the elements l of the row k which are past the diagonal

        // save element r(k,l) in rrkl
        double rrkl = r (k, l);
        if (r (k, l) > 1 - tol) {
          // if r(k,l)>1
          // set rrkl = 0.99
          rrkl = 1 - tol;
        }
        if (r (k, l) < tol - 1) {
          // if r(k,l)<-1
          // set rrkl=-0.99
          rrkl = tol - 1;
        }

        // computes the probability for two normal variates X and Y whose
        // correlation is rrkl, that -x(k) <= X and -x(l) <= Y.
        z (k, l) = bivnor (-x (k), -x (l), rrkl);
        // make the matrix symmetric
        z (l, k) = z (k, l);
      }
    }

    // reserve memory for an m-by-m-Matrix of doubles called tt
    Eigen::MatrixXd tt (m, m);
    // make tt the m-by-m-Zero-Matrix
    tt.setZero ();
    // set the diagonal entries of tt to the entries of P (cdf-values of x(k))
    tt.diagonal () = P;
    // add to z (the above defined matrix of bivariate probabilities) the
    // matrix tt of diagonal univariate probabilities
    z = z + tt;

    if (z (0, 1) > tol) {
      // if z(0,1) > tol ( z(0,1) = biv_normal_cdf(-x(0), -x(1), r(k,l)) ),
      // (0,1) is
      // first off-diagonal element
      // set lcond = std::log (z(0,1)) in MATLAB 1,2 --> 0,1, so first row,
      // second
      // element
      lcond = std::log (z (0, 1));
    } else {
      // if the first off diagonal element entry is less than toe predefined
      // tollerance
      // set lcond = std::log (tol)
      lcond = std::log (tol);
      tol_which (0) = 1.0;
    }
    Eigen::VectorXd ones_m_P (m);
    ones_m_P.setZero ();
    ones_m_P = ones - P;

    // carry on with iterating over k
    for (int k = 2; k < m; k++) {
      // initialize at k=2, since first two are dealt with.
      Eigen::VectorXd omega12 (k);
      omega12.setZero ();
      omega12 = Om.block (0, k, k, 1);

      Eigen::MatrixXd omega11 (k, k);
      omega11.setZero ();
      omega11 = Om.block (0, 0, k, k);

      Eigen::MatrixXd iom (k, k);
      iom.setZero ();
      iom = omega11.inverse ();
      double condk = 0;
      double lcondk = 0;

      Eigen::VectorXd he_v (k);
      he_v.setZero ();
      he_v = iom * ones_m_P.block (0, 0, k, 1);

      condk = P (k) + omega12.dot (he_v);

      if (condk > tol) {
        // check if condk = std_normal_cdf(x(k)) + tempcon is greater than the
        // tollarance
        // if so, lcondk = std::log (condk)
        lcondk = std::log (condk);
      } else {
        // else, lcondk = std::log (tol)
        lcondk = std::log (tol);
        tol_which (k) = 1.0 * k + 1.0;
      }
      // Update lcond, which will give the approximation of
      // std::log ( Pr(X_1 < x(1), X_2 < x(2) ) * \prod_{k=3}^m Pr(X_k < x(k) |
      // X_1 < x(1), ..., X_{k-1} < x(k-1) ) )
      // This will be one part of the output!

      lcond = lcond + lcondk;
      // lcond = lcondk;
    }
  } else {
    lcond = std::log (P (0));
  }
  return lcond;
}

// Gradient of the Solow-Joe covariance approximation.
//
// x: Evaluation point.
// r: Correlation matrix.
//
// Returns a list of derivative matrices.

Rcpp::List dOmega (Eigen::VectorXd x, Eigen::MatrixXd r) {
  // save the number of rows of x in m
  int m = x.rows ();

  // Define Vector zeros such that
  // ones(k) = 0

  // reserve memory for a vector of doubles of length m called zeroes
  Eigen::VectorXd zeros (m);
  // set all entries of zeroes to 0
  zeros.setZero ();

  // Define Matrix Mones such that
  // Mones(k) = 1

  // reserve memory for an m-by-m-Matrix of doubles called Mones
  Eigen::MatrixXd Mones (m, m);
  // set all entries of the m-by-m-Matrix Mones to 1
  Mones.setOnes ();

  // Define Matrix Mones such that
  // Mones(k) = 1

  // reserve memory for an m-by-m-Matrix of doubles called Mones
  Eigen::MatrixXd Mzeros (m, m);
  // set all entries of the m-by-m-Matrix Mones to 1
  Mzeros.setZero ();

  // Define Vector ones such that
  // ones(k) = 1

  // reserve memory for vector of length m called ones with double entries
  Eigen::VectorXd ones (m);
  // set all entries of ones to 1
  ones.setOnes ();

  // Define Vector P such that
  // P(k) = Pr(X_k < x_k)

  // reserve memory for the vector P of length m with double entries
  Eigen::VectorXd P (m);
  for (int k = 0; k < m; k++) {
    // go through all entries k of P
    // compute the sandard normal cumulative distribution function at x=x(k)
    // and save it in P(k)
    P (k) = std_normal_cdf (x (k));
  }

  // calculate PDF for later usage.
  // reserve memory for the vector P of length m with double entries
  Eigen::VectorXd p (m);
  for (int k = 0; k < m; k++) {
    // go through all entries k of P
    // compute the sandard normal cumulative distribution function at x=x(k)
    // and save it in P(k)
    p (k) = std_normal_pdf (x (k));
  }

  // Define Matrix z such that
  // z(k,l) = Pr(X_k < x_k, X_l < x_l)
  // z(k,k) = Pr(X_k < x_k)

  // reserve memory for an m-by-m matrix of doubles and call it z
  Eigen::MatrixXd z (m, m);
  // set all entries of z to zero
  z.setZero ();
  for (int k = 0; k < (m); k++) {
    // go trough the rows k of z
    for (int l = (k + 1); l < (m); l++) {
      // go trough the elements l of the row k which are past the diagonal

      // save element r(k,l) in rrkl
      double rrkl = r (k, l);
      if (r (k, l) > 1 - tol) {
        // if r(k,l)>1
        // set rrkl = 0.99
        rrkl = 1 - tol;
      }
      if (r (k, l) < tol - 1) {
        // if r(k,l)<-1
        // set rrkl=-0.99
        rrkl = tol - 1;
      }

      // computes the probability for two normal variates X and Y whose
      // correlation is rrkl, that -x(k) <= X and -x(l) <= Y.
      z (k, l) = bivnor (-x (k), -x (l), rrkl);
      // make the matrix symmetric
      z (l, k) = z (k, l);
    }
  }

  // reserve memory for an m-by-m-Matrix of doubles called tt
  Eigen::MatrixXd tt (m, m);
  // make tt the m-by-m-Zero-Matrix
  tt.setZero ();
  // set the diagonal entries of tt to the entries of P (cdf-values of x(k))
  tt.diagonal () = P;
  // add to z (the above defined matrix of bivariate probabilities) the matrix
  // tt of diagonal univariate probabilities
  z = z + tt;

  // Define Matrix z1 such that
  // z1(k,l) = Pr(X_k < x(k)) * Pr(X_l < x(l))

  // reserve memory for an m-by-m-Matrix of doubles called z1
  Eigen::MatrixXd z1 (m, m);
  // define z1(k,l) = P(k)*P(l)
  z1 = P * P.transpose ();

  // calculate Omega
  Eigen::MatrixXd Omega (m, m);

  Omega = z - z1;

  // go for derivatives

  // Define Matrix rho1 such that
  // rho1(k,l) = cor(X_k, X_l) = r(k,l) for k!=l
  // rho1(k,k) = 0

  // reserve memory for an m-by-m-Matrix of doubles called rho1
  Eigen::MatrixXd rho1 (m, m);
  // set rho1 equal to r
  rho1 = r;
  // set the diagonal entries of rho1 to 0
  rho1.diagonal () = zeros;

  // Define Matrix rho2 such that
  // rho2(k,l) = sqrt(1-r(k,l)^2) for k!=l
  // rho2(k,k) = 1

  // reserve memory for an m-by-m-Matrix of doubles called rho2
  Eigen::MatrixXd rho2 (m, m);
  rho2 =
    // calculate elementwise rho2 = sqrt(1-rho1^2) and interpret it as a
    // matrix
    (Mones.array () - rho1.array () * rho1.array ()).array ().sqrt ().matrix ();

  // Define Matrix muygx such that
  // muygx(k,l) = (x(l) - x(k) * r(k,l)) / (sqrt(1-r(k,l)^2)) for k!=l
  // muygx(k,k) = x(k)*(1-r(k,k)) = 0

  // reserve memory for an m-by-m-Matrix of doubles called muygx
  Eigen::MatrixXd muygx (m, m);
  muygx = (((ones * x.transpose ()
               - (r.array () * (x * ones.transpose ()).array ()).matrix ()).array ())
             // diagonal entries of muygx(i,i) = (1-r(i,i))*x(i)=0 and other elements
             // muygx(i,j) = (x(j)-x(i)*r(i,j))/(sqrt(1-r(i,j)^2))
             / rho2.array ()).matrix ();

             // Define Matrix nmuygx such that
             // nmuygx(k,l) = Pr(X_k < (x(l)-x(k)*r(k,l))/(sqrt(1-r(k,l)^2)) ) for k!=l
             // nmuygx(k,k) = Pr(X_k < x(k)*(1-r(k,k)) ) = 1/2

             // reserve memory for an m-by-m-Matrix of doubles called nmuygx
             Eigen::MatrixXd nmuygx (m, m);
             for (int k = 0; k < m; k++) {
               // for each row k of nmuygx
               for (int l = 0; l < m; l++) {
                 // for each element l of row k of nmuygx
                 nmuygx (k, l) = std_normal_cdf (muygx (k, l));
               }
             }

             // Define Matrix pdfmuygx such that
             // pdfmuygx(k,l) = \phi( (x(l) - x(k)*r(k,l))/(sqrt(1-r(k,l)^2)) ) for k!=l
             // pdfmuygx(k,k) = \phi( x(k)*(1-r(k,k)) ) = \phi( 0 )

             // reserve memory for an m-by-m-Matrix of doubles called pdfmuygx
             Eigen::MatrixXd pdfmuygx (m, m);
             for (int k = 0; k < m; k++) {
               // for each row k of pdfmuygx
               for (int l = 0; l < m; l++) {
                 // for each element l of row k of pdfmuygx
                 //
                 pdfmuygx (k, l) = std_normal_pdf (muygx (k, l));
               }
             }

             // Define Matrix pdiag such that
             // pdiag(k,k) = \phi(x(k))
             // pdiag(k,l) = 0 for k!=l

             // reserve memory for an m-by-m-Matrix of doubles called pdiag
             Eigen::MatrixXd pdiag (m, m);
             // set all elements of pdiag to 0
             pdiag.setZero ();
             // set the diagonal elements of pdiag to p = std_normal_pdf(x(k))
             pdiag.diagonal () = p;

             // Define Matrix dPhixy such that
             // dPhixy(k,l) = \phi(x(k)) * Pr(X < muygx(k,l))

             // reserve memory for an m-by-m-Matrix of doubles called dPhixy
             Eigen::MatrixXd dPhixy (m, m);
             // set dPhixy(i,j) = std_normal_pdf(x(i)) * std_normal_cdf(muygx(i,j))
             dPhixy = pdiag * nmuygx;

             // Define Matrix phiPhi such that
             // phiPhi(k,l) = \phi(x(k)) * Pr(X_l < x(l) )

             // reserve memory for an m-by-m-Matrix of doubles called phiPhi
             Eigen::MatrixXd phiPhi (m, m);
             // set phiPhi(i,j) = std_normal_pdf(x(i)) * std_normal_cdf(x(j))
             phiPhi = p * P.transpose ();

             // Define Matrix dCov such that
             // dCov(k,l) = dPhixy(k,l) - phiPhi(k,l)
             // = \phi(x(k)) * Pr(X < muygx(k,l)) - \phi(x(k)) * Pr(X_l < x(l))
             // = \phi(x(k)) * ( Pr(X < muygx(k,l)) - Pr(X_l < x(l)) ) for k!=l
             // dCov(k,k) = \phi(x(k)) * (1 - 2 * Pr(X_k < x(k)) ) updated later

             // reserve memory for an m-by-m-Matrix of doubles called dCov
             Eigen::MatrixXd dCov (m, m);
             // set dCov(i,j) = std_normal_pdf(x(i)) * ( std_normal_cdf(muygx(i,j)) -
             // std_normal_cdf(x(j)) )
             dCov = phiPhi - dPhixy;

             // Define Matrix dCovDiag such that
             // dCovDiag(k,k) = \phi(x(k)) * ( Pr(X < muygx(k,k)) - Pr(X_k < x(k)) )

             // reserve memory for an m-by-m-Matrix of doubles called dCovDiag
             Eigen::MatrixXd dCovDiag (m, m);
             // set all elements of dCovDiag to 0
             dCovDiag.setZero ();
             // set diagonal elements of dCovDiag to diagonal elements of dCov
             dCovDiag.diagonal () = dCov.diagonal ();

             // Define Matrix dCovDiag1 such that
             // dCovDiag1(k,k) = \phi(x(k)) * ( Pr(X < muygx(k,k)) - Pr(X_k < x(k)) )

             // reserve memory for an m-by-m-Matrix of doubles called dCovDiag1
             Eigen::MatrixXd dCovDiag1 (m, m);
             // set all elements of dCovDiag1 to 0
             dCovDiag1.setZero ();
             dCovDiag1.diagonal () = p.array ()
               // set diagonal elements of dCovDiag1 to dCovDiag1(i,i) =
               // std_normal_pdf(x(i)) * (1 - 2 * std_normal_cdf(x(i)))
               * (ones.array () - (ones.array () * 2) * P.array ());

             // Update Matrix dCov such that the diagonal is
             // dCov(k,k) = \phi(x(k)) * (1 - 2 * Pr(X_k < x(k)) )
             // set the diagonal elements of dCov to std_normal_pdf(x(i)) * (1 - 2 *
             // std_normal_cdf(x(i)))
             dCov = dCov - dCovDiag + dCovDiag1;
             // dCov(k,l) = d/(d x(k)) Cov(I_k,I_l)

             // derivatives with respect to correlations.
             // reserve memory for an m-by-m-Matrix of doubles called dPhi2
             Eigen::MatrixXd dPhi2 (m, m);
             // set dPhi2(i,j) = std_normal_pdf(x(i)) * std_normal_pdf(muygx(i,j)) /
             // sqrt(1-rho1(i,j)^2) with rho1(i,j) = r(i,j)
             dPhi2 = ((pdiag * pdfmuygx).array () / rho2.array ()).matrix ();
             // set the diagonal elements of dPhi2 to 0
             dPhi2.diagonal () = zeros;

             // calculate dOmega
             // number of derivatives
             int nd;
             nd = m + m * (m - 1) / 2;
             Rcpp::List dOmega (nd);
             // derivatives with respect to kk1
             // cycle over kk1 entries
             for (int k = 0; k < (m); k++) {
               Eigen::MatrixXd dom = Mzeros;
               if ((x (k) > -bound_max) && (x (k) < bound_max)) {
                 // adjust k-th row of gradient
                 for (int j = 0; j < (m); j++) {
                   dom (k, j) = dPhixy (k, j) - p (k) * P (j);
                   dom (j, k) = dPhixy (k, j) - p (k) * P (j);
                 }
                 dom (k, k) = p (k) * (1 - 2 * P (k));
               }
               dOmega[k] = dom;
             }

             // derivatives with respect to correlation entries.
             Eigen::VectorXd colI = trilCols (m);
             Eigen::VectorXd rowI = trilRows (m);
             // cycle over entries
             for (int j = 0; j < nd - m; j++) {
               int rc = rowI (j);
               int cc = colI (j);
               // corr r_{i,j} influences only element i,j.
               Eigen::MatrixXd dom = Mzeros;
               if ((r (rc, cc) > tol - 1) && (r (rc, cc) < 1 - tol)) {
                 dom (rc, cc) = dPhi2 (rc, cc);
                 dom (cc, rc) = dom (rc, cc);
               }
               dOmega[m + j] = dom;
             }

             // return the matrix Omega
             return dOmega;
}

//' @describeIn SJ Gradient computation
//' @export
//'
// [[Rcpp::export]]

Eigen::VectorXd SJ_grad (Eigen::VectorXd x, Eigen::MatrixXd r) {
  double lcond = 0;
  // save the number of rows of x in m
  int m = x.rows ();

  Eigen::VectorXd tol_which (m);
  tol_which.setZero ();

  // reserve memory for a vector of doubles of length m called zeroes
  Eigen::VectorXd zeros (m);
  // set all entries of zeroes to 0
  zeros.setZero ();

  // reserve memory for an m-by-m-Matrix of doubles called Mones
  Eigen::MatrixXd Mones (m, m);
  // set all entries of the m-by-m-Matrix Mones to 1
  Mones.setOnes ();

  // reserve memory for an m-by-m-Matrix of doubles called Mones
  Eigen::MatrixXd Mzeros (m, m);
  // set all entries of the m-by-m-Matrix Mones to 1
  Mzeros.setZero ();

  // reserve memory for vector of length m called ones with double entries
  Eigen::VectorXd ones (m);
  // set all entries of ones to 1
  ones.setOnes ();

  // Define Vector P such that P(k) = Pr(X_k < x_k)

  // reserve memory for the vector P of length m with double entries
  Eigen::VectorXd P (m);
  P.setZero ();

  for (int k = 0; k < m; k++) {
    // go through all entries k of P
    // compute the standard normal cumulative distribution function at x=x(k)
    // and save it in P(k)
    P (k) = std_normal_cdf (x (k));
  }

  // calculate PDF for later usage.
  // reserve memory for the vector P of length m with double entries
  Eigen::VectorXd p (m);
  p.setZero ();
  for (int k = 0; k < m; k++) {
    // go through all entries k of P
    // compute the standard normal cumulative distribution function at x=x(k)
    // and save it in P(k)
    p (k) = std_normal_pdf (x (k));
  }

  // Define Matrix z such that
  // z(k,l) = Pr(X_k < x_k, X_l < x_l)
  // z(k,k) = Pr(X_k < x_k)

  // calculate Omega
  Eigen::MatrixXd Om (m, m);
  Om.setZero ();
  Om = Omega (x, r);

  // define gradient and start calculations.
  int nd;
  nd = m + m * (m - 1) / 2;
  Eigen::VectorXd dlcond (nd);
  dlcond.setZero ();

  if (m > 1) {

    Rcpp::List dom (nd);
    dom = dOmega (x, r);

    // calc initials
    // reserve memory for an m-by-m matrix of doubles and call it z
    Eigen::MatrixXd z (m, m);
    // set all entries of z to zero
    z.setZero ();
    for (int k = 0; k < (m); k++) {
      // go trough the rows k of z
      for (int l = (k + 1); l < (m); l++) {
        // go trough the elements l of the row k which are past the diagonal

        // save element r(k,l) in rrkl
        double rrkl = r (k, l);
        if (r (k, l) > 1 - tol) {
          // if r(k,l)>1
          // set rrkl = 0.99
          rrkl = 1 - tol;
        }
        if (r (k, l) < tol - 1) {
          // if r(k,l)<-1
          // set rrkl=-0.99
          rrkl = tol - 1;
        }

        // computes the probability for two normal variates X and Y whose
        // correlation is rrkl, that -x(k) <= X and -x(l) <= Y.
        z (k, l) = bivnor (-x (k), -x (l), rrkl);
        // make the matrix symmetric
        z (l, k) = z (k, l);
      }
    }

    // reserve memory for an m-by-m-Matrix of doubles called tt
    Eigen::MatrixXd tt (m, m);
    // make tt the m-by-m-Zero-Matrix
    tt.setZero ();
    // set the diagonal entries of tt to the entries of P (cdf-values of x(k))
    tt.diagonal () = P;
    // add to z (the above defined matrix of bivariate probabilities) the
    // matrix tt of diagonal univariate probabilities
    z = z + tt;

    if (z (0, 1) > tol) {
      // if z(0,1) > tol ( z(0,1) = biv_normal_cdf (-x(0), -x(1), r(k,l)) ),
      // (0,1)
      // is first off-diagonal element
      // set lcond = std::log (z(0,1)) in MATLAB 1,2 --> 0,1, so first row,
      // second
      // element
      lcond = std::log (z (0, 1));
    } else {
      // if the first off diagonal element entry is less than toe predefined
      // tollerance
      // set lcond = std::log (tol)
      lcond = std::log (tol);
    }

    Eigen::VectorXd ones_m_P (m);
    ones_m_P.setZero ();
    ones_m_P = ones - P;

    for (int cur = 0; cur < m; cur++) {
      // cycle over derivatives with respect to xk
      Eigen::MatrixXd dom_cur (m, m);
      //
      dom_cur.setZero ();
      dom_cur = dom (cur);

      if (cur == 0) {
        // first entry -> derive
        double r12 = r (1, 0);
        if (z (0, 1) > tol) {
          dlcond (0) = dlcond (0)
          + (std_normal_pdf (x (0))
               * normal_cdf (x (1) - x (0) * r12, 1 - r12 * r12))
               / z (0, 1);
        }
      }
      if (cur == 1) {
        // second entry -> derive
        double r12 = r (1, 0);
        if (z (0, 1) > tol) {
          dlcond (1) = dlcond (1)
          + (std_normal_pdf (x (1))
               * normal_cdf (x (0) - x (1) * r12, 1 - r12 * r12))
               / z (0, 1);
        }
      }

      for (int k = 2; k < m; k++) {
        // k runs to calculate condk.
        Eigen::VectorXd omega12 (k);
        omega12.setZero ();
        omega12 = Om.block (0, k, k, 1);
        Eigen::MatrixXd omega11 (k, k);
        omega11.setZero ();
        omega11 = Om.block (0, 0, k, k);
        Eigen::MatrixXd iom (k, k);
        iom.setZero ();
        iom = omega11.inverse ();
        double condk = 0;

        Eigen::VectorXd he_v (k);
        he_v.setZero ();
        he_v = iom * ones_m_P.block (0, 0, k, 1);

        Eigen::VectorXd oio (k);
        oio.setZero ();
        oio = iom * omega12;

        condk = P (k) + omega12.dot (he_v);
        if (condk > tol) {
          // check whether condk = Phi(x(k)) + tempcon is greater than tol
          if (cur == 0) {
            lcond += std::log (condk);
          }
          // no truncation -< contribution to derivative // if so, lcondk =
          // std::log
          // (condk)
          double dcondk = 0;
          if (k == cur) {
            dcondk = std_normal_pdf (x (k));
          }
          Eigen::VectorXd dcb1 (k);
          dcb1.setZero ();
          dcb1 = dom_cur.block (0, k, k, 1);

          Eigen::MatrixXd dcb2 (k, k);
          dcb2.setZero ();
          dcb2 = dom_cur.block (0, 0, k, k);

          dcondk = dcondk + dcb1.dot (he_v);
          dcondk = dcondk - (oio.transpose () * (dcb2 * he_v));
          if (k > cur) {
            dcondk = dcondk - oio (cur) * p (cur);
          }
          dcondk = dcondk / condk;
          dlcond (cur) = dcondk + dlcond (cur);
        } else {
          if (cur == 0) {
            lcond += std::log (tol);
          }
        }

      }
    }
    // after the first two, the initial contrib does not matter.
    // same for derivs with respect to rij
    for (int cur = 0; cur < nd - m; cur++) {
      // cycle over derivatives with respect to xk
      Eigen::MatrixXd dom_cur (m, m);
      dom_cur.setZero ();

      // derivatives after the derivs w.r.t. xk.
      dom_cur = dom (cur + m);
      if (cur == 0) {
        // first entry -> this appears in initial term.
        double r12 = r (1, 0);
        if (z (0, 1) > tol) {
          dlcond (m) = (std_normal_pdf (x (0))
                          * normal_pdf (x (1) - x (0) * r12, 1 - r12 * r12))
          / z (0, 1);
        }
      }

      // after the first iteration continue as above.
      for (int k = 2; k < m; k++) {
        // k runs to calculate condk.
        Eigen::VectorXd omega12 (k);
        omega12.setZero ();
        omega12 = Om.block (0, k, k, 1);
        Eigen::MatrixXd omega11 (k, k);
        omega11.setZero ();
        omega11 = Om.block (0, 0, k, k);
        Eigen::MatrixXd iom (k, k);
        iom.setZero ();
        iom = omega11.inverse ();

        double condk;
        condk = 0;

        Eigen::VectorXd he_v (k);
        he_v.setZero ();
        he_v = iom * ones_m_P.block (0, 0, k, 1);

        Eigen::VectorXd oio (k);
        oio.setZero ();
        oio = iom * omega12;

        condk = P (k) + omega12.dot (he_v);
        if (condk > tol) {
          // check whether condk = Phi(x(k)) + tempcon is greater than tol
          // no truncation -< contribution to derivative // if so, lcondk =
          // std::log(condk)
          double dcondk = 0;

          Eigen::VectorXd dcb1 (k);
          dcb1.setZero ();
          dcb1 = dom_cur.block (0, k, k, 1);

          Eigen::MatrixXd dcb2 (k, k);
          dcb2.setZero ();
          dcb2 = dom_cur.block (0, 0, k, k);

          dcondk = dcb1.dot (he_v);
          dcondk = dcondk - (oio.transpose () * (dcb2 * he_v));

          dlcond (cur + m) = dcondk / condk + dlcond (cur + m);

        }

      }
    }
  } else {
    dlcond (0) = p (0) / P (0);
    lcond = std::log (P[0]);
  }

  // add log of approximated probability to the gradient vector

  Eigen::VectorXd output_final (1 + dlcond.size ());
  output_final << dlcond, lcond;
  return output_final;
}

// Hessian of the Solow-Joe covariance approximation.
//
// Computes second derivatives of the Solow-Joe covariance approximation.
//
// x: Evaluation point.
// r: Correlation matrix.
//
// Returns a nested list of derivative matrices.

Rcpp::List d2Omega (Eigen::VectorXd x, Eigen::MatrixXd r) {
  // save the number of rows of x in m
  int m = x.rows ();

  // reserve memory for a vector of doubles of length m called zeroes
  Eigen::VectorXd zeros (m);
  // set all entries of zeroes to 0
  zeros.setZero ();

  // reserve memory for an m-by-m-Matrix of doubles called Mones
  Eigen::MatrixXd Mones (m, m);
  // set all entries of the m-by-m-Matrix Mones to 1
  Mones.setOnes ();

  // reserve memory for an m-by-m-Matrix of doubles called Mones
  Eigen::MatrixXd Mzeros (m, m);
  // set all entries of the m-by-m-Matrix Mones to 1
  Mzeros.setZero ();

  // reserve memory for vector of length m called ones with double entries
  Eigen::VectorXd ones (m);
  // set all entries of ones to 1
  ones.setOnes ();

  // Define Vector P such that P(k) = Pr(X_k < x_k)
  // reserve memory for the vector P of length m with double entries
  Eigen::VectorXd P (m);
  for (int k = 0; k < m; k++) {
    // go through all entries k of P
    // compute the sandard normal cumulative distribution function at x=x(k)
    // and save it in P(k)
    P (k) = std_normal_cdf (x (k));
  }

  // calculate PDF for later usage.
  // reserve memory for the vector P of length m with double entries
  Eigen::VectorXd p (m);
  for (int k = 0; k < m; k++) {
    // go through all entries k of P
    // compute the sandard normal cumulative distribution function at x=x(k)
    // and save it in P(k)
    p (k) = std_normal_pdf (x (k));
  }

  // / include the dOmega part in order to obtain the relevant expressions!!

  // Define Matrix z such that
  // z(k,l) = Pr(X_k < x_k, X_l < x_l)
  // z(k,k) = Pr(X_k < x_k)

  // reserve memory for an m-by-m matrix of doubles and call it z
  Eigen::MatrixXd z (m, m);
  // set all entries of z to zero
  z.setZero ();
  for (int k = 0; k < (m); k++) {
    // go trough the rows k of z
    for (int l = (k + 1); l < (m); l++) {
      // go trough the elements l of the row k which are past the diagonal

      // save element r(k,l) in rrkl
      double rrkl = r (k, l);
      if (r (k, l) > 1 - tol) {
        // if r(k,l)>1
        // set rrkl = 0.99
        rrkl = 1 - tol;
      }
      if (r (k, l) < tol - 1) {
        // if r(k,l)<-1
        // set rrkl=-0.99
        rrkl = 1 - tol;
      }

      // computes the probability for two normal variates X and Y whose
      // correlation is rrkl, that -x(k) <= X and -x(l) <= Y.
      z (k, l) = bivnor (-x (k), -x (l), rrkl);
      // make the matrix symmetric
      z (l, k) = z (k, l);
    }
  }

  // reserve memory for an m-by-m-Matrix of doubles called tt
  Eigen::MatrixXd tt (m, m);
  // make tt the m-by-m-Zero-Matrix
  tt.setZero ();
  // set the diagonal entries of tt to the entries of P (cdf-values of x(k))
  tt.diagonal () = P;
  // add to z (the above defined matrix of bivariate probabilities) the matrix
  // tt of diagonal univariate probabilities
  z = z + tt;

  // Define Matrix z1 such that
  // z1(k,l) = Pr(X_k < x(k)) * Pr(X_l < x(l))

  // reserve memory for an m-by-m-Matrix of doubles called z1
  Eigen::MatrixXd z1 (m, m);
  // define z1(k,l) = P(k)*P(l)
  z1 = P * P.transpose ();

  // calculate Omega
  Eigen::MatrixXd Omega (m, m);

  Omega = z - z1;

  // go for derivatives

  // Define Matrix rho1 such that
  // rho1(k,l) = cor(X_k, X_l) = r(k,l) for k!=l
  // rho1(k,k) = 0

  // reserve memory for an m-by-m-Matrix of doubles called rho1
  Eigen::MatrixXd rho1 (m, m);
  // set rho1 equal to r
  rho1 = r;
  // set the diagonal entries of rho1 to 0
  rho1.diagonal () = zeros;

  // Define Matrix rho2 such that
  // rho2(k,l) = sqrt(1-r(k,l)^2) for k!=l
  // rho2(k,k) = 1

  // reserve memory for an m-by-m-Matrix of doubles called rho2
  Eigen::MatrixXd rho2 (m, m);
  rho2 =
    // calculate elementwise rho2 = sqrt(1-rho1^2) and interpret it as a
    // matrix
    (Mones.array () - rho1.array () * rho1.array ()).array ().sqrt ().matrix ();

  // Define Matrix muygx such that
  // muygx(k,l) = (x(l) - x(k) * r(k,l)) / (sqrt(1-r(k,l)^2)) for k!=l
  // muygx(k,k) = x(k)*(1-r(k,k)) = 0

  // reserve memory for an m-by-m-Matrix of doubles called muygx
  Eigen::MatrixXd muygx (m, m);
  muygx = (((ones * x.transpose ()
               - (r.array () * (x * ones.transpose ()).array ()).matrix ()).array ())
             // diagonal entries of muygx(i,i) = (1-r(i,i))*x(i)=0 and other elements
             // muygx(i,j) = (x(j)-x(i)*r(i,j))/(sqrt(1-r(i,j)^2))
             / rho2.array ()).matrix ();

             // Define Matrix nmuygx such that
             // nmuygx(k,l) = Pr(X_k < (x(l)-x(k)*r(k,l))/(sqrt(1-r(k,l)^2)) ) for k!=l
             // nmuygx(k,k) = Pr(X_k < x(k)*(1-r(k,k)) ) = 1/2

             // reserve memory for an m-by-m-Matrix of doubles called nmuygx
             Eigen::MatrixXd nmuygx (m, m);
             for (int k = 0; k < m; k++) {
               // for each row k of nmuygx
               for (int l = 0; l < m; l++) {
                 // for each element l of row k of nmuygx
                 nmuygx (k, l) = std_normal_cdf (muygx (k, l));
               }
             }

             // Define Matrix pdfmuygx such that
             // pdfmuygx(k,l) = \phi( (x(l) - x(k)*r(k,l))/(sqrt(1-r(k,l)^2)) ) for k!=l
             // pdfmuygx(k,k) = \phi( x(k)*(1-r(k,k)) ) = \phi( 0 )

             // reserve memory for an m-by-m-Matrix of doubles called pdfmuygx
             Eigen::MatrixXd pdfmuygx (m, m);
             for (int k = 0; k < m; k++) {
               // for each row k of pdfmuygx
               for (int l = 0; l < m; l++) {
                 // for each element l of row k of pdfmuygx
                 //
                 pdfmuygx (k, l) = std_normal_pdf (muygx (k, l));
               }
             }

             // Define Matrix pdiag such that
             // pdiag(k,k) = \phi(x(k))
             // pdiag(k,l) = 0 for k!=l

             // reserve memory for an m-by-m-Matrix of doubles called pdiag
             Eigen::MatrixXd pdiag (m, m);
             // set all elements of pdiag to 0
             pdiag.setZero ();
             // set the diagonal elements of pdiag to p = std_normal_pdf(x(k))
             pdiag.diagonal () = p;

             // Define Matrix dPhixy such that
             // dPhixy(k,l) = \phi(x(k)) * Pr(X < muygx(k,l))

             // reserve memory for an m-by-m-Matrix of doubles called dPhixy
             Eigen::MatrixXd dPhixy (m, m);
             // set dPhixy(i,j) = std_normal_pdf(x(i)) * std_normal_cdf(muygx(i,j))
             dPhixy = pdiag * nmuygx;

             // Define Matrix phiPhi such that
             // phiPhi(k,l) = \phi(x(k)) * Pr(X_l < x(l) )

             // reserve memory for an m-by-m-Matrix of doubles called phiPhi
             Eigen::MatrixXd phiPhi (m, m);
             // set phiPhi(i,j) = std_normal_pdf(x(i)) * std_normal_cdf(x(j))
             phiPhi = p * P.transpose ();

             // Define Matrix dCov such that
             // dCov(k,l) = dPhixy(k,l) - phiPhi(k,l)
             // = \phi(x(k)) * Pr(X < muygx(k,l)) - \phi(x(k)) * Pr(X_l < x(l))
             // = \phi(x(k)) * ( Pr(X < muygx(k,l)) - Pr(X_l < x(l)) ) for k!=l
             // dCov(k,k) = \phi(x(k)) * (1 - 2 * Pr(X_k < x(k)) ) updated later

             // reserve memory for an m-by-m-Matrix of doubles called dCov
             Eigen::MatrixXd dCov (m, m);
             // set dCov(i,j) = std_normal_pdf(x(i)) * ( std_normal_cdf(muygx(i,j)) -
             // std_normal_cdf(x(j)) )
             dCov = phiPhi - dPhixy;

             // Define Matrix dCovDiag such that
             // dCovDiag(k,k) = \phi(x(k)) * ( Pr(X < muygx(k,k)) - Pr(X_k < x(k)) )

             // reserve memory for an m-by-m-Matrix of doubles called dCovDiag
             Eigen::MatrixXd dCovDiag (m, m);
             // set all elements of dCovDiag to 0
             dCovDiag.setZero ();
             // set diagonal elements of dCovDiag to diagonal elements of dCov
             dCovDiag.diagonal () = dCov.diagonal ();

             // Define Matrix dCovDiag1 such that
             // dCovDiag1(k,k) = \phi(x(k)) * ( Pr(X < muygx(k,k)) - Pr(X_k < x(k)) )

             // reserve memory for an m-by-m-Matrix of doubles called dCovDiag1
             Eigen::MatrixXd dCovDiag1 (m, m);
             // set all elements of dCovDiag1 to 0
             dCovDiag1.setZero ();
             dCovDiag1.diagonal () = p.array ()
               // set diagonal elements of dCovDiag1 to dCovDiag1(i,i) =
               // std_normal_pdf(x(i)) * (1 - 2 * std_normal_cdf(x(i)))
               * (ones.array () - (ones.array () * 2) * P.array ());

             // Update Matrix dCov such that the diagonal is
             // dCov(k,k) = \phi(x(k)) * (1 - 2 * Pr(X_k < x(k)) )
             // set the diagonal elements of dCov to std_normal_pdf(x(i)) * (1 - 2 *
             // std_normal_cdf(x(i)))
             dCov = dCov - dCovDiag + dCovDiag1;
             // dCov(k,l) = d/(d x(k)) Cov(I_k,I_l)

             // derivatives with respect to correlations.
             // reserve memory for an m-by-m-Matrix of doubles called dPhi2
             Eigen::MatrixXd dPhi2 (m, m);
             // set dPhi2(i,j) = std_normal_pdf(x(i)) * std_normal_pdf(muygx(i,j)) /
             // sqrt(1-rho1(i,j)^2) with rho1(i,j) = r(i,j)
             dPhi2 = ((pdiag * pdfmuygx).array () / rho2.array ()).matrix ();
             // set the diagonal elements of dPhi2 to 0
             dPhi2.diagonal () = zeros;

             // calculate dOmega
             // number of derivatives
             int nd;
             nd = m + m * (m - 1) / 2;
             Rcpp::List dOmega (nd);
             // derivatives with respect to kk1
             // cycle over kk1 entries
             for (int k = 0; k < (m); k++) {
               Eigen::MatrixXd dom = Mzeros;
               // Bound guard was disabled in the original implementation:
               // if ((x(k) > -6) && (x(k) < 6))
               // adjust k-th row of gradient
               for (int j = 0; j < (m); j++) {
                 dom (k, j) = dPhixy (k, j) - p (k) * P (j);
                 dom (j, k) = dPhixy (k, j) - p (k) * P (j);
               }
               dom (k, k) = p (k) * (1 - 2 * P (k));

               dOmega[k] = dom;
             }

             // derivatives with respect to correlation entries.
             Eigen::VectorXd colI = trilCols (m);
             Eigen::VectorXd rowI = trilRows (m);
             // cycle over entries
             for (int j = 0; j < nd - m; j++) {
               int rc = rowI (j);
               int cc = colI (j);
               // corr r_{i,j} influences only element i,j.
               Eigen::MatrixXd dom = Mzeros;
               if ((r (rc, cc) > tol - 1) && (r (rc, cc) < 1 - tol)) {
                 dom (rc, cc) = dPhi2 (rc, cc);
                 dom (cc, rc) = dom (rc, cc);
               }
               dOmega[m + j] = dom;
             }

             // initialize d2Omega as a list of lists.
             Rcpp::List d2Omega (nd);

             // cases:
             Eigen::MatrixXd pcond (m, m);
             pcond = pdfmuygx.array () / rho2.array ();

             // first derive with respect to element in kk1.
             for (int k = 0; k < m; k++) {
               // list structure also in second direction
               Rcpp::List d2om (nd);

               Eigen::MatrixXd dPhi2r (m, m);
               dPhi2r = dPhi2.array () * r.array ();

               for (int j = 0; j < m; j++) {
                 // initialize to zero matrix.
                 Eigen::MatrixXd d2omh = Mzeros;
                 // second with respect to another element in kk1
                 if (j != k) {
                   // two different elements
                   d2omh (k, j) = dPhi2 (j, k) - p (j) * p (k);
                   d2omh (j, k) = d2omh (k, j);
                 }
                 if (j == k) {
                   // two times the same element
                   Eigen::MatrixXd dom_h = dOmega[j];
                   d2omh = dom_h * (-1) * x (j);
                   for (int c = 0; c < m; c++) {
                     d2omh (k, c) = d2omh (k, c) - dPhi2r (k, c);
                     d2omh (c, k) = d2omh (k, c);
                   }
                   d2omh (j, j) = dom_h (j, j) * (-1) * x (j) - 2 * p (j) * p (j);
                 }

                 // if ((x(k)>-6)&&(x(k)<6)&&(x(j)>-6)&&(x(j)<6)){
                 d2om[j] = d2omh;
                 // }
               }

               // second with respect to element in corr.
               // cycle over entries
               for (int j = 0; j < nd - m; j++) {
                 int rc = rowI (j);
                 int cc = colI (j);
                 // corr r_{i,j} influences only element i,j.
                 Eigen::MatrixXd d2omh = Mzeros;
                 if (k == rc) {
                   // non-zero derivative only, if kk1 elements corresponds to row index.
                   if ((r (rc, cc) > -0.99) && (r (rc, cc) < 0.99)) {
                     double rcur = r (rc, cc);
                     double onemr2 = 1 - rcur * rcur;
                     d2omh (rc, cc) = dPhi2 (rc, cc)
                       * (x (cc) * r (rc, cc) - x (rc)) / (onemr2);
                     d2omh (cc, rc) = d2omh (rc, cc);
                   }
                 }
                 if (k == cc) {
                   // non-zero derivative only, if kk1 elements corresponds to row index.
                   if ((r (rc, cc) > -0.99) && (r (rc, cc) < 0.99)) {
                     double rcur = r (rc, cc);
                     double onemr2 = 1 - rcur * rcur;
                     d2omh (rc, cc) = dPhi2 (rc, cc)
                       * (x (rc) * r (rc, cc) - x (cc)) / (onemr2);
                     d2omh (cc, rc) = d2omh (rc, cc);
                   }
                 }
                 d2om[m + j] = d2omh;
               }
               d2Omega[k] = d2om;

             }
             // then with respect to element in corr matrix.
             for (int j = 0; j < nd - m; j++) {
               // first with respect to an entry in the correlation matrix.

               // list structure also in second direction
               Rcpp::List d2om (nd);
               int rc = rowI (j);
               int cc = colI (j);
               // cycle over entries (derivs with respect to xk)
               for (int k = 0; k < m; k++) {

                 // corr r_{i,j} influences only element i,j.
                 Eigen::MatrixXd d2omh = Mzeros;
                 if (k == rc) {
                   // non-zero derivative only, if kk1 elements corresponds to row index.
                   if ((r (rc, cc) > -0.99) && (r (rc, cc) < 0.99)) {
                     double rcur = r (rc, cc);
                     double onemr2 = 1 - rcur * rcur;

                     d2omh (rc, cc) = dPhi2 (rc, cc)
                       * (x (cc) * r (rc, cc) - x (rc)) / (onemr2);
                     d2omh (cc, rc) = d2omh (rc, cc);
                   }
                 }
                 if (k == cc) {
                   // non-zero derivative only, if kk1 elements corresponds to row index.
                   if ((r (rc, cc) > -0.99) && (r (rc, cc) < 0.99)) {
                     double rcur = r (rc, cc);
                     double onemr2 = 1 - rcur * rcur;

                     d2omh (rc, cc) = dPhi2 (rc, cc)
                       * (x (rc) * r (rc, cc) - x (cc)) / (onemr2);
                     d2omh (cc, rc) = d2omh (rc, cc);
                   }
                 }
                 d2om[k] = d2omh;
               }

               d2Omega[m + j] = d2om;

               // finally two times entries in corr matrix.
               for (int jj = 0; jj < nd - m; jj++) {
                 Eigen::MatrixXd d2omh = Mzeros;
                 // this is zero unless we pick the same
                 if (j == jj) {
                   double rcur = r (rc, cc);
                   double onemr2 = 1 - rcur * rcur;

                   if ((r (rc, cc) > tol - 1) && (r (rc, cc) < 1 - tol)) {
                     double he = x (rc) * x (cc) * onemr2
                     - (x (rc) * x (rc) + x (cc) * x (cc)
                          - 2 * r (rc, cc) * x (rc) * x (cc)) * r (rc, cc);
                          d2omh (rc, cc) = dPhi2 (rc, cc)
                            * (he / (onemr2 * onemr2) + r (rc, cc) / onemr2);
                          d2omh (cc, rc) = d2omh (rc, cc);
                   }
                 }

                 d2om[m + jj] = d2omh;
               }
               d2Omega[m + j] = d2om;
             }

             // return the matrix Omega
             return d2Omega;
}

//' @describeIn SJ Hessian computation
//' @export
//'
// [[Rcpp::export]]

Eigen::MatrixXd SJ_hess (Eigen::VectorXd x, Eigen::MatrixXd r) {
  // save the number of rows of x in m
  int m = x.rows ();

  Eigen::VectorXd tol_which (m);
  tol_which.setZero ();

  // reserve memory for a vector of doubles of length m called zeroes
  Eigen::VectorXd zeros (m);
  // set all entries of zeroes to 0
  zeros.setZero ();

  // reserve memory for an m-by-m-Matrix of doubles called Mones
  Eigen::MatrixXd Mones (m, m);
  // set all entries of the m-by-m-Matrix Mones to 1
  Mones.setOnes ();

  // reserve memory for an m-by-m-Matrix of doubles called Mones
  Eigen::MatrixXd Mzeros (m, m);
  // set all entries of the m-by-m-Matrix Mones to 1
  Mzeros.setZero ();

  // reserve memory for vector of length m called ones with double entries
  Eigen::VectorXd ones (m);
  // set all entries of ones to 1
  ones.setOnes ();

  // Define Vector P such that P(k) = Pr(X_k < x_k)
  // reserve memory for the vector P of length m with double entries
  Eigen::VectorXd P (m);
  P.setZero ();

  for (int k = 0; k < m; k++) {
    // go through all entries k of P
    // compute the sandard normal cumulative distribution function at x=x(k)
    // and save it in P(k)
    P (k) = std_normal_cdf (x (k));
  }

  // calculate PDF for later usage.
  // reserve memory for the vector P of length m with double entries
  Eigen::VectorXd p (m);
  p.setZero ();
  for (int k = 0; k < m; k++) {
    // go through all entries k of P
    // compute the sandard normal cumulative distribution function at x=x(k)
    // and save it in P(k)
    p (k) = std_normal_pdf (x (k));
  }

  // Define Matrix z such that
  // z(k,l) = Pr(X_k < x_k, X_l < x_l)
  // z(k,k) = Pr(X_k < x_k)

  // define gradient and start calculations.
  int nd;
  nd = m + m * (m - 1) / 2;
  Eigen::MatrixXd Hess (nd, nd);
  Hess.setZero ();

  Eigen::MatrixXd H_d (nd, nd);
  H_d.setZero ();

  if (m > 1) {

    // calc initials
    // reserve memory for an m-by-m matrix of doubles and call it z
    Eigen::MatrixXd z (m, m);
    // set all entries of z to zero
    z.setZero ();
    for (int k = 0; k < (m); k++) {
      // go trough the rows k of z
      for (int l = (k + 1); l < (m); l++) {
        // go trough the elements l of the row k which are past the diagonal

        // save element r(k,l) in rrkl
        double rrkl = r (k, l);
        if (r (k, l) > 1 - tol) {
          // if r(k,l)>1
          // set rrkl = 0.99
          rrkl = 1 - tol;
        }
        if (r (k, l) < tol - 1) {
          // if r(k,l)<-1
          // set rrkl=-0.99
          rrkl = tol - 1;
        }

        // computes the probability for two normal variates X and Y whose
        // correlation is rrkl, that -x(k) <= X and -x(l) <= Y.
        z (k, l) = bivnor (-x (k), -x (l), rrkl);
        // make the matrix symmetric
        z (l, k) = z (k, l);
      }
    }

    // reserve memory for an m-by-m-Matrix of doubles called tt
    Eigen::MatrixXd tt (m, m);
    // make tt the m-by-m-Zero-Matrix
    tt.setZero ();
    // set the diagonal entries of tt to the entries of P (cdf-values of x(k))
    tt.diagonal () = P;
    // add to z (the above defined matrix of bivariate probabilities) the
    // matrix tt of diagonal univariate probabilities
    z = z + tt;

    Eigen::VectorXd ones_m_P (m);
    ones_m_P.setZero ();
    ones_m_P = ones - P;
    // calculate Omega
    Eigen::MatrixXd Om (m, m);

    // calc Omega
    Om.setZero ();
    Om = Omega (x, r);

    // calc dOmega
    Rcpp::List dOm (nd);
    dOm = dOmega (x, r);
    // calc d2Omega (a list of lists)
    Rcpp::List d2Om (nd);

    d2Om = d2Omega (x, r);
    // now Omega, dOmega and d2Omega is set -> run through the calculations
    for (int cur1 = 0; cur1 < m; cur1++) {
      // // cycle over first derivatives with respect to xk -> cur1
      Eigen::MatrixXd dom_cur1 (m, m);
      //
      dom_cur1.setZero ();
      dom_cur1 = dOm (cur1);

      Rcpp::List d2cur1 = d2Om (cur1);

      for (int cur2 = cur1; cur2 < m; cur2++) {
        // only upper triangular -> rest will be fixed later on.
        // cycle over second derivs with respect to xj -> cur2
        Eigen::MatrixXd dom_cur2 (m, m);
        //
        dom_cur2.setZero ();
        dom_cur2 = dOm (cur2);

        // second derivative of Omega.
        Eigen::MatrixXd d2cur1cur2 (m, m);
        d2cur1cur2.setZero ();

        d2cur1cur2 = d2cur1 (cur2);

        double r12 = r (0, 1);

        // deriv of first term.
        Eigen::MatrixXd d2Phi (3, 3);
        d2Phi.setZero ();
        d2Phi = Hess_cdf (x (0), x (1), r12);
        Eigen::VectorXd gr (3);
        gr = grad_cdf (x (0), x (1), r12);

        if ((cur1 == 0) && (cur2 < 2) && (z (0, 1) > tol)) {
          Hess (cur1, cur2) = d2Phi (cur1, cur2) / z (0, 1)
          - gr (cur1) * gr (cur2) / (z (0, 1) * z (0, 1));
        }
        if ((cur1 == 1) && (cur2 == 1) && (z (0, 1) > tol)) {
          Hess (cur1, cur2) = d2Phi (cur1, cur2) / z (0, 1)
          - gr (cur1) * gr (cur2) / (z (0, 1) * z (0, 1));
        }

        // then derivative of other terms.

        for (int k = 2; k < m; k++) {
          // k runs to calculate condk.
          Eigen::VectorXd omega12 (k);
          omega12.setZero ();
          omega12 = Om.block (0, k, k, 1);
          Eigen::MatrixXd omega11 (k, k);
          omega11.setZero ();
          omega11 = Om.block (0, 0, k, k);
          Eigen::MatrixXd iom (k, k);
          iom.setZero ();
          iom = omega11.inverse ();

          double condk;
          condk = 0;

          Eigen::VectorXd he_v (k);
          he_v.setZero ();
          he_v = iom * ones_m_P.block (0, 0, k, 1);

          Eigen::VectorXd oio (k);
          oio.setZero ();
          oio = iom * omega12;

          condk = P (k) + omega12.dot (he_v);
          if (condk > tol) {
            // check if condk = std_normal_cdf(x(k)) + tempcon is greater than
            // the
            // tollarance
            // no truncation -> contribution to derivative // if so, lcondk =
            // std::log (condk)
            double d2condk = 0;

            // get the various derivatives.
            Eigen::VectorXd dcb1_a (k);
            dcb1_a.setZero ();
            dcb1_a = dom_cur1.block (0, k, k, 1);

            Eigen::VectorXd dcb1_b (k);
            dcb1_b.setZero ();
            dcb1_b = dom_cur2.block (0, k, k, 1);

            Eigen::MatrixXd dcb2_a (k, k);
            dcb2_a.setZero ();
            dcb2_a = dom_cur1.block (0, 0, k, k);

            Eigen::MatrixXd dcb2_b (k, k);
            dcb2_b.setZero ();
            dcb2_b = dom_cur2.block (0, 0, k, k);

            Eigen::VectorXd d2P (k);
            d2P.setZero ();
            d2P = zeros.block (0, 0, k, 1);
            if ((cur1 < k) && (cur1 == cur2)) {
              d2P (cur1) = p (cur1) * (x (cur1));
            }

            Eigen::VectorXd dP_a (k);
            dP_a.setZero ();
            if (cur1 < k) {
              dP_a (cur1) = p (cur1) * (-1);
            }

            Eigen::VectorXd dP_b (k);
            dP_b.setZero ();
            if (cur2 < k) {
              dP_b (cur2) = p (cur2) * (-1);
            }

            // first derivatives with respect to the two entries
            double dcondk_a = 0;
            if (k == cur1) {
              dcondk_a = std_normal_pdf (x (k));
            }
            dcondk_a = dcondk_a + dcb1_a.transpose () * he_v;
            dcondk_a = dcondk_a
            - (oio.transpose () * (dcb2_a * he_v));
            if (k > cur1) {
              dcondk_a = dcondk_a - oio (cur1) * p (cur1);
            }

            double dcondk_b = 0;
            if (k == cur2) {
              dcondk_b = std_normal_pdf (x (k));
            }
            dcondk_b = dcondk_b + dcb1_b.transpose () * he_v;
            dcondk_b = dcondk_b
            - (oio.transpose () * (dcb2_b * he_v));
            if (k > cur2) {
              dcondk_b = dcondk_b - oio (cur2) * p (cur2);
            }

            // now collect all the different terms.
            // term Phi(x_k)
            if ((cur1 == k) && (cur1 == cur2)) {
              d2condk = -p (k) * x (k);
            }

            // // term alpha derived w.r.t. cur1.
            Eigen::VectorXd dalpha_a (k);
            dalpha_a.setZero ();
            dalpha_a = iom * dcb1_a;

            Eigen::VectorXd d2c1c2 (k);
            d2c1c2.setZero ();
            d2c1c2 = d2cur1cur2.block (0, k, k, 1);

            d2condk = d2condk + d2c1c2.transpose () * he_v
            - dalpha_a.transpose () * dcb2_b * he_v
            + dalpha_a.transpose () * dP_b;
            // Om^(-1) derived with r.t.alpha.
            Eigen::MatrixXd d2cb2_ab (k, k);
            d2cb2_ab = iom * dcb2_b * iom * dcb2_a * iom
            - iom * d2cur1cur2.block (0, 0, k, k) * iom
            + iom * dcb2_a * iom * dcb2_b * iom;

            Eigen::VectorXd alpha_dOminv (k);
            alpha_dOminv.setZero ();
            alpha_dOminv = -iom * dcb2_a * iom * omega12;

            Eigen::VectorXd ompk (k);
            ompk = ones_m_P.block (0, 0, k, 1);

            d2condk = d2condk
            - dcb1_b.transpose () * iom * dcb2_a * he_v
            + omega12.transpose () * d2cb2_ab * ompk
            + alpha_dOminv.transpose () * dP_b;

            // finally ones_minus_P derived with respect to cur1.
            d2condk = d2condk + dcb1_b.transpose () * iom * dP_a
            - omega12.transpose () * iom * dcb2_b * iom * dP_a
            + oio.transpose () * d2P;

            Hess (cur1, cur2) = d2condk / condk
            - dcondk_a * dcondk_b / (condk * condk)
              + Hess (cur1, cur2);

          }

        }

      }

      for (int curr = 0; curr < nd - m; curr++) {
        // derivative first with respect to xk, then with respect to rij.
        double r12;
        r12 = r (0, 1);

        Eigen::MatrixXd dom_curr (m, m);
        //
        dom_curr.setZero ();
        dom_curr = dOm (m + curr);

        // second derivative of Omega.
        Eigen::MatrixXd d2cur1cur2 (m, m);
        d2cur1cur2.setZero ();
        d2cur1cur2 = d2cur1 (m + curr);

        // deriv of first term.
        Eigen::MatrixXd d2Phi (3, 3);
        d2Phi.setZero ();
        d2Phi = Hess_cdf (x (0), x (1), r12);
        Eigen::VectorXd gr (3);
        gr = grad_cdf (x (0), x (1), r12);

        if ((cur1 < 2) && (curr == 0) && (z (0, 1) > tol)) {
          Hess (cur1, m) = d2Phi (cur1, 2) / z (0, 1)
          - gr (cur1) * gr (2) / (z (0, 1) * z (0, 1));
        }

        // for the rest only the cross derivs count
        for (int k = 2; k < m; k++) {
          // k runs to calculate condk.
          Eigen::VectorXd omega12 (k);
          omega12.setZero ();
          omega12 = Om.block (0, k, k, 1);
          Eigen::MatrixXd omega11 (k, k);
          omega11.setZero ();
          omega11 = Om.block (0, 0, k, k);
          Eigen::MatrixXd iom (k, k);
          iom.setZero ();
          iom = omega11.inverse ();

          double condk;
          condk = 0;

          Eigen::VectorXd he_v (k);
          he_v.setZero ();
          he_v = iom * ones_m_P.block (0, 0, k, 1);

          Eigen::VectorXd oio (k);
          oio.setZero ();
          oio = iom * omega12;

          condk = P (k) + omega12.transpose () * he_v;
          if (condk > tol) {
            // check if condk = std_normal_cdf(x(k)) + tempcon is greater than
            // the
            // tolerance
            // no truncation -> contribution to derivative // if so, lcondk =
            // std::log (condk)
            double d2condk = 0;

            // get the various derivatives.
            Eigen::VectorXd dcb1_a (k);
            dcb1_a.setZero ();
            dcb1_a = dom_cur1.block (0, k, k, 1);

            Eigen::VectorXd dcb1_b (k);
            dcb1_b.setZero ();
            dcb1_b = dom_curr.block (0, k, k, 1);

            Eigen::MatrixXd dcb2_a (k, k);
            dcb2_a.setZero ();
            dcb2_a = dom_cur1.block (0, 0, k, k);

            Eigen::MatrixXd dcb2_b (k, k);
            dcb2_b.setZero ();
            dcb2_b = dom_curr.block (0, 0, k, k);

            Eigen::VectorXd dP_a (k);
            dP_a.setZero ();
            if (cur1 < k) {
              dP_a (cur1) = p (cur1) * (-1);
            }

            // first derivatives with respect to the two entries
            double dcondk_a = 0;
            if (k == cur1) {
              dcondk_a = std_normal_pdf (x (k));
            }
            dcondk_a = dcondk_a + dcb1_a.transpose () * he_v;
            dcondk_a = dcondk_a
            - (oio.transpose () * (dcb2_a * he_v));
            if (k > cur1) {
              dcondk_a = dcondk_a - oio (cur1) * p (cur1);
            }

            double dcondk_b = 0;
            dcondk_b = dcb1_b.transpose () * he_v;
            dcondk_b = dcondk_b
            - (oio.transpose () * (dcb2_b * he_v));

            // now collect all the different terms.
            // // term alpha derived w.r.t. cur1.
            Eigen::VectorXd dalpha_a (k);
            dalpha_a.setZero ();
            dalpha_a = iom * dcb1_a;

            Eigen::VectorXd d2c1c2 (k);
            d2c1c2.setZero ();
            d2c1c2 = d2cur1cur2.block (0, k, k, 1);

            d2condk = d2condk + d2c1c2.transpose () * he_v
            - dalpha_a.transpose () * dcb2_b * he_v;
            // Om^(-1) derived with r.t.alpha.
            Eigen::MatrixXd d2cb2_ab (k, k);
            d2cb2_ab = iom * dcb2_b * iom * dcb2_a * iom
            - iom * d2cur1cur2.block (0, 0, k, k) * iom
            + iom * dcb2_a * iom * dcb2_b * iom;

            // VectorXd alpha_dOminv(k);
            // alpha_dOminv.setZero();
            // alpha_dOminv = -iom * dcb2_a * iom * omega12;

            Eigen::VectorXd ompk (k);
            ompk = ones_m_P.block (0, 0, k, 1);

            d2condk = d2condk
            - dcb1_b.transpose () * iom * dcb2_a * he_v
            + omega12.transpose () * d2cb2_ab * ompk;

            // finally ones_minus_P derived with respect to cur1.
            d2condk = d2condk + dcb1_b.transpose () * iom * dP_a
            - omega12.transpose () * iom * dcb2_b * iom * dP_a;

            Hess (cur1, m + curr) = d2condk / condk
            - (dcondk_a * dcondk_b) / (condk * condk)
              + Hess (cur1, m + curr);

          }
        }
      }

    }

    for (int curr = 0; curr < nd - m; curr++) {
      // derivatives with respect to pairs of correlation entries.
      // Again only upper triangular part.
      // // cycle over first derivatives with respect to xk -> cur1
      Eigen::MatrixXd dom_curr (m, m);
      //
      dom_curr.setZero ();
      dom_curr = dOm (curr + m);

      Rcpp::List d2curr (nd);
      d2curr = d2Om (curr + m);

      for (int curd = curr; curd < nd - m; curd++) {

        if ((curr == 0) && (curd == 0) && (z (0, 1) > tol)) {
          // only the first corr influences the first term.
          double r12;
          r12 = r (0, 1);

          // deriv of first term.
          Eigen::MatrixXd d2Phi (3, 3);
          d2Phi.setZero ();
          d2Phi = Hess_cdf (x (0), x (1), r12);
          Eigen::VectorXd gr (3);
          gr = grad_cdf (x (0), x (1), r12);

          // cross derivatives are always zero for first term.
          Hess (m, m) = d2Phi (2, 2) / z (0, 1)
            - gr (2) * gr (2) / (z (0, 1) * z (0, 1));
        }
        // cycle over second derivs with respect to second corr. -> curd
        Eigen::MatrixXd dom_curd (m, m);
        //
        dom_curd.setZero ();
        dom_curd = dOm (curd + m);

        // second derivative of Omega.
        Eigen::MatrixXd d2currcurd (m, m);
        d2currcurd.setZero ();
        d2currcurd = d2curr (curd + m);

        // then derivative of other terms.
        for (int k = 2; k < m; k++) {

          // Rcout << "k:" << k <<std::endl;

          // k runs to calculate condk.
          Eigen::VectorXd omega12 (k);
          omega12.setZero ();
          omega12 = Om.block (0, k, k, 1);
          Eigen::MatrixXd omega11 (k, k);
          omega11.setZero ();
          omega11 = Om.block (0, 0, k, k);
          Eigen::MatrixXd iom (k, k);
          iom.setZero ();
          iom = omega11.inverse ();

          double condk;
          condk = 0;

          Eigen::VectorXd he_v (k);
          he_v.setZero ();
          he_v = iom * ones_m_P.block (0, 0, k, 1);

          Eigen::VectorXd oio (k);
          oio.setZero ();
          oio = iom * omega12;

          condk = P (k) + omega12.transpose () * he_v;

          if (condk > tol) {
            // check if condk = std_normal_cdf(x(k)) + tempcon is greater than
            // the
            // tolerance
            // no truncation -> contribution to derivative // if so, lcondk =
            // std::log (condk)
            double d2condk = 0;

            // get the various derivatives.
            Eigen::VectorXd dcb1_a (k);
            dcb1_a.setZero ();
            dcb1_a = dom_curr.block (0, k, k, 1);

            Eigen::VectorXd dcb1_b (k);
            dcb1_b.setZero ();
            dcb1_b = dom_curd.block (0, k, k, 1);

            Eigen::MatrixXd dcb2_a (k, k);
            dcb2_a.setZero ();
            dcb2_a = dom_curr.block (0, 0, k, k);

            Eigen::MatrixXd dcb2_b (k, k);
            dcb2_b.setZero ();
            dcb2_b = dom_curd.block (0, 0, k, k);

            // first derivatives with respect to the two entries
            double dcondk_a = 0;

            dcondk_a = dcondk_a + dcb1_a.transpose () * he_v;
            dcondk_a = dcondk_a
            - (oio.transpose () * (dcb2_a * he_v));

            double dcondk_b = 0;
            dcondk_b = dcondk_b + dcb1_b.transpose () * he_v;
            dcondk_b = dcondk_b
            - (oio.transpose () * (dcb2_b * he_v));

            // now collect all the different terms.

            // term alpha derived w.r.t. cur1.
            Eigen::VectorXd dalpha_a (k);
            dalpha_a.setZero ();
            dalpha_a = iom * dcb1_a;

            Eigen::VectorXd d2c1c2 (k);
            d2c1c2.setZero ();
            d2c1c2 = d2currcurd.block (0, k, k, 1);

            d2condk = d2condk + d2c1c2.transpose () * he_v
            - dalpha_a.transpose () * dcb2_b * he_v;
            // Om^(-1) derived with r.t.alpha.
            Eigen::MatrixXd d2cb2_ab (k, k);
            d2cb2_ab = iom * dcb2_b * iom * dcb2_a * iom
            - iom * d2currcurd.block (0, 0, k, k) * iom
            + iom * dcb2_a * iom * dcb2_b * iom;

            Eigen::VectorXd alpha_dOminv (k);
            alpha_dOminv.setZero ();
            alpha_dOminv = -iom * dcb2_a * iom * omega12;

            Eigen::VectorXd ompk (k);
            ompk = ones_m_P.block (0, 0, k, 1);

            d2condk = d2condk
            - dcb1_b.transpose () * iom * dcb2_a * he_v
            + omega12.transpose () * d2cb2_ab * ompk;

            Hess (curr + m, curd + m) = d2condk / condk
            - dcondk_a * dcondk_b / (condk * condk)
              + Hess (curr + m, curd + m);

          }
        }
      }

    }

    // now symmetrize the Hessian
    H_d = Hess + Hess.transpose ();
    H_d.diagonal () = Hess.diagonal ();
  } else {
    double pP = p (0) / P (0);
    H_d (0, 0) = -pP * x (0) - pP * pP;
  }

  return H_d;
}
