#include <RcppEigen.h>

#include <cmath>

#include "SJ.h"
#include "macml_grad_Hessian.h"
#include "lin_alg.h"

//' Calculate Approximated Composite Marginal Log-Likelihood
//'
//' @param theta \[`numeric()`\]\cr
//' Parameter vector.
//'
//' @param design_diff \[`list`\]\cr
//' Differenced design list as returned by
//' `choice_data$build_design(..., difference = TRUE)`.
//'
//' @param mod_list \[`list`\]\cr
//' Details of the model specification as returned by `mod$mod_list`.
//'
//' @param cml_pair_type \[`integer(1)`\]\cr
//' Pair construction type:
//' - `0` for full pairwise,
//' - `1` for adjacent pairwise,
//' - `2` for cyclical adjacent pairwise.
//'
//' @param calc_hess \[`logical(1)`\]\cr
//' Set to `TRUE` to calculate the Hessian and `FALSE` otherwise.
//'
//' @return
//' A numeric vector containing the negative log-likelihood with gradient and
//' Hessian attributes.
//'
//' @keywords internal
//'
//' @noRd

Rcpp::NumericVector ll_macml (
  Eigen::VectorXd theta,
  Rcpp::List design_diff,
  Rcpp::List mod_list,
  int cml_pair_type,
  bool calc_hess
) {

  // initialize output
  Rcpp::NumericVector out;

  int N = design_diff.length();
  int length_theta_b = mod_list["length_theta_b"];
  int length_theta_O = mod_list["length_theta_O"];
  int length_theta_L = mod_list["length_theta_L"];
  Rcpp::IntegerVector b_count = mod_list["b_count"];
  int length_RE = b_count[1];
  int alt = mod_list["alt"];
  int length_theta = length_theta_b + length_theta_O + length_theta_L;

  int hess = 0;
  if (calc_hess){
    hess = 1;
  }

  double ll = 0;

  Rcpp::List syst(7);
  syst = syst_from_mod (theta, mod_list);

  // Recover the vector b
  Eigen::VectorXd b = syst["b"];
  int lthbb = b.rows ();

  // Recover the Sigma matrix
  Eigen::MatrixXd Sigma_chol_vech = syst["Sigma_chol_vech"];
  Eigen::MatrixXd Sigma = syst["Sigma"];
  int lth_Sigma_chol_vech = Sigma_chol_vech.rows();

  // Recover the matrix Omega
  Eigen::MatrixXd Omega_chol_vech = syst["Omega_chol_vech"];
  Eigen::MatrixXd Omega = syst["Omega"];
  int lth_Omega_chol_vech = Omega_chol_vech.rows();

  // ------------------------------------------------------------//
  // SETUP ONE-OFF COMPUTATIONS FOR GRADIENT AND HESSIAN
  // ------------------------------------------------------------//

  // Only for pairwise with equal number of alternatives for all deciders.
  Eigen::MatrixXd D_JOm = duplmat (length_RE);

  // Compute matrix dimensions from vech lengths.
  int m_l = -0.5 + std::sqrt (0.25 + 2 * lth_Sigma_chol_vech);
  int m_o = -0.5 + std::sqrt (0.25 + 2 * lth_Omega_chol_vech);

  Eigen::MatrixXd L_JcholL = elimmat (m_l);
  Eigen::MatrixXd D_JcholL = duplmat (m_l);
  Eigen::MatrixXd K_JcholL = commmat (m_l, m_l);

  Eigen::MatrixXd L_JcholO = elimmat (m_o);
  Eigen::MatrixXd D_JcholO = duplmat (m_o);
  Eigen::MatrixXd K_JcholO = commmat (m_o, m_o);

  // This is used later to compute the Hessian of the Cholesky, which does not
  // depend on the specific value of the parameters.
  Eigen::VectorXd HoO = Omega_chol_vech.array () + 1;
  Eigen::VectorXd HlL = Sigma_chol_vech.array () + 1;

  Eigen::MatrixXd HJ2 = J_bvOvSig_bLOLSig (
    HoO, HlL, lthbb, L_JcholO, K_JcholO, D_JcholO, L_JcholL,
    K_JcholL, D_JcholL
  );

  // ------------------------------------------------------------//
  // END of setup
  // ------------------------------------------------------------//

  // ------------------------------------------------------------//
  // INITIALIZE GRADIENT OUTPUT
  // ------------------------------------------------------------//

  Eigen::MatrixXd Hess (length_theta, length_theta);
  Hess.setZero ();
  Eigen::MatrixXd Hess_approx (length_theta, length_theta);
  Hess_approx.setZero ();
  Eigen::MatrixXd Hess_approx2 (length_theta, length_theta);
  Hess_approx2.setZero ();
  Eigen::VectorXd grad (length_theta);
  grad.setZero ();

  double ll_n = 0;

  for (int i_n = 0; i_n < N; i_n++) { //loop over the sample

    // Initialize log-likelihood, gradient, and Hessian contributions for the
    // current decision maker.
    ll_n = 0.0;

    Eigen::VectorXd grad_n (length_theta);
    grad_n.setZero ();

    Eigen::MatrixXd Hess_n (length_theta, length_theta);
    Hess_n.setZero ();

    Eigen::MatrixXd Hess_approx_n (length_theta, length_theta);
    Hess_approx_n.setZero ();

    // Unpack the design data for the current decision maker.
    Rcpp::List data_n = design_diff[i_n];
    Rcpp::List X_n = data_n["X"];
    Eigen::VectorXd y_n = data_n["y"];

    // ------------------------------------------------------------//
    // CML SETUP
    // ------------------------------------------------------------//

    // ------------------------------------------------------------//
    // START: Build the pairs structure
    // ------------------------------------------------------------//

    // determine the panel length
    int Tp_n = y_n.rows ();


    // Define a matrix that allows to specify the pairs (possible extension to
    // higher order CMLs)
    Eigen::MatrixXd pairs_n;
    int lth_pairs_n = 0;
    if (Tp_n>1){
      pairs_n = cal_pairs_matrix (Tp_n, cml_pair_type);
      lth_pairs_n = pairs_n.rows();
    }

    // ------------------------------------------------------------//
    // END: Build the pairs structure
    // ------------------------------------------------------------//

    if (lth_pairs_n == 0) { // Tp_n = 1, only one decision taken

      // Current decision maker has only one choice occasion.
      int alt2 = (alt - 1);

      // setup X, M, M_Sig.
      Eigen::MatrixXd X (alt2, lthbb);
      X.setZero ();

      Eigen::MatrixXd M_Sig (alt2, alt2);
      M_Sig.setZero ();

      // ------------------------------------------------------------//
      // Data Prep and Lin Utilities
      // ------------------------------------------------------------//
      //

      // Extract the chosen alternative for the current choice occasion.
      int y1;
      y1 = y_n[0];

      // Extract the covariance matrix for the current choice occasion.
      X = X_n[0];

      Eigen::MatrixXd M1 = makeM ((alt), (y1 - 1));
      Eigen::MatrixXd M1t = M1.transpose ();
      Eigen::MatrixXd Sigma_norm_1 = M1 * Sigma * M1t;

      // ------------------------------------------------------------//
      // Covariance Matrix
      // ------------------------------------------------------------//

      M_Sig = Sigma_norm_1;

      // find the NA elements

      Eigen::MatrixXd X_red = X;
      Eigen::MatrixXd M_Sig_red = M_Sig;

      Rcpp::IntegerVector ind = NA_ind (X);

      int cur = ind.size ();
      if ((cur == 1) && (ind[0] == -1)) {
        continue;
      }

      for (int ij = 0; ij < cur; ij++) {
        X_red.row (ij) = X.row (ind[ij]);
        M_Sig_red.row (ij) = M_Sig.row (ind[ij]);
      }

      for (int ij = 0; ij < cur; ij++) {
        M_Sig_red.col (ij) = M_Sig_red.col (ind[ij]);
      }

      // unstandardized lin-utilties (upper limits of the MVNCDF)
      Eigen::VectorXd x_upper = (-1)
      * (X_red.block (0, 0, cur, X_red.cols ()) * b);
      //
      // compute covariance matrix Omega of the random effects for this pair of
      // observations using the observed effects
      Eigen::MatrixXd X_RE = X_red.block (0, 0, cur, length_RE);
      Eigen::MatrixXd X_RE_t = X_RE.transpose ();
      //


      // // Compute combined covariance matrix
      // // Changed Sig to Lambda for combined Covariance
      Eigen::MatrixXd Lambda = M_Sig_red.block (0, 0, cur, cur);
      if (length_RE != 0) {
        Lambda = Lambda + X_RE * Omega * X_RE_t;
      }

      Eigen::MatrixXd lambda = Lambda.diagonal ().array ().sqrt ();
      Eigen::MatrixXd lambda_t = lambda.transpose ();

      // ------------------------------------------------------------//
      //  Normalize upper limits and covariance matrix (to correlation matrix)
      // ------------------------------------------------------------//
      //
      Eigen::VectorXd x_norm = x_upper.array () / lambda.array ();
      Eigen::MatrixXd lambda_li_lj = lambda * lambda_t;
      Eigen::MatrixXd Lambda_cor =
      (Lambda.array () / lambda_li_lj.array ()).matrix ();

      // ------------------------------------------------------------//
      // Likelihood and Gradient
      // ------------------------------------------------------------//

      //double lp;
      //     Eigen::VectorXd grad_output;
      //     grad_output.setZero ();
      //     int lth_output;
      //
      //     grad_output = SJ_grad (x_norm, Lambda_cor);
      //
      //     lth_output = grad_output.size ();
      //     ll_n -= grad_output[lth_output - 1];

      // ------------------------------------------------------------//
      // Gradient
      // ------------------------------------------------------------//
      //
      // Eigen::VectorXd grad_MVNCDF_approx;
      // grad_MVNCDF_approx = grad_output.head (lth_output - 1);
      //
      // Eigen::VectorXd grad_output;
      // grad_output.setZero ();
      // int lth_output;

      Rcpp::List outn = cal_grad_Hess(
        x_norm, Lambda_cor, Lambda, mod_list, X_red, cur, alt, ind, y_n, D_JOm,
        syst, L_JcholO, K_JcholO, D_JcholO, L_JcholL, K_JcholL, D_JcholL, HJ2,
        hess
      );
      double lln = outn["ll"];
      ll_n -= lln;

      // ------------------------------------------------------------//
      // Gradient
      // ------------------------------------------------------------//
      //

      // initialize gradient contribution of the current pair
      Eigen::VectorXd grad_n_pair (length_theta);
      grad_n_pair = outn["gradient"];

      // add gradient of the current pair to the gradient of the decision maker
      grad_n -= grad_n_pair;  // added the weight

      // ------------------------------------------------------------//
      // Hessian
      // ------------------------------------------------------------//

      Hess_approx2 += grad_n_pair * grad_n_pair.transpose ();

      Eigen::MatrixXd Hess_n_pair (length_theta, length_theta);
      Hess_n_pair.setZero();

      if (hess == 1) {

        // update Hessian
        Hess_n_pair = outn["hessian"];
        Hess_n += Hess_n_pair;

      }

    }


    // normalisation factor
    double norm_factor = 1.0;
    double Tpn = (double)Tp_n;
    if (Tp_n > 1) {

      if (cml_pair_type == 0) {
        norm_factor = 1.0 / (Tpn - 1.0);
      }
      if (cml_pair_type == 1) {
        norm_factor = Tp_n == 2 ? 1.0 : 0.5; // Tp_n == 2 yields one pair
      }
      if (cml_pair_type == 2) {
        norm_factor = Tpn / (2.0 * (Tpn - 1.0));
      }
    }


    // NOW LOOP OVER THE PAIRS MATRIX to compute the CML
    for (int i_pair = 0; i_pair < lth_pairs_n; i_pair++) {

      // only calculate the pair if the weight is !=0
      if (pairs_n (i_pair, 2) != 0) {
        // ------------------------------------------------------------//
        // Data Prep and Lin Utilities
        // ------------------------------------------------------------//

        // extract the chosen alternatives for the current pair of choice
        // occasions
        int y1, y2;
        y1 = y_n[pairs_n (i_pair, 0)];
        y2 = y_n[pairs_n (i_pair, 1)];

        Eigen::VectorXd y12 (2);
        y12 (0) = y1;
        y12 (1) = y2;

        // Extract covariance matrices for the current pair.
        Eigen::MatrixXd X1, X2;
        X1 = X_n[pairs_n (i_pair, 0)];
        X2 = X_n[pairs_n (i_pair, 1)];

        // Joint matrix for pairs
        // This keeps flexibility for the number of alternatives in each period.
        int alt1 = X1.rows ();
        // do not expect this flexiblity beyond this point
        // TODO: if we want to allow different numbers of choice alternatives
        // for different time points, we would also have to know which
        // choice alternatives have to be excluded to cut the Lambda matrix and
        // the theta vector accordingly
        int alt2 = X2.rows ();
        Eigen::MatrixXd X_pair_norm ((alt1 + alt2), lthbb);
        X_pair_norm.setZero ();
        X_pair_norm << X1, X2;

        // ------------------------------------------------------------//
        // Covariance Matrix
        // ------------------------------------------------------------//

        // Use M1 and M2 to obtain Sigma after normalizing for the chosen
        // alternative.
        Eigen::MatrixXd M1 = makeM ((alt1 + 1), (y1 - 1));
        Eigen::MatrixXd M2 = makeM ((alt2 + 1), (y2 - 1));
        Eigen::MatrixXd M1t = M1.transpose ();
        Eigen::MatrixXd M2t = M2.transpose ();

        //       << "Sig:" << std::endl;
        Eigen::MatrixXd Sigma_norm_1 = M1 * Sigma * M1t;
        Eigen::MatrixXd Sigma_norm_2 = M2 * Sigma * M2t;

        // compute joint Covariance matrix for the current pair
        Eigen::MatrixXd Sigma_norm_pair ((alt1 + alt2), (alt1 + alt2));
        Sigma_norm_pair.setZero ();
        Sigma_norm_pair.block (0, 0, alt1, alt1) = Sigma_norm_1;
        Sigma_norm_pair.block (alt1, alt1, alt2, alt2) = Sigma_norm_2;

        // find the NA elements

        Eigen::MatrixXd X_red = X_pair_norm;
        Eigen::MatrixXd M_Sig_red = Sigma_norm_pair;

        Rcpp::IntegerVector ind = NA_ind (X_red);

        int cur = ind.size ();
        if ((cur == 1) && (ind[0] == -1)) {
          continue;
        }

        for (int ij = 0; ij < cur; ij++) {
          X_red.row (ij) = X_pair_norm.row (ind[ij]);
          M_Sig_red.row (ij) = Sigma_norm_pair.row (ind[ij]);
        }

        for (int ij = 0; ij < cur; ij++) {
          M_Sig_red.col (ij) = M_Sig_red.col (ind[ij]);
        }

        // // unstandardized lin-utilties (upper limits of the MVNCDF)

        Eigen::VectorXd x_upper = (-1)
        * (X_red.block (0, 0, cur, X_pair_norm.cols ()) * b);

        // Compute Omega of the random effects for this pair of observations
        // using the observed effects.
        Eigen::MatrixXd X_RE = X_red.block (0, X_red.cols() - length_RE, cur, length_RE);
        Eigen::MatrixXd X_RE_t = X_RE.transpose ();


        // Compute combined covariance matrix
        // Changed Sig to Lambda for combined Covariance
        Eigen::MatrixXd Lambda = M_Sig_red.block (0, 0, cur, cur);
        if (length_RE != 0) {
          Lambda = Lambda + X_RE * Omega * X_RE_t;
        }

        Eigen::MatrixXd lambda = Lambda.diagonal ().array ().sqrt ();
        Eigen::MatrixXd lambda_t = lambda.transpose ();

        // ------------------------------------------------------------//
        // Normalize upper limits and covariance matrix (to correlation matrix)
        // ------------------------------------------------------------//

        Eigen::VectorXd x_norm = x_upper.array () / lambda.array ();
        Eigen::MatrixXd lambda_li_lj = lambda * lambda_t;
        Eigen::MatrixXd Lambda_cor = (Lambda.array ()
        / lambda_li_lj.array ()).matrix ();

        // ------------------------------------------------------------//
        // Likelihood and Gradient
        // ------------------------------------------------------------//

        Rcpp::List outn = cal_grad_Hess(
          x_norm, Lambda_cor, Lambda, mod_list, X_red, cur, alt, ind, y12,
          D_JOm, syst, L_JcholO, K_JcholO, D_JcholO, L_JcholL, K_JcholL,
          D_JcholL, HJ2, hess
        );
        double lln = outn["ll"];
        ll_n -= (lln *  pairs_n (i_pair, 2))*norm_factor;

        // ------------------------------------------------------------//
        // Gradient
        // ------------------------------------------------------------//
        //

        // initialize gradient contribution of the current pair
        Eigen::VectorXd grad_n_pair (length_theta);
        grad_n_pair = outn["gradient"];

        // Add gradient of the current pair to the current decision maker.
        grad_n -= grad_n_pair * pairs_n (i_pair, 2) * norm_factor;  // added the weight

        // ------------------------------------------------------------//
        // Hessian
        // ------------------------------------------------------------//

        Hess_approx2 +=
          grad_n_pair * grad_n_pair.transpose () * pairs_n (i_pair, 2) * norm_factor;

        Eigen::MatrixXd Hess_n_pair (
          length_theta_b + length_theta_O + length_theta_L,
          length_theta_b + length_theta_O + length_theta_L
        );
        Hess_n_pair.setZero ();

        if (hess == 1) {

          // update Hessian
          Hess_n_pair = outn["hessian"];
          Hess_n += Hess_n_pair * pairs_n (i_pair, 2)*norm_factor;

        } //end if hessian

      } // end if-weight!=0

    } //end pairs-loop

    grad += grad_n; //it is already negative in the pairs loop
    ll += ll_n; //it is already negative in the pairs loop

    Hess_approx_n = grad_n * grad_n.transpose ();
    // This is positive throughout because it is only used this way.
    Hess_approx += Hess_approx_n;

    if (hess == 1) {
      Hess += Hess_n;
    }

  } //end N-loop

  // ------------------------------------------------------------//
  // Format the output
  // ------------------------------------------------------------//
  out = ll;
  out.attr ("gradient") = (grad.array ());
  out.attr ("hessian_approx") = (Hess_approx2.array ());
  out.attr ("Vhat") = (Hess_approx.array ());

  if (hess == 1) {
    out.attr ("hessian") = (Hess.array ());
  }

  return out;

}

//' Calculate Full Probit Log-Likelihood
//'
//' @param theta \[`numeric()`\]\cr
//' Parameter vector.
//'
//' @param design_diff \[`list`\]\cr
//' Differenced design list as returned by
//' `choice_data$build_design(..., difference = TRUE)`.
//'
//' @param mod_list \[`list`\]\cr
//' Details of the model specification as returned by `mod$mod_list`.
//'
//' @param calc_hess \[`logical(1)`\]\cr
//' Set to `TRUE` to calculate the Hessian and `FALSE` otherwise.
//'
//' @return
//' A numeric vector containing the negative log-likelihood with gradient and
//' Hessian attributes.
//'
//' @keywords internal
//'
//' @noRd

Rcpp::NumericVector ll_probit (
  Eigen::VectorXd theta,
  Rcpp::List design_diff,
  Rcpp::List mod_list,
  bool calc_hess
) {

  // ------------------------------------------------------------//
  // UNPACK MODEL SPEC
  // ------------------------------------------------------------//

  int N = design_diff.length();
  int length_theta_b = mod_list["length_theta_b"];
  int length_theta_O = mod_list["length_theta_O"];
  int length_theta_L = mod_list["length_theta_L"];
  Rcpp::IntegerVector b_count = mod_list["b_count"];
  int length_RE = b_count[1];
  int alt = mod_list["alt"];
  int length_theta = length_theta_b + length_theta_O + length_theta_L;

  int hess = 0;
  if (calc_hess){
    hess = 1;
  }
  // ------------------------------------------------------------//
  // INITIALIZE OUTPUT
  // ------------------------------------------------------------//

  Rcpp::NumericVector out;
  double ll = 0;

  // ------------------------------------------------------------//
  // RECOVER PARAMETERS
  // ------------------------------------------------------------//
  Rcpp::List syst(7);
  syst = syst_from_mod (theta, mod_list);

  // Recover the vector b
  Eigen::VectorXd b = syst["b"];
  int lthbb = b.rows ();

  // Recover the Sigma matrix
  Eigen::MatrixXd Sigma_chol_vech = syst["Sigma_chol_vech"];
  Eigen::MatrixXd Sigma = syst["Sigma"];
  int lth_Sigma_chol_vech = Sigma_chol_vech.rows();

  // Recover the matrix Omega
  Eigen::MatrixXd Omega_chol_vech = syst["Omega_chol_vech"];
  Eigen::MatrixXd Omega = syst["Omega"];
  int lth_Omega_chol_vech = Omega_chol_vech.rows();

  // ------------------------------------------------------------//
  // SETUP ONE-OFF COMPUTATIONS FOR GRADIENT AND HESSIAN
  // ------------------------------------------------------------//

  // Only for pairwise with equal number of alternatives for all deciders.
  Eigen::MatrixXd D_JOm = duplmat (length_RE);

  // Compute matrix dimensions from vech lengths.
  int m_l = -0.5 + std::sqrt (0.25 + 2 * lth_Sigma_chol_vech);
  int m_o = -0.5 + std::sqrt (0.25 + 2 * lth_Omega_chol_vech);

  Eigen::MatrixXd L_JcholL = elimmat (m_l);
  Eigen::MatrixXd D_JcholL = duplmat (m_l);
  Eigen::MatrixXd K_JcholL = commmat (m_l, m_l);

  Eigen::MatrixXd L_JcholO = elimmat (m_o);
  Eigen::MatrixXd D_JcholO = duplmat (m_o);
  Eigen::MatrixXd K_JcholO = commmat (m_o, m_o);

  // This is used later to compute the Hessian of the Cholesky, which does not
  // depend on the specific value of the parameters.
  Eigen::VectorXd HoO = Omega_chol_vech.array () + 1;
  Eigen::VectorXd HlL = Sigma_chol_vech.array () + 1;

  Eigen::MatrixXd HJ2 = J_bvOvSig_bLOLSig (
    HoO, HlL, lthbb, L_JcholO, K_JcholO, D_JcholO, L_JcholL,
    K_JcholL, D_JcholL
  );

  // ------------------------------------------------------------//
  // END of setup
  // ------------------------------------------------------------//

  // ------------------------------------------------------------//
  // INITIALIZE GRADIENT OUTPUT
  // ------------------------------------------------------------//

  Eigen::MatrixXd Hess (length_theta, length_theta);
  Hess.setZero ();
  Eigen::MatrixXd Hess_approx2 (length_theta, length_theta);
  Hess_approx2.setZero ();
  Eigen::VectorXd grad (length_theta);
  grad.setZero ();
  for (int i_n = 0; i_n < N; i_n++) { //loop over the sample

    // Initialize log-likelihood, gradient, and Hessian contributions for the
    // current decision maker.

    Eigen::VectorXd grad_n (length_theta);
    grad_n.setZero ();

    Eigen::MatrixXd Hess_n (length_theta, length_theta);
    Hess_n.setZero ();

    // Unpack the design data for the current decision maker.
    Rcpp::List data_n = design_diff[i_n];
    Rcpp::List X_n = data_n["X"];
    Eigen::VectorXd y_n = data_n["y"];

    int Tp_n = y_n.rows ();
    int alt2 = Tp_n * (alt - 1);

    // setup X, M, M_Sig.
    Eigen::MatrixXd X (alt2, lthbb);
    X.setZero ();

    Eigen::MatrixXd M_Sig (alt2, alt2);
    M_Sig.setZero ();

    // NOW LOOP OVER THE observations MATRIX
    for (int jj = 0; jj < Tp_n; jj++) { // loop jj over observations.

      // ------------------------------------------------------------//
      // Data Prep and Lin Utilities
      // ------------------------------------------------------------//
      //
      // extract the chosen alternatives for the current pair of choice
      // occasions
      int y1;
      y1 = y_n[jj];

      // extract the covariance matrix for the current pair of choice occasions
      Eigen::MatrixXd X1;
      X1 = X_n[jj];

      X.block ((jj) * (alt - 1), 0, alt - 1, lthbb) = X1;

      Eigen::MatrixXd M1 = makeM ((alt), (y1 - 1));
      Eigen::MatrixXd M1t = M1.transpose ();
      Eigen::MatrixXd Sigma_norm_1 = M1 * Sigma * M1t;

      // ------------------------------------------------------------//
      // Covariance Matrix
      // ------------------------------------------------------------//

      M_Sig.block ((jj) * (alt - 1), (jj) * (alt - 1), alt - 1, alt - 1) =
      Sigma_norm_1;
    } // end jj loop over observations.

    // find the NA elements

    //Rcpp::NumericVector x1c = X.block(0,0,X.rows(),1);
    Eigen::MatrixXd X_red = X;
    Eigen::MatrixXd M_Sig_red = M_Sig;

    Rcpp::IntegerVector ind = NA_ind (X);

    int cur = ind.size ();
    if ((cur == 1) && (ind[0] == -1)) {
      continue;
    }

    for (int ij = 0; ij < cur; ij++) {
      X_red.row (ij) = X.row (ind[ij]);
      M_Sig_red.row (ij) = M_Sig.row (ind[ij]);
    }

    for (int ij = 0; ij < cur; ij++) {
      M_Sig_red.col (ij) = M_Sig_red.col (ind[ij]);
    }

    // // unstandardized lin-utilties (upper limits of the MVNCDF)

    Eigen::VectorXd x_upper = (-1) * (X_red.block (0, 0, cur, X_red.cols ()) * b);
    //
    // // compute covariance matrix Omega of the random effects for this pair
    // of observations using the observed effects
    Eigen::MatrixXd X_RE = X_red.block (0, X_red.cols() - length_RE , cur, length_RE);
    Eigen::MatrixXd X_RE_t = X_RE.transpose ();

    // Compute combined covariance matrix
    // Changed Sig to Lambda for combined Covariance
    Eigen::MatrixXd Lambda = M_Sig_red.block (0, 0, cur, cur);
    if (length_RE != 0) {
      Lambda = Lambda + X_RE * Omega * X_RE_t;
    }

    Eigen::MatrixXd lambda = Lambda.diagonal ().array ().sqrt ();
    Eigen::MatrixXd lambda_t = lambda.transpose ();

    // ------------------------------------------------------------//
    //  Normalize upper limits and covariance matrix (to correlation matrix)
    // ------------------------------------------------------------//
    //
    Eigen::VectorXd x_norm = x_upper.array () / lambda.array ();
    Eigen::MatrixXd lambda_li_lj = lambda * lambda_t;
    Eigen::MatrixXd Lambda_cor =
    (Lambda.array () / lambda_li_lj.array ()).matrix ();

    // ------------------------------------------------------------//
    // Likelihood and Gradient
    // ------------------------------------------------------------//


    Rcpp::List outn = cal_grad_Hess(
      x_norm, Lambda_cor, Lambda, mod_list, X_red, cur, alt, ind, y_n, D_JOm,
      syst, L_JcholO, K_JcholO, D_JcholO, L_JcholL, K_JcholL, D_JcholL, HJ2,
      hess
    );
    double ll_n = outn["ll"];

    // ------------------------------------------------------------//
    // Gradient
    // ------------------------------------------------------------//
    //

    // initialize gradient contribution of the current pair
    grad_n = outn["gradient"];

    // ------------------------------------------------------------//
    // Hessian
    // ------------------------------------------------------------//

    Hess_approx2 += grad_n * grad_n.transpose ();

    Hess_n.setZero();

    if (hess == 1) {
      // update Hessian
      Hess_n = outn["hessian"];
      Hess += Hess_n;
    } //end if hessian

    ///// old part
    ll -= ll_n;
    grad -= grad_n;

  } //end N-loop

  // ------------------------------------------------------------//
  // Format the output
  // ------------------------------------------------------------//
  out = ll;
  out.attr ("gradient") = (grad.array ());
  out.attr ("hessian_approx") = (Hess_approx2.array ());
  out.attr ("Vhat") = (Hess_approx2.array ());

  if (hess == 1) {
    out.attr ("hessian") = (Hess.array ());
  }

  return out;

}

//' Calculate Approximated Composite Marginal or Full Probit Log-Likelihood
//'
//' @param theta \[`numeric()`\]\cr
//' Parameter vector.
//'
//' @param design_diff \[`list`\]\cr
//' Differenced design list as returned by
//' `choice_data$build_design(..., difference = TRUE)`.
//'
//' @param mod_list \[`list`\]\cr
//' Details of the model specification as returned by `mod$mod_list`.
//'
//' @param cml_type \[`character(1)`\]\cr
//' Defines the type of composite marginal likelihood used. Possible values:
//' - `"probit"` for full probit,
//' - `"FP"` for full pairwise,
//' - `"AP"` for adjacent pairwise,
//' - `"AC"` for cyclical adjacent pairwise.
//'
//' @param calc_hess \[`logical(1)`\]\cr
//' If `TRUE`, the returned value also has an analytical `"hessian"` attribute.
//' The outer-product approximation in `"hessian_approx"` is returned
//' regardless of this setting.
//'
//' @return
//' A numeric vector containing the negative log-likelihood with `"gradient"`,
//' `"hessian_approx"`, and `"Vhat"` attributes. If `calc_hess = TRUE`, it also
//' has an analytical `"hessian"` attribute.
//'
//' @export
//'
// [[Rcpp::export]]

Rcpp::NumericVector ll_cml (
  Eigen::VectorXd theta,
  Rcpp::List design_diff,
  Rcpp::List mod_list,
  std::string cml_type,
  bool calc_hess = false
) {

  // initialize output
  Rcpp::NumericVector out;
  int cml_pair_type = 0;

  if (cml_type == "probit") {
    out = ll_probit(theta, design_diff, mod_list, calc_hess);
  } else if (cml_type == "FP") {
    cml_pair_type = 0;
    out = ll_macml(theta, design_diff, mod_list, cml_pair_type, calc_hess);
  } else if (cml_type == "AP") {
    cml_pair_type = 2;
    out = ll_macml(theta, design_diff, mod_list, cml_pair_type, calc_hess);
  } else if (cml_type == "AC") {
    cml_pair_type = 1;
    out = ll_macml(theta, design_diff, mod_list, cml_pair_type, calc_hess);
  } else {
    Rcpp::stop(
      "Invalid cml_type '%s'. Use one of 'probit', 'FP', 'AP', or 'AC'.",
      cml_type.c_str()
    );
  }

  return out;
}

//' Calculate Approximated Probit Choice Probability
//'
//' @param X \[`matrix`\]\cr
//' Design matrix for one choice occasion, differenced with respect to the
//' chosen alternative, see `choice_data$build_design(..., difference = TRUE)`.
//'
//' @param y \[`integer(1)`\]\cr
//' Index of chosen alternative.
//'
//' @param b \[`numeric()`\]\cr
//' Coefficient vector.
//'
//' @param Sigma \[`matrix`\]\cr
//' Error covariance matrix.
//'
//' @param Omega \[`matrix`\]\cr
//' Random-effect covariance matrix.
//'
//' @return
//' A scalar log-probability for the chosen alternative.
//'
// [[Rcpp::export]]

double calc_choice_probs (
  Eigen::MatrixXd X,
  int y,
  Eigen::VectorXd b,
  Eigen::MatrixXd Sigma,
  Eigen::MatrixXd Omega
) {
  double out = 0;

  int alt = X.rows()+1;
  int alt2 = (alt - 1);
  Eigen::MatrixXd M_Sig (alt2, alt2);
  M_Sig.setZero ();

  Eigen::MatrixXd M1 = makeM ((alt), (y - 1));
  Eigen::MatrixXd M1t = M1.transpose ();
  Eigen::MatrixXd Sigma_norm_1 = M1 * Sigma * M1t;

  // ------------------------------------------------------------//
  // Covariance Matrix
  // ------------------------------------------------------------//

  M_Sig = Sigma_norm_1;

  // find the NA elements
  Eigen::MatrixXd X_red = X;
  Eigen::MatrixXd M_Sig_red = M_Sig;

  Rcpp::IntegerVector ind = NA_ind (X);

  int cur = ind.size ();
  if ((cur == 1) && (ind[0] == -1)) {
    return 0;
  }

  for (int ij = 0; ij < cur; ij++) {
    X_red.row (ij) = X.row (ind[ij]);
    M_Sig_red.row (ij) = M_Sig.row (ind[ij]);
  }

  for (int ij = 0; ij < cur; ij++) {
    M_Sig_red.col (ij) = M_Sig_red.col (ind[ij]);
  }

  // unstandardized lin-utilties (upper limits of the MVNCDF)
  Eigen::VectorXd x_upper =
    (-1) * (X_red.block (0, 0, cur, X_red.cols ()) * b);

  //
  // compute covariance matrix Omega of the random effects for this pair of
  // observations using the observed effects
  int length_RE = Omega.rows ();

  Eigen::MatrixXd X_RE = X_red.block (0, X_red.cols()-length_RE , cur, length_RE);
  Eigen::MatrixXd X_RE_t = X_RE.transpose ();
  //

  // // Compute combined covariance matrix
  // // Changed Sig to Lambda for combined Covariance
  Eigen::MatrixXd Lambda = M_Sig_red.block (0, 0, cur, cur);
  if (length_RE != 0) {
    Lambda = Lambda + X_RE * Omega * X_RE_t;
  }

  Eigen::MatrixXd lambda = Lambda.diagonal ().array ().sqrt ();
  Eigen::MatrixXd lambda_t = lambda.transpose ();

  // ------------------------------------------------------------//
  //  Normalize upper limits and covariance matrix (to correlation matrix)
  // ------------------------------------------------------------//
  //
  Eigen::VectorXd x_norm = x_upper.array () / lambda.array ();
  Eigen::MatrixXd lambda_li_lj = lambda * lambda_t;
  Eigen::MatrixXd Lambda_cor =
  (Lambda.array () / lambda_li_lj.array ()).matrix ();

  out = SJ (x_norm, Lambda_cor);

  return out;
}

//' Predict Probit Choice Probabilities
//'
//' @param theta \[`numeric()`\]\cr
//' Parameter vector.
//'
//' @param design_single \[`list`\]\cr
//' Design list (not differenced) for one decider, see
//' `choice_data$build_design(..., difference = FALSE)`.
//'
//' @param mod_list \[`list`\]\cr
//' Details of the model specification as returned by `mod$mod_list`.
//'
//' @return
//' A matrix of predicted choice probabilities according to the model. The
//' panel dimension is ignored.
//'
//' @export
//'
// [[Rcpp::export]]

Eigen::MatrixXd pred_probit_approx (
  Eigen::VectorXd theta,
  Rcpp::List design_single,
  Rcpp::List mod_list
) {

  // ------------------------------------------------------------//
  // UNPACK MODEL SPEC
  // ------------------------------------------------------------//

  int alt = mod_list["alt"];

  // ------------------------------------------------------------//
  // Extract the single design list and calculate number of choice occasions.
  // ------------------------------------------------------------//

  // Save the design list of the current decision maker.
  Rcpp::List data_n = design_single;
  Rcpp::List X_n = data_n["X"];
  Eigen::VectorXd y_n = data_n["y"];
  int Tp_n = y_n.rows ();

  // ------------------------------------------------------------//
  // INITIALIZE OUTPUT
  // ------------------------------------------------------------//

  Eigen::MatrixXd prob_predict (Tp_n, alt);
  prob_predict.setConstant (-100.0);
  // ------------------------------------------------------------//
  // RECOVER PARAMETERS
  // ------------------------------------------------------------//

  Rcpp::List syst(7);
  syst = syst_from_mod (theta, mod_list);

  // Recover the vector b
  Eigen::VectorXd b = syst["b"];

  // Recover the Sigma matrix
  Eigen::MatrixXd Sigma = syst["Sigma"];

  // Recover the matrix Omega
  Eigen::MatrixXd Omega = syst["Omega"];

  // ------------------------------------------------------------//
  // END of setup
  // ------------------------------------------------------------//

  // ------------------------------------------------------------//
  // CYCLE over choices to calculate the  probs.
  // ------------------------------------------------------------//
  for (int i_choice_1 = 0; i_choice_1 < Tp_n; i_choice_1++) {

    // extract the covariance matrix for the current pair of choice occasions
    Eigen::MatrixXd X1 = X_n[i_choice_1];
    Eigen::MatrixXd X_red = X1;

    Rcpp::IntegerVector ind = NA_ind (X1);
    int cur = ind.size ();
    if ((cur == 1) && (ind[0] == -1)) {
      continue;
    }

    // cycle over chosen alternative
    for (int ijc = 0; ijc < cur; ijc++) {

      int ij = ind[ijc];

      // diff the matrix by the chosen alternative
      Eigen::MatrixXd X1work = X_red;
      Eigen::VectorXd X1_chosen = X1work.row ((ij));

      removeRow (X1work, ij); //cpp indexing, choices start at 1

      X1work.rowwise () -= X1_chosen.transpose ();

      // X1work is prepared for probability calculation.
      double ll_n_est = 0.0;

      ll_n_est += std::exp(
        calc_choice_probs(X1work, ij + 1, b, Sigma, Omega)
      );

      // take logarithm again
      prob_predict ((i_choice_1), ij) = std::log (ll_n_est);

    } //end y_c loop over alternatives

  } //end i_choice_1 loop over choice decisions.

  // ------------------------------------------------------------//
  // undo logarithm and normalize
  // ------------------------------------------------------------//
  prob_predict = prob_predict.array ().exp ();
  Eigen::VectorXd prob_sum = prob_predict.array ().rowwise ().sum ();

  for (int y_c = 0; y_c < alt; y_c++) {
    prob_predict.col(y_c) = prob_predict.col(y_c).array () / prob_sum.array ();
  }
  // ------------------------------------------------------------//
  // Format the output
  // ------------------------------------------------------------//

  return prob_predict;

}
