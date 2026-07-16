// prox_admm_optimized_fixed.cpp
// CORRECTED: Matches R ADMM exactly with warm-starting and comprehensive diagnostics

#include <RcppEigen.h>
#include <cmath>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::SimplicialLLT;
using Eigen::Map;

typedef Eigen::Triplet<double> T;

//' Build k-th order difference matrix (sparse)
SparseMatrix<double> build_diff_matrix(int n, int k) {
  SparseMatrix<double> D(n, n);
  D.setIdentity();
  
  for (int order = 0; order <= k; order++) {
    int rows = D.rows() - 1;
    int cols = D.cols();
    
    std::vector<T> triplets;
    triplets.reserve(rows * 2);
    
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        double val = D.coeff(i+1, j) - D.coeff(i, j);
        if (std::abs(val) > 1e-14) {
          triplets.push_back(T(i, j, val));
        }
      }
    }
    
    SparseMatrix<double> D_new(rows, cols);
    D_new.setFromTriplets(triplets.begin(), triplets.end());
    D = D_new;
  }
  
  D.makeCompressed();
  return D;
}

//' Soft-thresholding operator
VectorXd soft_threshold(const VectorXd& x, double lambda) {
  VectorXd result(x.size());
  for (int i = 0; i < x.size(); i++) {
    if (x(i) > lambda) {
      result(i) = x(i) - lambda;
    } else if (x(i) < -lambda) {
      result(i) = x(i) + lambda;
    } else {
      result(i) = 0.0;
    }
  }
  return result;
}

VectorXd box_project(const VectorXd& x, const VectorXd& lower, const VectorXd& upper) {
  VectorXd result(x.size());
  for (int i = 0; i < x.size(); i++) {
    double val = x(i);
    if (!std::isinf(lower(i))) {
      val = std::max(val, lower(i));
    }
    if (!std::isinf(upper(i))) {
      val = std::min(val, upper(i));
    }
    result(i) = val;
  }
  return result;
}

//' Proximal ADMM with Warm-Starting and Comprehensive Diagnostics
//' 
//' @param w Langevin proposal point
//' @param lambda Penalty parameter
//' @param gamma Proximal scaling parameter
//' @param lower_bound Lower box constraints
//' @param upper_bound Upper box constraints
//' @param k Trend filtering order (2 for piecewise quadratic)
//' @param rho ADMM penalty parameter
//' @param max_iter Maximum ADMM iterations
//' @param tol Convergence tolerance
//' @param theta_init Warm start for theta (NULL for cold start)
//' @param v_init Warm start for v (z variable)
//' @param w_init_param Warm start for w variable
//' @param u1_init Warm start for dual u1
//' @param u2_init Warm start for dual u2
//' @param verbose Print iteration info
//' @param print_every Print diagnostics every N iterations (0 = no printing)
// [[Rcpp::export]]
Rcpp::List prox_admm_optimized_fixed(
    const Eigen::VectorXd& w,
    double lambda,
    double gamma,
    const Eigen::VectorXd& lower_bound,
    const Eigen::VectorXd& upper_bound,
    int k = 2,
    double rho = 1.0,
    int max_iter = 5000,
    double tol = 1e-5,
    Rcpp::Nullable<Eigen::VectorXd> theta_init = R_NilValue,
    Rcpp::Nullable<Eigen::VectorXd> v_init = R_NilValue,
    Rcpp::Nullable<Eigen::VectorXd> w_init_param = R_NilValue,
    Rcpp::Nullable<Eigen::VectorXd> u1_init = R_NilValue,
    Rcpp::Nullable<Eigen::VectorXd> u2_init = R_NilValue,
    bool verbose = false,
    int print_every = 100
) {
  
  int n = w.size();
  double rho1 = rho;
  double rho2 = rho;
  
  if (lower_bound.size() != n || upper_bound.size() != n) {
    Rcpp::stop("lower_bound and upper_bound must have same length as w");
  }
  
  // Build D matrix
  SparseMatrix<double> D = build_diff_matrix(n, k);
  int m = D.rows();
  
  // Build system matrix: A = rho1*D^T*D + rho2*I
  SparseMatrix<double> DtD = D.transpose() * D;
  SparseMatrix<double> A = rho1 * DtD;
  
  for (int i = 0; i < n; i++) {
    A.coeffRef(i, i) += rho2;
  }
  A.makeCompressed();
  
  // Factor the matrix ONCE
  SimplicialLLT<SparseMatrix<double>> chol;
  chol.compute(A);
  
  if (chol.info() != Eigen::Success) {
    Rcpp::stop("Matrix factorization failed");
  }
  
  // Initialize variables
  VectorXd theta = theta_init.isNotNull() ? Rcpp::as<VectorXd>(theta_init) : w;
  VectorXd z = v_init.isNotNull() ? Rcpp::as<VectorXd>(v_init) : (D * theta);
  VectorXd w_var = w_init_param.isNotNull() ? Rcpp::as<VectorXd>(w_init_param) : theta;
  VectorXd u1 = u1_init.isNotNull() ? Rcpp::as<VectorXd>(u1_init) : VectorXd::Zero(m);
  VectorXd u2 = u2_init.isNotNull() ? Rcpp::as<VectorXd>(u2_init) : VectorXd::Zero(n);
  
  double r_primal_1, r_primal_2, r_dual_1, r_dual_2;
  double theta_L1, z_L1, w_L1, w_theta_L1;
  VectorXd z_old, w_var_old, theta_old;
  VectorXd w_init = w_var;
  
  // ADMM LOOP
  int iter;
  for (iter = 1; iter <= max_iter; iter++) {
    
    if (iter % 100 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    theta_old = theta;
    z_old = z;
    w_var_old = w_var;
    
    // Theta update
    VectorXd Dt_z_minus_u1 = D.transpose() * (z - u1);
    VectorXd RHS = rho1 * Dt_z_minus_u1 + rho2 * (w_var - u2);
    theta = chol.solve(RHS);
    theta_L1 = (theta - theta_old).cwiseAbs().sum();
    
    // Z update
    VectorXd Dtheta = D * theta;
    VectorXd z_tilde = Dtheta + u1;
    double threshold = lambda / rho1;
    z = soft_threshold(z_tilde, threshold);
    z_L1 = (z - z_old).cwiseAbs().sum();
    
    // W_var update
    double alpha = (1.0 / gamma) + rho2;
    VectorXd c_vec = (1.0 / gamma) * w + rho2 * (theta + u2);
    VectorXd w_tilde = c_vec / alpha;
    w_var = box_project(w_tilde, lower_bound, upper_bound);
    w_L1 = (w_var - w_var_old).cwiseAbs().sum();
    w_theta_L1 = (w_var - theta).cwiseAbs().sum();
    
    // Dual updates
    u1 += Dtheta - z;
    u2 += theta - w_var;
    
    // Convergence check
    r_primal_1 = (Dtheta - z).norm();
    r_primal_2 = (theta - w_var).norm();
    r_dual_1 = (rho1 * D.transpose() * (z - z_old)).norm();
    r_dual_2 = (rho2 * (w_var - w_init)).norm();
    
    double max_residual = std::max({r_primal_1, r_primal_2, r_dual_1, r_dual_2});
    
    if (max_residual < tol) {
      break;
    }
    
    w_init = w_var;
    
    // Print progress (only if verbose AND print_every > 0)
    if (verbose && print_every > 0 && iter % print_every == 0) {
      Rcpp::Rcout << "      [ADMM " << iter << "]"
                  << " Resid: r_p1=" << r_primal_1 
                  << " r_p2=" << r_primal_2 
                  << " r_d1=" << r_dual_1
                  << " r_d2=" << r_dual_2
                  << " | L1: theta=" << theta_L1
                  << " z=" << z_L1
                  << " w=" << w_L1
                  << " w-theta=" << w_theta_L1 << std::endl;
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("theta") = theta,
    Rcpp::Named("v") = z,
    Rcpp::Named("w") = w_var,
    Rcpp::Named("u1") = u1,
    Rcpp::Named("u2") = u2,
    Rcpp::Named("r_primal_1") = r_primal_1,
    Rcpp::Named("r_primal_2") = r_primal_2,
    Rcpp::Named("r_dual_1") = r_dual_1,
    Rcpp::Named("r_dual_2") = r_dual_2,
    Rcpp::Named("theta_L1") = theta_L1,
    Rcpp::Named("z_L1") = z_L1,
    Rcpp::Named("w_L1") = w_L1,
    Rcpp::Named("w_theta_L1") = w_theta_L1,
    Rcpp::Named("iterations") = iter,
    Rcpp::Named("converged") = (iter <= max_iter)
  );
}
