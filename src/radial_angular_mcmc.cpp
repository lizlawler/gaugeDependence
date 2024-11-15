#include "gauge_functions.h"
#include "helpers.h"
#include "loglik.h"
#include "priors.h"
#include <RcppArmadillo.h>
#include <string>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
Rcpp::List adaptive_mh_rad_ang_vol(const arma::vec& radii, const arma::vec& angles, const arma::vec& threshold, 
                                   const arma::vec& sum_term, const arma::vec& sqrt_term, const int& dim,
                                   const arma::rowvec& starting_theta_rad,
                                   const double& starting_theta_ang,
                                   std::string rad_gauge_type,
                                   std::string ang_gauge_type,
                                   Rcpp::Nullable<Rcpp::NumericMatrix> prop_Sigma_radii_ = R_NilValue,
                                   Rcpp::Nullable<double> prop_Sigma_angles_ = R_NilValue,
                                   const int& n_updates = 10000, const int& update_freq = 100,
                                   const int& n_burnin = 0, const int& n_thin = 1,
                                   const bool& adapt_cov = false,
                                   const double& r_opt = 0.234, 
                                   const double& c0 = 10, const double& c1 = 0.8, 
                                   const int& K = 10) {
  
  double eps = 0.001;
  int k = 3, // iteration offest
    p_rad = starting_theta_rad.size();
  
  // Initialize sigma_m to rule of thumb
  double sigma_m_rad = pow(2.4 / (double)p_rad, 2), sigma_m_ang = pow(2.4, 2);
  
  // Radii proposal check
  mat prop_Sigma_radii;
  if (prop_Sigma_radii_.isNotNull()) {
    prop_Sigma_radii = Rcpp::as<arma::mat>(prop_Sigma_radii_);
    if (prop_Sigma_radii.size() != p_rad * p_rad || !prop_Sigma_radii.is_symmetric()) {
      Rcpp::Rcout << "Proposal matrix is invalid, using identity matrix." << std::endl;
      prop_Sigma_radii = arma::eye(p_rad, p_rad);  // Use identity matrix
    }
  } else {
    Rcpp::Rcout << "Proposal matrix for radial parameters not provided, using identity matrix." << std::endl;
    prop_Sigma_radii = eye(p_rad, p_rad);  // Use identity matrix
  }

  // Angular proposal check
  double prop_Sigma_angles = prop_Sigma_angles_.isNotNull() 
    ? Rcpp::as<double>(prop_Sigma_angles_)
      : (Rcpp::Rcout << "Proposal variance for angular parameter not provided, using identity.\n", 1.0);

  // Initialize cholesky decomp of proposal covariance matrix (and square root of variance for angles)
  mat prop_rad_C = chol(prop_Sigma_radii);
  double prop_ang_C = sqrt(prop_Sigma_angles);
  
  // Trace setup and initialization
  mat trace_rad = mat(n_updates, p_rad, fill::zeros);
  vec sigma_m_rad_trace = vec(n_updates, fill::zeros), 
    jump_rad_trace = vec(n_updates, fill::zeros);
  
  vec trace_ang = vec(n_updates, fill::zeros), 
    sigma_m_ang_trace = vec(n_updates, fill::zeros),
    jump_ang_trace = vec(n_updates, fill::zeros);
  
  trace_rad.row(0) = starting_theta_rad;
  trace_ang(0) = starting_theta_ang;
  
  // Initialize MH
  vec theta_rad = starting_theta_rad.t();
  double theta_ang = starting_theta_ang;
  double prior_rad = prior_fn_radii(theta_rad), prior_ang = prior_fn_angles(theta_ang);
  double likelihood_rad = radial_cens_loglik(radii, theta_rad, threshold, angles, rad_gauge_type),
    likelihood_ang = angular_loglik(angles, theta_ang, sum_term, sqrt_term, dim, ang_gauge_type);
  
  // Main loop of MH
  for(int iter = 1; iter < n_updates; ++iter) {
    // Radial update
    vec theta_rad_star = theta_rad + sigma_m_rad * (randn<rowvec>(p_rad) * prop_rad_C).t();
    double prior_rad_star = prior_fn_radii(theta_rad_star);
    if(prior_rad_star != -datum::inf) {
      double likelihood_rad_star = radial_cens_loglik(radii, theta_rad_star, threshold, angles, rad_gauge_type);
      if ((prior_rad_star + likelihood_rad_star - prior_rad - likelihood_rad) > log(randu())) {
        theta_rad = theta_rad_star;
        prior_rad = prior_rad_star;
        likelihood_rad = likelihood_rad_star;
        jump_rad_trace(iter) = 1;
      }
    }
  
    // Angular update
    double theta_ang_star = theta_ang + sigma_m_ang * (randn() * prop_ang_C);
    double prior_ang_star = prior_fn_angles(theta_ang_star);
    if(prior_ang_star != -datum::inf) {
      double likelihood_ang_star = angular_loglik(angles, theta_ang_star, sum_term, sqrt_term, dim, ang_gauge_type);
      if ((prior_ang_star + likelihood_ang_star - prior_ang - likelihood_ang) > log(randu())) {
        theta_ang = theta_ang_star;
        prior_ang = prior_ang_star;
        likelihood_ang = likelihood_ang_star;
        jump_ang_trace(iter) = 1;
      }
    }
  
    // Adaptive update
    if(((iter + 1) % K) == 0) {
      double gamma1 = c0 / pow((((iter + 1) / K) + k), c1), 
        gamma2 = 1 / pow((((iter + 1) / K) + k), c1);
      
      sigma_m_rad = exp(log(sigma_m_rad) + 
        gamma1 * (mean(jump_rad_trace.subvec(((iter + 1) - K), iter)) - r_opt));
      sigma_m_ang = exp(log(sigma_m_ang) + 
        gamma1 * (mean(jump_ang_trace.subvec(((iter + 1) - K), iter)) - r_opt));
      
      if(adapt_cov) {
        prop_Sigma_radii += gamma2 * (cov(trace_rad.rows(iter - K + 1, iter)) - prop_Sigma_radii);
        prop_Sigma_angles += gamma2 * (var(trace_ang.subvec(iter - K + 1, iter)) - prop_Sigma_angles);
        
        // Attempt Cholesky decomp and adjust if necessary
        try {
          prop_rad_C = chol(prop_Sigma_radii);
        } catch (std::runtime_error& e) {
           Rcpp::Rcout << "Proposal covariance matrix became non-positive definite, adjusting with epsilon." << std::endl;
          prop_Sigma_radii += eps * eye(p_rad, p_rad);  
          prop_rad_C = chol(prop_Sigma_radii);
        }
      
        prop_rad_C = chol(prop_Sigma_radii);
        prop_ang_C = sqrt(prop_Sigma_angles);
      }
    } // End adaptation phase
  
    // Update trace/sample objects
    trace_rad.row(iter) = theta_rad.t();
    trace_ang(iter) = theta_ang;
    sigma_m_rad_trace(iter) = sigma_m_rad;
    sigma_m_ang_trace(iter) = sigma_m_ang;
    
    // Report progress to console
    if ((iter + 1) % update_freq == 0) {
      int progress = (100 * (iter + 1)) / n_updates;  // Calculate percentage progress
      Rcpp::Rcout << "\rProgress: |" 
                  << std::string(progress / 2, '-')  // Add `-` for progress
                  << std::string(50 - (progress / 2), ' ')  // Add spaces for remaining part
                  << "| " << progress << "% completed." << std::flush;
      Rcpp::checkUserInterrupt();  // Allow user interruption
    }
  }

  // Handle burn-in removal; ensure burnin <= n_updates
  int start = std::min(n_burnin, n_updates);  
  
  // Handle thinning; account for when total iterations isn't an integer multiple of thinning interval
  uvec thinning_idx = (n_updates % n_thin == 0) 
    ? regspace<uvec>(start, n_thin, n_updates - 1) 
      : (Rcpp::Rcout << "\nWarning: n_updates is not an integer multiple of n_thin, returning all iterations.", regspace<uvec>(start, n_updates - 1)); 
  
  // Re-formatting the way results are saved for easy use with 'coda' package in R
  // Create a matrix for the parameter traces; may need to adjust "5" to a different number
  mat result_matrix(thinning_idx.size(), 5);
  
  // Fill the matrix with data from each trace
  result_matrix.col(0) = trace_rad.submat(thinning_idx, uvec{0});  
  result_matrix.col(1) = trace_rad.submat(thinning_idx, uvec{1});  
  result_matrix.col(2) = trace_ang.elem(thinning_idx);         
  result_matrix.col(3) = sigma_m_rad_trace.elem(thinning_idx); 
  result_matrix.col(4) = sigma_m_ang_trace.elem(thinning_idx); 
  
  // Convert arma::mat to Rcpp::NumericMatrix for column name assignment
  Rcpp::NumericMatrix r_result_matrix = Rcpp::wrap(result_matrix);
  
  // Assign column names (again, for use later in 'coda')
  Rcpp::CharacterVector col_names = {"alpha", "dep_r", "dep_w", "sigma_m_radii", "sigma_m_angles"};
  Rcpp::colnames(r_result_matrix) = col_names;
  
  // Return list with the samples and acceptance probabilities
  return Rcpp::List::create(
    Rcpp::Named("samples") = r_result_matrix,
    Rcpp::Named("acc_prob_radii") = mean(jump_rad_trace),
    Rcpp::Named("acc_prob_angles") = mean(jump_ang_trace)
  );
}
