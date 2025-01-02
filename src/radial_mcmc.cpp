#include "gauge_functions.h"
#include "star_vol.h"
#include "loglik.h"
#include "priors.h"
#include <RcppArmadillo.h>
#include <string>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
Rcpp::List radial_adaptive_mh(const arma::vec& radii, const arma::vec& r0w, const arma::vec& angles, 
                              const arma::rowvec& starting_theta,
                              const std::string& likelihood_type,
                              const std::string& gauge_type,
                              Rcpp::Nullable<Rcpp::NumericMatrix> prop_Sigma_ = R_NilValue,
                              const int& n_updates = 10000, const int& update_freq = 100,
                              const int& n_burnin = 0, const int& n_thin = 1,
                              const bool& adapt_cov = false,
                              const double& r_opt = 0.234, 
                              const double& c0 = 10, const double& c1 = 0.8, 
                              const int& K = 10) {
  
  // Grab appropriate log likelihood function
  radial_loglik_fn loglik_fn = get_loglik_function(likelihood_type);
  
  double eps = 0.001;
  int k = 3, // iteration offest
    p = starting_theta.size();
  
  // Initialize sigma_m to rule of thumb
  double sigma_m = pow(2.4 / (double)p, 2);
  
  // Radii proposal check
  mat prop_Sigma;
  if (prop_Sigma_.isNotNull()) {
    prop_Sigma = Rcpp::as<arma::mat>(prop_Sigma_);
    if (prop_Sigma.size() != p * p || !prop_Sigma.is_symmetric()) {
      Rcpp::Rcout << "Proposal matrix is invalid, using identity matrix." << std::endl;
      prop_Sigma = arma::eye(p, p);  // Use identity matrix
    }
  } else {
    Rcpp::Rcout << "Proposal matrix for radial parameters not provided, using identity matrix." << std::endl;
    prop_Sigma = eye(p, p);  // Use identity matrix
  } 
  
  // Initialize cholesky decomp of proposal covariance matrix
  mat prop_C = chol(prop_Sigma);
  
  // Trace setup and initialization
  mat trace = mat(n_updates, p, fill::zeros);
  vec sigma_m_trace = vec(n_updates, fill::zeros), 
    jump_trace = vec(n_updates, fill::zeros);
  
  trace.row(0) = starting_theta;
  
  // Initialize MH
  vec theta = starting_theta.t();
  double prior = prior_fn_radii(theta);
  double likelihood = loglik_fn(radii, r0w, angles,
                                theta(0), theta.subvec(1, theta.n_elem - 1), 
                                gauge_type);
  
  // Main loop of MH
  for(int iter = 1; iter < n_updates; ++iter) {
    // Radial update
    vec theta_star = theta + sigma_m * (randn<rowvec>(p) * prop_C).t();
    double prior_star = prior_fn_radii(theta_star);
    if(prior_star != -datum::inf) {
      double likelihood_star = loglik_fn(radii, r0w, angles,
                                         theta_star(0), theta_star.subvec(1, theta_star.n_elem - 1),
                                         gauge_type);
      if ((prior_star + likelihood_star - prior - likelihood) > log(randu())) {
        theta = theta_star;
        prior = prior_star;
        likelihood = likelihood_star;
        jump_trace(iter) = 1;
      }
    }
    
    // Adaptive update
    if(((iter + 1) % K) == 0) {
      double gamma1 = c0 / pow((((iter + 1) / K) + k), c1), 
        gamma2 = 1 / pow((((iter + 1) / K) + k), c1);
      
      sigma_m = exp(log(sigma_m) + 
        gamma1 * (mean(jump_trace.subvec(((iter + 1) - K), iter)) - r_opt));
      
      if(adapt_cov) {
        prop_Sigma += gamma2 * (cov(trace.rows(iter - K + 1, iter)) - prop_Sigma);
        
        // Attempt Cholesky decomp and adjust if necessary
        try {
          prop_C = chol(prop_Sigma);
        } catch (std::runtime_error& e) {
          Rcpp::Rcout << "Proposal covariance matrix became non-positive definite, adjusting with epsilon." << std::endl;
          prop_Sigma += eps * eye(p, p);  
          prop_C = chol(prop_Sigma);
        }
        
        prop_C = chol(prop_Sigma); 
      }
    } // End adaptation phase
    
    // Update trace/sample objects
    trace.row(iter) = theta.t();
    sigma_m_trace(iter) = sigma_m;
    
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
  mat result_matrix(thinning_idx.size(), p + 1);
  Rcpp::CharacterVector col_names;
  
  // Fill the matrix with data from each trace
  result_matrix.col(0) = trace.submat(thinning_idx, uvec{0});  
  result_matrix.col(1) = trace.submat(thinning_idx, uvec{1});
  if(p > 2) { // account for Dirichlet scenario
    result_matrix.col(2) = trace.submat(thinning_idx, uvec{2});
    result_matrix.col(3) = sigma_m_trace.elem(thinning_idx); 
    col_names = {"alpha", "theta1", "theta2", "sigma_m"};
  } else {
    result_matrix.col(2) = sigma_m_trace.elem(thinning_idx); 
    col_names = {"alpha", "dep", "sigma_m"};
  }
  
  // Convert arma::mat to Rcpp::NumericMatrix for column name assignment
  Rcpp::NumericMatrix r_result_matrix = Rcpp::wrap(result_matrix);
  
  // Assign column names (again, for use later in 'coda')
  Rcpp::colnames(r_result_matrix) = col_names;
  
  // Return list with the samples and acceptance probabilities
  return Rcpp::List::create(
    Rcpp::Named("samples") = r_result_matrix,
    Rcpp::Named("acc_prob") = mean(jump_trace)
  );
} 

/*** R
# # data <- RcppSimdJson::fload("../gaugeDep/data/gauss/high_100.json")
# w <- data$W
# r <- data$R
# r0w <- data$r0_w
# idx <- data$idx
# trunc_results <- radial_adaptive_mh(radii = r[idx], r0w = r0w[idx], angles = w[idx],
#                               starting_theta = c(rgamma(1, 4, 2), runif(1)),
#                               # starting_theta = c(rgamma(1, 4, 2), abs(rt(1, 4,ncp = 0)), abs(rt(1, 4,ncp = 0))),
#                               likelihood_type = "trunc",
#                               gauge_type = "gauss",
#                               n_updates = 15000, update_freq = 250,adapt_cov = TRUE)
# # 
# cens_results <- radial_adaptive_mh(radii = r, r0w = r0w, angles = w,
#                                     starting_theta = c(rgamma(1, 4, 2), runif(1)),
#                                     likelihood_type = "cens",
#                                     gauge_type = "gauss",
#                                     n_updates = 15000, update_freq = 250,adapt_cov = TRUE)
# # # 
# plot(trunc_results$samples[,"dep"], type = "l")
# plot(density(trunc_results$samples[,"dep"]))
# plot(cens_results$samples[,"dep"], type = "l")
# plot(density(cens_results$samples[,"dep"]))
*/