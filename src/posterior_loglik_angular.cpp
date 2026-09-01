// =============================================================================
// posterior_loglik_angular.cpp
//
// Computes the posterior pointwise log-likelihood matrix for the angular
// density across all MCMC iterations. Exported to R as calc_loglik_angular(),
// called by the R wrapper angular_loglik().
//
// For each MCMC iteration (row of posterior_params), evaluates the angular
// log-density at each observed angle. Returns an (n_iter x n_obs) matrix
// for use with loo::loo() or loo::stacking_weights().
// =============================================================================

#include "gauge_functions.h"
#include "star_vol.h"
#include "lpdf.h"
#include "priors.h"
#include <RcppArmadillo.h>
#include <string>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
arma::mat calc_loglik_angular(const arma::vec& angles,
                              const arma::mat& grid_x, const int& dim,
                              const arma::mat& posterior_params,
                              const std::string& gauge_type) {
  int n_iters = posterior_params.n_rows,
      n_obs   = angles.n_elem;

  mat post_loglik = mat(n_iters, n_obs, fill::zeros);

  for (int iter = 0; iter < n_iters; ++iter) {
    // Extract current iteration's parameter vector (column vector)
    vec current_params = posterior_params.row(iter).t();
    post_loglik.row(iter) = angular_lpdf(angles, current_params,
                                         grid_x, dim, gauge_type).t();
  }

  return post_loglik;
}
