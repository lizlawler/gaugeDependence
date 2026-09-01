// =============================================================================
// posterior_loglik_radial.cpp
//
// Computes the posterior pointwise log-likelihood matrix for the radial model
// across all MCMC iterations. Exported to R as calc_loglik_radial(), called
// by the R wrapper radial_loglik().
//
// For each MCMC iteration (row of posterior_params), evaluates the radial
// log-density at each observation. The first column of posterior_params is
// the Gamma shape parameter alpha; subsequent columns are the gauge dependence
// parameter(s). Returns an (n_iter x n_obs) matrix for use with loo::loo()
// or loo::stacking_weights().
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
arma::mat calc_loglik_radial(const arma::vec& radii, const arma::vec& r0w,
                             const arma::vec& angles,
                             const arma::mat& posterior_params,
                             const std::string& likelihood_type,
                             const std::string& gauge_type) {
  // Select the appropriate pointwise log-density function
  radial_lpdf_fn lpdf_fn = get_lpdf_function(likelihood_type);

  int n_iters = posterior_params.n_rows,
      n_obs   = radii.n_elem;

  mat post_loglik = mat(n_iters, n_obs, fill::zeros);

  for (int iter = 0; iter < n_iters; ++iter) {
    vec current_params = posterior_params.row(iter).t();
    // First element is alpha (shape); remainder are dependence parameter(s)
    post_loglik.row(iter) = lpdf_fn(radii, r0w, angles,
                                    current_params(0),
                                    current_params.subvec(1, current_params.n_elem - 1),
                                    gauge_type).t();
  }

  return post_loglik;
}
