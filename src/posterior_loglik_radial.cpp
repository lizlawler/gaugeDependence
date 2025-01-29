#include "gauge_functions.h"
#include "star_vol.h"
#include "lpdf.h"
#include "priors.h"
#include <RcppArmadillo.h>
#include <string>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
arma::mat calc_loglik_radial(const arma::vec& radii, const arma::vec& r0w, const arma::vec& angles, 
                             const arma::mat& posterior_params,
                             const std::string& likelihood_type,
                             const std::string& gauge_type) {
  
  // Grab appropriate log likelihood function
  radial_lpdf_fn lpdf_fn = get_lpdf_function(likelihood_type);
  
  int n_iters = posterior_params.n_rows,
    n_obs = radii.n_elem;
  
  // Container setup and initialization
  mat post_loglik = mat(n_iters, n_obs, fill::zeros);
  
  // Loop through all iterations to generate posterior loglikelihood for each iteration
  for(int iter = 0; iter < n_iters; ++iter) {
    vec current_params = posterior_params.row(iter).t();
    post_loglik.row(iter) = lpdf_fn(radii, r0w, angles,
                    current_params(0), current_params.subvec(1, current_params.n_elem - 1),
                    gauge_type).t();
  }
  
  // Return vector of pointwise posterior log likelihoods
  return post_loglik;
}  

/*** R
# post_params_gauss <- qs::qread("../gaugeDep/samplers/rcpp/radial_mcmc_fits/gauss/gauss_cens_low_20.qs")
# data <- RcppSimdJson::fload("../gaugeDep/data/gauss/low_20.json")
# w <- data$W
# r <- data$R
# r0w <- data$r0_w
# idx <- data$idx
# result_gauss <- calc_loglik_radial(radii = r[idx], r0w = r0w[idx], angles = w[idx],
#                                    posterior_params = post_params_gauss$samples[,1:2],
#                                    likelihood_type = "trunc",
#                                    gauge_type = "gauss")
*/