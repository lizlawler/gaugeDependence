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
    n_obs = angles.n_elem;
  
  // Container setup and initialization
  mat post_loglik = mat(n_iters, n_obs, fill::zeros);
  
  // Loop through all iterations to generate posterior loglikelihood for each iteration
  for(int iter = 0; iter < n_iters; ++iter) {
    vec current_params = posterior_params.row(iter).t();
    post_loglik.row(iter) = angular_lpdf(angles, 
                    current_params, 
                    grid_x, dim, 
                    gauge_type).t();
  }
  
  // Return vector of pointwise posterior log likelihoods
  return post_loglik;
}  

/*** R
# post_params <- qs::qread("../gaugeDep/samplers/rcpp/angular_vol_mcmc_fits/gauss/gauss_high_1.qs")
# data <- RcppSimdJson::fload("../gaugeDep/data/gauss/low_20.json")
# w <- data$W
# x1 <- seq(0,1,length.out = 100)
# grid_x <- expand.grid(x1, x1)
# result <- calc_loglik_angular(angles = w,
#                               grid_x = as.matrix(grid_x),
#                               dim = 2,
#                               posterior_params = as.matrix(post_params$samples[,"dep"]),
#                               gauge_type = "gauss")
*/