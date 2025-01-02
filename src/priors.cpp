#include "priors.h"
#include <RcppArmadillo.h>

// Radial parameters priors
double prior_fn_radii(const arma::vec& params) {
  if(params.n_elem <= 2) {
    return R::dgamma(params(0), 4, 0.5, true) + R::dunif(params(1), 0, 1, true);
  } else { // accommodate multiple parameters in dirichlet gauge
    if(params(1) < 0 || params(2) < 0) {
      return -arma::datum::inf;
    } else {
      return R::dgamma(params(0), 4, 0.5, true) + R::dt(params(1)/4, 4, true) + R::dt(params(2)/2, 4, true);
    }
  }
}

// Angular parameter priors
double prior_fn_angles(const arma::vec& params) {
  if(params.n_elem == 1) {
    return R::dunif(params(0), 0, 1, true);
  } else { // accommodate multiple parameters in dirichlet gauge
    if(params(0) < 0 || params(1) < 0) {
      return -arma::datum::inf;
    } else {
      return R::dt(params(0)/4, 4, true) + R::dt(params(1)/2, 4, true);
    }
  }
}