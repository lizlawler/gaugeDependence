#include "loglik.h"
#include "helpers.h"
#include <RcppArmadillo.h>

// Radial parameters priors
double prior_fn_radii(arma::vec const& params) {
  return R::dgamma(params(0), 4, 0.5, true) + R::dunif(params(1), 0, 1, true);
}

// ANgular parameter prior
double prior_fn_angles(double const& params) {
  return R::dunif(params, 0, 1, true);
}