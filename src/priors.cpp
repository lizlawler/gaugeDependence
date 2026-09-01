// =============================================================================
// priors.cpp
//
// Log-prior functions for the radial and angular dependence parameters.
//
// Priors (chosen to be weakly informative):
//   alpha (Gamma shape):    Gamma(4, 2)   — mean 2, positive support
//   dep (non-Dirichlet):    Uniform(0, 1) for angular;
//                           Gamma(4, 2) + Uniform(0, 1) for radial
//   Dirichlet theta1:       half-t(4) scaled by 4 (positive support, heavy-tailed)
//   Dirichlet theta2:       half-t(4) scaled by 2 (positive support, heavy-tailed)
//
// Returns -Inf for parameter values outside the support (used in MH acceptance
// step to immediately reject proposals that violate constraints).
// =============================================================================

#include "priors.h"
#include <RcppArmadillo.h>

// Log-prior for the radial parameter vector: params = c(alpha, dep[, theta2])
// For non-Dirichlet gauges: params has length 2, c(alpha, dep)
// For Dirichlet gauge:      params has length 3, c(alpha, theta1, theta2)
double prior_fn_radii(const arma::vec& params) {
  if (params.n_elem <= 2) {
    // alpha ~ Gamma(4, 0.5),  dep ~ Uniform(0, 1)
    return R::dgamma(params(0), 4, 0.5, true) + R::dunif(params(1), 0, 1, true);
  } else {
    // Dirichlet gauge: theta1, theta2 must be positive
    if (params(1) < 0 || params(2) < 0) return -arma::datum::inf;
    // alpha ~ Gamma(4, 0.5),  theta1 ~ half-t(4) * 4,  theta2 ~ half-t(4) * 2
    return R::dgamma(params(0), 4, 0.5, true)
         + R::dt(params(1) / 4, 4, true)
         + R::dt(params(2) / 2, 4, true);
  }
}

// Log-prior for the angular parameter vector: params = c(dep) or c(theta1, theta2)
// For non-Dirichlet gauges: params has length 1, dep ~ Uniform(0, 1)
// For Dirichlet gauge:      params has length 2, theta1 and theta2 ~ half-t(4)
double prior_fn_angles(const arma::vec& params) {
  if (params.n_elem == 1) {
    // dep ~ Uniform(0, 1)
    return R::dunif(params(0), 0, 1, true);
  } else {
    // Dirichlet gauge: theta1, theta2 must be positive
    if (params(0) < 0 || params(1) < 0) return -arma::datum::inf;
    // theta1 ~ half-t(4) * 4,  theta2 ~ half-t(4) * 2
    return R::dt(params(0) / 4, 4, true) + R::dt(params(1) / 2, 4, true);
  }
}
