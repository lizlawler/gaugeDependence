// =============================================================================
// priors.h
//
// Declarations for the log-prior functions for the radial and angular
// dependence parameters. Priors are:
//   alpha (radial shape): Gamma(4, 2)
//   dep (non-Dirichlet):  Uniform(0, 1) for angular; Gamma(4, 2) + Uniform(0,1) for radial
//   Dirichlet (theta1, theta2): half-t(4) scaled by 4 and 2 respectively
// =============================================================================
#ifndef PRIOR_H
#define PRIOR_H

#include <RcppArmadillo.h>

// Log-prior for the radial parameter vector (alpha, dep[, theta2])
double prior_fn_radii(const arma::vec& params);

// Log-prior for the angular parameter vector (dep[, theta2] for Dirichlet)
double prior_fn_angles(const arma::vec& params);

#endif  // PRIOR_H
