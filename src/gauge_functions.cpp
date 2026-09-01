// =============================================================================
// gauge_functions.cpp
//
// Implementations of the six bivariate gauge function families. Each gauge
// function g(w1, w2; theta) defines the unit star body {x : g(x) <= 1}, which
// determines the shape of the extremal dependence. The radial model is
// R | W ~ Gamma(2, g(W; theta)) above the threshold r0(W).
//
// All six functions are exported to R via [[Rcpp::export]] and appear in
// RcppExports.R. They also serve as internal components of the MCMC samplers
// and log-likelihood calculations via the get_gauge_function() dispatcher.
//
// Gauge families implemented:
//   gauss       - Gaussian (asymptotic independence); dep in (0, 1)
//   inv_log     - Inverted logistic (asymptotic dependence); dep in (0, 1)
//   rectangular - Rectangular / logistic-max; dep in (0, 1)
//   logistic    - Logistic (Gumbel); dep in (0, 1]
//   asym_log    - Asymmetric logistic; dep in (0, 1]
//   dirichlet   - Dirichlet; dep = c(theta1, theta2) with theta1, theta2 > 0
// =============================================================================

#include "gauge_functions.h"
#include <RcppArmadillo.h>
#include <string>

// -----------------------------------------------------------------------------
// Gauge function definitions
// -----------------------------------------------------------------------------

// Gaussian gauge: g(w1, w2; rho) = (w1 + w2 - 2*rho*sqrt(w1*w2)) / (1 - rho^2)
// Corresponds to Gaussian copula with correlation rho. dep in (-1, 1);
// in practice restricted to (0, 1) via the prior.
// [[Rcpp::export]]
arma::vec gauss_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep) {
  double dep_value = dep(0);
  return (w1 + w2 - 2 * dep_value * sqrt(w1 % w2)) / (1 - pow(dep_value, 2));
}

// Inverted logistic gauge: g(w1, w2; r) = (w1^(1/r) + w2^(1/r))^r
// Corresponds to inverted logistic (Gumbel survival) dependence. dep in (0, 1].
// [[Rcpp::export]]
arma::vec inv_log_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep) {
  double r_inv = 1 / dep(0);
  return pow((pow(w1, r_inv) + pow(w2, r_inv)), dep(0));
}

// Rectangular (logistic-max) gauge
// dep in (0, 1]; dep = 1 gives independence, dep -> 0 gives complete dependence.
// [[Rcpp::export]]
arma::vec rectangular_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep) {
  double dep_value = dep(0);
  arma::mat w_max(w1.size(), 3);
  w_max.col(0) = (w1 - w2) / dep_value;
  w_max.col(1) = (w2 - w1) / dep_value;
  w_max.col(2) = (w1 + w2) / (2 - dep_value);
  return max(w_max, 1);
}

// Logistic (Gumbel) gauge: g(w1, w2; r) = r^{-1} * max(w1,w2) + (1 - r^{-1}) * min(w1,w2)
// dep in (0, 1]; dep = 1 gives complete dependence (independence limit as dep -> 0).
// [[Rcpp::export]]
arma::vec logistic_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep) {
  double r_inv = 1 / dep(0);
  arma::mat w1_and_2(w1.size(), 2);
  w1_and_2.col(0) = w1;
  w1_and_2.col(1) = w2;
  return r_inv * max(w1_and_2, 1) + (1 - r_inv) * min(w1_and_2, 1);
}

// Asymmetric logistic gauge: min(w1 + w2, logistic_gauge(w1, w2; dep))
// Interpolates between independence (w1+w2) and complete dependence.
// dep in (0, 1].
// [[Rcpp::export]]
arma::vec asym_log_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep) {
  double r_inv = 1 / dep(0);
  arma::mat w1_and_2(w1.size(), 2);
  w1_and_2.col(0) = w1;
  w1_and_2.col(1) = w2;
  return min(w1 + w2, (r_inv * max(w1_and_2, 1) + (1 - r_inv) * min(w1_and_2, 1)));
}

// Dirichlet gauge: g(w1, w2; theta1, theta2) = (1+theta1+theta2)*max(w1,w2) - theta1*w1 - theta2*w2
// Two-parameter family; theta1, theta2 > 0 with half-t(4) priors.
// dep = c(theta1, theta2).
// [[Rcpp::export]]
arma::vec dirichlet_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep) {
  double theta1 = dep(0), theta2 = dep(1);
  arma::mat w1_and_2(w1.size(), 2);
  w1_and_2.col(0) = w1;
  w1_and_2.col(1) = w2;
  return (1 + theta1 + theta2) * max(w1_and_2, 1) - (theta1 * w1 + theta2 * w2);
}

// -----------------------------------------------------------------------------
// String-to-function dispatcher
// -----------------------------------------------------------------------------

// Returns a function pointer to the gauge function corresponding to type_str.
// Throws an error via Rcpp::stop() for unrecognised strings.
gauge_function get_gauge_function(const std::string& type_str) {
  if      (type_str == "gauss")       return &gauss_gauge;
  else if (type_str == "inv_log")     return &inv_log_gauge;
  else if (type_str == "rectangular") return &rectangular_gauge;
  else if (type_str == "logistic")    return &logistic_gauge;
  else if (type_str == "asym_log")    return &asym_log_gauge;
  else if (type_str == "dirichlet")   return &dirichlet_gauge;
  else Rcpp::stop("Unknown gauge type: " + type_str +
    ". Valid options: gauss, inv_log, rectangular, logistic, asym_log, dirichlet.");
}
