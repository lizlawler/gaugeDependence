// =============================================================================
// gauge_functions.h
//
// Declarations for the six bivariate gauge function families and the
// string-to-function dispatcher used throughout the C++ sampler code.
// All gauge functions share the signature:
//   arma::vec gauge_fn(w1, w2, dep)
// where w1, w2 are angles in (0,1) with w1 + w2 = 1, and dep is the
// (possibly multi-element) dependence parameter vector.
// =============================================================================
#ifndef GAUGE_FUNCTIONS_H
#define GAUGE_FUNCTIONS_H

#include <RcppArmadillo.h>
#include <string>

// Type alias for a pointer to any gauge function
typedef arma::vec (*gauge_function)(const arma::vec&, const arma::vec&, const arma::vec&);

// Gauge function declarations
arma::vec gauss_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep);
arma::vec inv_log_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep);
arma::vec rectangular_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep);
arma::vec logistic_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep);
arma::vec asym_log_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep);
arma::vec dirichlet_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep);

// Dispatcher: maps a gauge type string to the corresponding function pointer
gauge_function get_gauge_function(const std::string& type_str);

#endif  // GAUGE_FUNCTIONS_H
