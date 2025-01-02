#ifndef GAUGE_FUNCTIONS_H
#define GAUGE_FUNCTIONS_H

#include <RcppArmadillo.h>
#include <string>

// Define the type for a gauge function pointer
typedef arma::vec (*gauge_function)(const arma::vec&, const arma::vec&, const arma::vec&);

// Function declarations
arma::vec gauss_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep);
arma::vec inv_log_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep);
arma::vec rectangular_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep);
arma::vec logistic_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep);
arma::vec asym_log_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep);
arma::vec dirichlet_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep);

// Function to map strings to gauge functions
gauge_function get_gauge_function(const std::string& type_str);

#endif  // GAUGE_FUNCTIONS_H
