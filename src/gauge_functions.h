#ifndef GAUGE_FUNCTIONS_H
#define GAUGE_FUNCTIONS_H

#include <RcppArmadillo.h>
#include <string> 

// Define the type for a gauge function pointer
typedef arma::vec (*gauge_function)(const arma::vec&, const double&);

// Function declarations
arma::vec gauss_gauge(const arma::vec& W, const double& dep);
arma::vec logistic_gauge(const arma::vec& W, const double& dep);

// Function to map strings to gauge functions
gauge_function get_gauge_function(const std::string& type_str);

#endif  // GAUGE_FUNCTIONS_H
