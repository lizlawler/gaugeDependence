#ifndef LOGLIK_H
#define LOGLIK_H

#include <RcppArmadillo.h>
#include "gauge_functions.h"

// Define the type for a gauge function pointer
typedef double (*radial_loglik_fn)(const arma::vec&, const arma::vec&, const arma::vec&,
                const double&, const arma::vec&,
                const std::string&);

// Log likelihood function declarations
double angular_loglik(const arma::vec& w1, const arma::vec& dep_par, 
                      const arma::mat& grid_x, const int& dim,
                      const std::string& ang_gauge_type);
double radial_cens_loglik(const arma::vec& r, const arma::vec& r0w, const arma::vec& w1,
                          const double& alpha, const arma::vec& dep_par,
                          const std::string& rad_gauge_type);
double radial_trunc_loglik(const arma::vec& r, const arma::vec& r0w, const arma::vec& w1,
                           const double& alpha, const arma::vec& dep_par,
                           const std::string& rad_gauge_type);

// Function to map strings to radial loglikelihood functions
radial_loglik_fn get_loglik_function(const std::string& type_str);

#endif  // LOGLIK_H