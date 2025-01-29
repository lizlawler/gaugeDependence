#ifndef LPDF_H
#define LPDF_H

#include <RcppArmadillo.h>
#include "gauge_functions.h"

// Define the type for a gauge function pointer
typedef arma::vec (*radial_lpdf_fn)(const arma::vec&, const arma::vec&, const arma::vec&,
                   const double&, const arma::vec&,
                   const std::string&);

// Log likelihood function declarations
arma::vec angular_lpdf(const arma::vec& w1, const arma::vec& dep_par, 
                       const arma::mat& grid_x, const int& dim,
                       const std::string& ang_gauge_type);
arma::vec radial_cens_lpdf(const arma::vec& r, const arma::vec& r0w, const arma::vec& w1,
                           const double& alpha, const arma::vec& dep_par,
                           const std::string& rad_gauge_type);
arma::vec radial_trunc_lpdf(const arma::vec& r, const arma::vec& r0w, const arma::vec& w1,
                            const double& alpha, const arma::vec& dep_par,
                            const std::string& rad_gauge_type);

// Function to map strings to radial loglikelihood functions
radial_lpdf_fn get_lpdf_function(const std::string& type_str);

#endif  // LOGLIK_H