// =============================================================================
// lpdf.h
//
// Declarations for the observation-wise (pointwise) log-density functions
// used in the posterior log-likelihood computation for LOO-CV. Unlike the
// functions in loglik.h (which return a scalar sum), these return an arma::vec
// with one entry per observation. get_lpdf_function dispatches by string.
// =============================================================================
#ifndef LPDF_H
#define LPDF_H

#include <RcppArmadillo.h>
#include "gauge_functions.h"

// Type alias for a pointer to a pointwise radial log-density function
typedef arma::vec (*radial_lpdf_fn)(const arma::vec&, const arma::vec&, const arma::vec&,
                   const double&, const arma::vec&, const std::string&);

// Angular pointwise log-density (vector over observations)
arma::vec angular_lpdf(const arma::vec& w1, const arma::vec& dep_par,
                       const arma::mat& grid_x, const int& dim,
                       const std::string& ang_gauge_type);

// Radial pointwise log-density variants
arma::vec radial_cens_lpdf(const arma::vec& r, const arma::vec& r0w, const arma::vec& w1,
                           const double& alpha, const arma::vec& dep_par,
                           const std::string& rad_gauge_type);
arma::vec radial_trunc_lpdf(const arma::vec& r, const arma::vec& r0w, const arma::vec& w1,
                            const double& alpha, const arma::vec& dep_par,
                            const std::string& rad_gauge_type);

// Dispatcher: maps "cens" or "trunc" to the corresponding function pointer
radial_lpdf_fn get_lpdf_function(const std::string& type_str);

#endif  // LPDF_H
