// =============================================================================
// loglik.h
//
// Declarations for the scalar (summed) log-likelihood functions used inside
// the MCMC samplers. The angular_loglik function computes the total log-
// likelihood for the star-shaped angular density; radial_cens_loglik and
// radial_trunc_loglik implement the censored and truncated Gamma radial
// likelihoods respectively. get_loglik_function dispatches by string.
// =============================================================================
#ifndef LOGLIK_H
#define LOGLIK_H

#include <RcppArmadillo.h>
#include "gauge_functions.h"

// Type alias for a pointer to a radial log-likelihood function
typedef double (*radial_loglik_fn)(const arma::vec&, const arma::vec&, const arma::vec&,
                const double&, const arma::vec&, const std::string&);

// Angular log-likelihood (scalar sum over all observations)
double angular_loglik(const arma::vec& w1, const arma::vec& dep_par,
                      const arma::mat& grid_x, const int& dim,
                      const std::string& ang_gauge_type);

// Radial log-likelihood variants
double radial_cens_loglik(const arma::vec& r, const arma::vec& r0w, const arma::vec& w1,
                          const double& alpha, const arma::vec& dep_par,
                          const std::string& rad_gauge_type);
double radial_trunc_loglik(const arma::vec& r, const arma::vec& r0w, const arma::vec& w1,
                           const double& alpha, const arma::vec& dep_par,
                           const std::string& rad_gauge_type);

// Dispatcher: maps "cens" or "trunc" to the corresponding function pointer
radial_loglik_fn get_loglik_function(const std::string& type_str);

#endif  // LOGLIK_H
