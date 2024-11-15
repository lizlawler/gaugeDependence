#ifndef LOGLIK_H
#define LOGLIK_H

#include <RcppArmadillo.h>


// Log likelihood function declarations
double angular_loglik(arma::vec const& W, double const& pars, 
                      arma::vec const& sum_term, arma::vec const& sqrt_term, int const& dim,
                      const std::string& ang_gauge_type);
double radial_cens_loglik(arma::vec const& R, arma::vec const& pars, 
                          arma::vec const& threshold, arma::vec const& W,
                          const std::string& rad_gauge_type);
double radial_top5_loglik(arma::vec const& R, arma::vec const& pars, 
                          arma::vec const& W,
                          const std::string& gauge_type);

#endif  // LOGLIK_H
