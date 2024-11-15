#include "loglik.h"
#include "helpers.h"
#include "gauge_functions.h"
#include <RcppArmadillo.h>

// Angular loglikelihood
double angular_loglik(arma::vec const& W, double const& pars, 
                      arma::vec const& sum_term, arma::vec const& sqrt_term, int const& dim,
                      const std::string& ang_gauge_type) {
  gauge_function ang_dens_fn = get_gauge_function(ang_gauge_type);
  int N = W.size();
  double L_volume = (ang_dens_fn == &logistic_gauge) ? pars : est_vol(sum_term, sqrt_term, pars);
  return(-(double)dim * sum(log(ang_dens_fn(W, pars))) - (double)N * (log((double)dim) + log(L_volume)));
}

// Radial censored loglikelihood
double radial_cens_loglik(arma::vec const& R, arma::vec const& pars, 
                          arma::vec const& threshold, arma::vec const& W,
                          const std::string& rad_gauge_type) {
  gauge_function rad_dens_fn = get_gauge_function(rad_gauge_type);
  int n = R.size();
  double alpha = pars(0);
  double dep = pars(1);
  arma::vec beta = rad_dens_fn(W, dep);
  double loglik = 0.0;
  for(int i = 0; i < n; i ++) {
    if(R(i) < threshold(i)) {
      loglik += R::pgamma(threshold(i), alpha, 1/beta(i), true, true);
    } else {
      loglik += R::dgamma(R(i), alpha, 1/beta(i), true);
    }
  }
  return(loglik);
} 

// Radial top 5% loglikelihood; excluding normalizing constant
double radial_top5_loglik(arma::vec const& R, arma::vec const& pars, 
                          arma::vec const& W,
                          const std::string& gauge_type) {
  gauge_function dens_fn = get_gauge_function(gauge_type);
  int n = R.size();
  double alpha = pars(0);
  double dep = pars(1);
  arma::vec beta = dens_fn(W, dep);
  double loglik = 0.0;
  for(int i = 0; i < n; i ++) {
      loglik += R::dgamma(R(i), alpha, 1/beta(i), true);
  }
  return(loglik);
}
