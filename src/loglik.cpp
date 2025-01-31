#include "loglik.h"
#include "star_vol.h"
#include "gauge_functions.h"
#include <RcppArmadillo.h>

// Angular likelihood
double angular_loglik(const arma::vec& w1, const arma::vec& dep_par, 
                      const arma::mat& grid_x, const int& dim,
                      const std::string& ang_gauge_type) {
  gauge_function ang_dens_gauge_fn = get_gauge_function(ang_gauge_type);
  int N = w1.size();
  const arma::vec& w2 = 1 - w1;
  double L_volume = (ang_dens_gauge_fn == &logistic_gauge) ? dep_par(0) : est_star_vol(grid_x, dep_par, ang_gauge_type);
  double total_loglik = -(double)dim * sum(log(ang_dens_gauge_fn(w1, w2, dep_par))) - (double)N * (log((double)dim) + log(L_volume));
  if(total_loglik == arma::datum::inf) {
    return -arma::datum::inf;
  } else {
    return total_loglik;
  }
}

// Radial censored likelihood
double radial_cens_loglik(const arma::vec& r, const arma::vec& r0w, const arma::vec& w1,
                          const double& alpha, const arma::vec& dep_par,
                          const std::string& rad_gauge_type) {
  gauge_function gauge_fn = get_gauge_function(rad_gauge_type);
  const int& N = r.size();
  const arma::vec& w2 = 1 - w1;
  arma::vec beta = gauge_fn(w1, w2, dep_par);
  double loglik = 0.0;
  for(int i = 0; i < N; i ++) {
    double scale_temp = 1 / beta(i);
    if(r(i) < r0w(i)) {
      loglik += R::pgamma(r0w(i), alpha, scale_temp, true, true);
    } else { 
      loglik += R::dgamma(r(i), alpha, scale_temp, true);
    } 
  }
  return(loglik);
}  

// Radial truncated likelihood, appropriately renormalized
double radial_trunc_loglik(const arma::vec& r, const arma::vec& r0w, const arma::vec& w1,
                           const double& alpha, const arma::vec& dep_par,
                           const std::string& rad_gauge_type) {
  gauge_function gauge_fn = get_gauge_function(rad_gauge_type);
  const int& N = r.size();
  const arma::vec& w2 = 1 - w1;
  arma::vec beta = gauge_fn(w1, w2, dep_par);
  double loglik = 0.0;
  for(int i = 0; i < N; i ++) {
    double scale_temp = 1 / beta(i);
    loglik += R::dgamma(r(i), alpha, scale_temp, true) - R::pgamma(r0w(i), alpha, scale_temp, false, true);
  }
  return(loglik);
}

// Function to map strings to loglikelihood functions
radial_loglik_fn get_loglik_function(const std::string& type_str) {
  if (type_str == "cens") {
    return &radial_cens_loglik;
  } else if (type_str == "trunc") {
    return &radial_trunc_loglik;
  } else {
    Rcpp::stop("Unknown likelihood: " + type_str);
  }
}
