// =============================================================================
// lpdf.cpp
//
// Observation-wise (pointwise) log-density functions for the angular and
// radial models. These return an arma::vec with one entry per observation,
// used in posterior_loglik_angular.cpp and posterior_loglik_radial.cpp to
// build the (n_iter x n_obs) log-likelihood matrix required by loo::loo().
//
// Mirrors loglik.cpp exactly in statistical content, differing only in that
// these functions return vectors rather than summing to a scalar.
// =============================================================================

#include "lpdf.h"
#include "star_vol.h"
#include "gauge_functions.h"
#include <RcppArmadillo.h>

// Angular pointwise log-density: log f_W(w_i) for each observation i
// f_W(w) = g(w, 1-w; theta)^{-dim} / (dim * vol(L))
arma::vec angular_lpdf(const arma::vec& w1, const arma::vec& dep_par,
                       const arma::mat& grid_x, const int& dim,
                       const std::string& ang_gauge_type) {
  gauge_function ang_dens_gauge_fn = get_gauge_function(ang_gauge_type);
  int N = w1.size();
  arma::vec pw_loglik(N);
  const arma::vec& w2 = 1 - w1;
  double L_volume = (ang_dens_gauge_fn == &logistic_gauge)
                      ? dep_par(0)
                      : est_star_vol(grid_x, dep_par, ang_gauge_type);
  arma::vec gauge_vec = ang_dens_gauge_fn(w1, w2, dep_par);
  for (int i = 0; i < N; ++i) {
    pw_loglik(i) = -(double)dim * log(gauge_vec(i)) - log((double)dim * L_volume);
  }
  return pw_loglik;
}

// Radial censored pointwise log-likelihood
// Sub-threshold observations contribute log CDF; exceedances contribute log density.
arma::vec radial_cens_lpdf(const arma::vec& r, const arma::vec& r0w, const arma::vec& w1,
                           const double& alpha, const arma::vec& dep_par,
                           const std::string& rad_gauge_type) {
  gauge_function gauge_fn = get_gauge_function(rad_gauge_type);
  const int& N = r.size();
  arma::vec pw_loglik(N);
  const arma::vec& w2 = 1 - w1;
  arma::vec beta = gauge_fn(w1, w2, dep_par);
  for (int i = 0; i < N; i++) {
    double scale_temp = 1 / beta(i);
    pw_loglik(i) = (r(i) < r0w(i))
                   ? R::pgamma(r0w(i), alpha, scale_temp, true, true)   // log CDF
                   : R::dgamma(r(i),   alpha, scale_temp, true);         // log density
  }
  return pw_loglik;
}

// Radial truncated pointwise log-density (renormalised to exceedance region)
arma::vec radial_trunc_lpdf(const arma::vec& r, const arma::vec& r0w, const arma::vec& w1,
                            const double& alpha, const arma::vec& dep_par,
                            const std::string& rad_gauge_type) {
  gauge_function gauge_fn = get_gauge_function(rad_gauge_type);
  const int& N = r.size();
  arma::vec pw_loglik(N);
  const arma::vec& w2 = 1 - w1;
  arma::vec beta = gauge_fn(w1, w2, dep_par);
  for (int i = 0; i < N; i++) {
    double scale_temp = 1 / beta(i);
    pw_loglik(i) = R::dgamma(r(i), alpha, scale_temp, true)
                 - R::pgamma(r0w(i), alpha, scale_temp, false, true);  // log CCDF
  }
  return pw_loglik;
}

// Dispatcher: maps likelihood type string to the corresponding function pointer
radial_lpdf_fn get_lpdf_function(const std::string& type_str) {
  if      (type_str == "cens")  return &radial_cens_lpdf;
  else if (type_str == "trunc") return &radial_trunc_lpdf;
  else Rcpp::stop("Unknown likelihood type: " + type_str + ". Valid options: cens, trunc.");
}
