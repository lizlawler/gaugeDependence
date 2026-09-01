// =============================================================================
// loglik.cpp
//
// Scalar (summed) log-likelihood functions used inside the MCMC samplers.
// These functions compute the total log-likelihood across all observations
// for a given parameter vector, used in the Metropolis-Hastings accept/reject
// step. For observation-wise pointwise log-likelihoods (needed for LOO-CV),
// see lpdf.cpp.
//
// Functions:
//   angular_loglik      - total log-likelihood for the star-shaped angular density
//   radial_cens_loglik  - total log-likelihood, censored Gamma radial model
//   radial_trunc_loglik - total log-likelihood, truncated Gamma radial model
//   get_loglik_function - dispatcher mapping "cens"/"trunc" to the above
// =============================================================================

#include "loglik.h"
#include "star_vol.h"
#include "gauge_functions.h"
#include <RcppArmadillo.h>

// Angular log-likelihood: sum_i [ -dim * log(g(w_i)) - log(dim * vol(L)) ]
// The star-body volume is estimated via Monte Carlo using grid_x (100x100 grid).
// Optimisation: for the logistic gauge, vol(L) = dep directly (exact).
double angular_loglik(const arma::vec& w1, const arma::vec& dep_par,
                      const arma::mat& grid_x, const int& dim,
                      const std::string& ang_gauge_type) {
  gauge_function ang_dens_gauge_fn = get_gauge_function(ang_gauge_type);
  int N = w1.size();
  const arma::vec& w2 = 1 - w1;
  double L_volume = (ang_dens_gauge_fn == &logistic_gauge)
                      ? dep_par(0)
                      : est_star_vol(grid_x, dep_par, ang_gauge_type);
  double total_loglik = -(double)dim * sum(log(ang_dens_gauge_fn(w1, w2, dep_par)))
                        - (double)N * (log((double)dim) + log(L_volume));
  // Guard against Inf (e.g. from log(0) when gauge evaluates to 0)
  return (total_loglik == arma::datum::inf) ? -arma::datum::inf : total_loglik;
}

// Radial censored log-likelihood
// Sub-threshold observations (r < r0(w)) contribute the Gamma CDF;
// exceedances (r >= r0(w)) contribute the Gamma log-density.
double radial_cens_loglik(const arma::vec& r, const arma::vec& r0w, const arma::vec& w1,
                          const double& alpha, const arma::vec& dep_par,
                          const std::string& rad_gauge_type) {
  gauge_function gauge_fn = get_gauge_function(rad_gauge_type);
  const int& N = r.size();
  const arma::vec& w2 = 1 - w1;
  arma::vec beta = gauge_fn(w1, w2, dep_par);  // rate = g(w; theta)
  double loglik = 0.0;
  for (int i = 0; i < N; i++) {
    double scale_temp = 1 / beta(i);  // Gamma parameterised by scale = 1/rate
    if (r(i) < r0w(i)) {
      loglik += R::pgamma(r0w(i), alpha, scale_temp, true, true);   // log CDF
    } else {
      loglik += R::dgamma(r(i),   alpha, scale_temp, true);          // log density
    }
  }
  return loglik;
}

// Radial truncated log-likelihood (renormalised to the exceedance region)
// All observations are assumed to be exceedances (r >= r0(w)).
// Denominator P(R > r0(w) | W=w) normalises the density to the tail region.
double radial_trunc_loglik(const arma::vec& r, const arma::vec& r0w, const arma::vec& w1,
                           const double& alpha, const arma::vec& dep_par,
                           const std::string& rad_gauge_type) {
  gauge_function gauge_fn = get_gauge_function(rad_gauge_type);
  const int& N = r.size();
  const arma::vec& w2 = 1 - w1;
  arma::vec beta = gauge_fn(w1, w2, dep_par);
  double loglik = 0.0;
  for (int i = 0; i < N; i++) {
    double scale_temp = 1 / beta(i);
    // log f(r | r > r0(w)) = log f(r) - log P(R > r0(w))
    loglik += R::dgamma(r(i), alpha, scale_temp, true)
            - R::pgamma(r0w(i), alpha, scale_temp, false, true);  // log CCDF
  }
  return loglik;
}

// Dispatcher: maps likelihood type string to the corresponding function pointer
radial_loglik_fn get_loglik_function(const std::string& type_str) {
  if      (type_str == "cens")  return &radial_cens_loglik;
  else if (type_str == "trunc") return &radial_trunc_loglik;
  else Rcpp::stop("Unknown likelihood type: " + type_str + ". Valid options: cens, trunc.");
}
