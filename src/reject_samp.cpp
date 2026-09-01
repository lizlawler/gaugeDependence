// =============================================================================
// reject_samp.cpp
//
// Rejection sampler for the star-shaped angular density
// f_W(w) = g(w, 1-w; theta)^{-dim} / (dim * vol(L))
//
// Uses a Uniform(0,1) proposal distribution with an envelope constant M
// chosen as the density evaluated at the mode (w = 0.5 for symmetric gauges).
// Samples are drawn one at a time until n accepted samples are obtained.
//
// Exported to R as reject_samp(), called by the R wrapper reject_sampling().
// =============================================================================

#include "star_vol.h"
#include "gauge_functions.h"
#include <RcppArmadillo.h>
#include <string>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
arma::vec reject_samp(const int& n,
                      const arma::vec& dep_par,
                      const arma::mat& grid_x,
                      const int& dim,
                      const std::string& ang_gauge_type) {
  vec x(n, fill::none);

  gauge_function ang_dens_gauge_fn = get_gauge_function(ang_gauge_type);

  // Compute star-body volume once (constant across all proposals)
  double L_volume = (ang_dens_gauge_fn == &logistic_gauge)
                      ? dep_par(0)
                      : est_star_vol(grid_x, dep_par, ang_gauge_type);
  double denom_dens = dim * L_volume;

  // Envelope constant M: supremum of f_W, evaluated at w = 0.5
  // (exact for symmetric gauges; conservative bound for asymmetric ones)
  double M = pow(ang_dens_gauge_fn({0.5}, {0.5}, dep_par)(0), -2) / denom_dens;
  Rcpp::Rcout << "Rejection sampling envelope M = " << M << std::endl;

  int i = 0;
  while (i < n) {
    double y = randu();          // proposal from Uniform(0, 1)
    double u = randu();          // uniform for accept/reject
    double prop_y = M;           // envelope density (Uniform scaled by M)
    double f_y = pow(ang_dens_gauge_fn({y}, {1.0 - y}, dep_par)(0), -2) / denom_dens;
    if (u <= f_y / prop_y) {
      x(i) = y;
      i++;
    }
  }
  return x;
}
