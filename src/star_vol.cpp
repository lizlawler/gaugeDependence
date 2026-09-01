// =============================================================================
// star_vol.cpp
//
// Monte Carlo estimator for the volume of the unit star body
// L = {x in [0,1]^2 : g(x) <= 1}.
//
// The star body volume is needed to normalise the angular density:
//   f_W(w) = g(w, 1-w)^{-dim} / (dim * vol(L))
//
// Volume is estimated as the proportion of grid_x points for which
// g(x1, x2; pars) <= 1. The calling R wrappers use a 100 x 100 uniform
// grid over [0,1]^2, giving a Monte Carlo estimate with ~10,000 grid points.
//
// Note: for the logistic gauge, the star-body volume equals the dependence
// parameter dep directly, so the estimation step is skipped in loglik.cpp
// and lpdf.cpp as an efficiency optimisation.
// =============================================================================

#include "star_vol.h"
#include "gauge_functions.h"
#include <RcppArmadillo.h>

// [[Rcpp::export]]
double est_star_vol(const arma::mat& grid_x,
                    const arma::vec& pars,
                    const std::string& ang_gauge_type) {
  gauge_function ang_dens_gauge_fn = get_gauge_function(ang_gauge_type);
  // Evaluate gauge at each grid point and return proportion with g(x) <= 1
  arma::vec gx = ang_dens_gauge_fn(grid_x.col(0), grid_x.col(1), pars);
  return mean(arma::conv_to<arma::vec>::from(gx <= 1));
}
