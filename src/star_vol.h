// =============================================================================
// star_vol.h
//
// Declaration for the Monte Carlo star-body volume estimator. The star body
// L = {x in R^d : g(x) <= 1} has volume estimated as the proportion of a
// uniform grid over [0,1]^2 for which the gauge function evaluates to <= 1.
// =============================================================================
#ifndef STAR_VOL_H
#define STAR_VOL_H

#include <RcppArmadillo.h>
#include "gauge_functions.h"

// Estimate the volume of the unit star body for a given gauge and parameters
double est_star_vol(const arma::mat& grid_x,
                    const arma::vec& pars,
                    const std::string& gauge_type);

#endif  // STAR_VOL_H
