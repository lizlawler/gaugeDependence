#ifndef STAR_VOL_H
#define STAR_VOL_H

#include <RcppArmadillo.h>
#include "gauge_functions.h"

// Function declarations
double est_star_vol(const arma::mat& grid_x,
                    const arma::vec& pars, 
                    const std::string& gauge_type);

#endif  // STAR_VOL_H
