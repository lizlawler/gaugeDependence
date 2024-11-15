#ifndef HELPERS_H
#define HELPERS_H

#include <RcppArmadillo.h>


// Function declarations
double est_vol(arma::vec const& sum_term, arma::vec const& sqrt_term, double const& pars);

#endif  // HELPERS_H
