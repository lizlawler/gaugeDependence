#ifndef PRIOR_H
#define PRIOR_H

#include <RcppArmadillo.h>


// Prior function declarations
double prior_fn_angles(double const& params);
double prior_fn_radii(arma::vec const& params);

#endif  // PRIOR_H
