#ifndef PRIOR_H
#define PRIOR_H

#include <RcppArmadillo.h>


// Prior function declarations
double prior_fn_radii(const arma::vec& params);
double prior_fn_angles(const arma::vec& params);


#endif  // PRIOR_H
