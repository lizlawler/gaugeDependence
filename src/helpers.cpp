#include "helpers.h"
#include <RcppArmadillo.h>

double est_vol(arma::vec const& sum_term, arma::vec const& sqrt_term, double const& pars) {
  arma::vec gx = (sum_term - 2 * pars * sqrt_term) / (double)(1 - pow(pars, 2));
  return mean(arma::conv_to<arma::vec>::from(gx <= 1));
}

