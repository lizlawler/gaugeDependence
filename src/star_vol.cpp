#include "star_vol.h"
#include "gauge_functions.h"
#include <RcppArmadillo.h>

//[[Rcpp::export]]
double est_star_vol(const arma::mat& grid_x,
                    const arma::vec& pars, 
                    const std::string& ang_gauge_type) {
  gauge_function ang_dens_gauge_fn = get_gauge_function(ang_gauge_type);
  arma::vec gx = ang_dens_gauge_fn(grid_x.col(0), grid_x.col(1), pars);
  return mean(arma::conv_to<arma::vec>::from(gx <= 1));
}

/*** R
# x1 <- seq(0, 1, length.out = 200)
# grid <- as.matrix(expand.grid(x1, x1))
# est_star_vol(grid, 0.2, "rectangular")
*/