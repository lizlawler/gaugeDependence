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
  
  // volume will stay consistent throughout, so eval once to be efficient
  gauge_function ang_dens_gauge_fn = get_gauge_function(ang_gauge_type);
  double L_volume = (ang_dens_gauge_fn == &logistic_gauge) ? dep_par(0) : est_star_vol(grid_x, dep_par, ang_gauge_type);
  double denom_dens = dim * L_volume;
  
  // envelope value; highest point of density
  double M = pow(ang_dens_gauge_fn({0.5}, {0.5}, dep_par)(0), -2) / denom_dens;
  Rcpp::Rcout << "M value is: " << M << std::endl;
  
  int i = 0;
  while(i < n) {
    double y = randu();
    double u = randu();
    double prop_y = R::dunif(y, 0, 1, false) * M;
    double f_y = pow(ang_dens_gauge_fn({y}, {1.0 - y}, dep_par)(0), -2) / denom_dens;
    double ratio = f_y / prop_y;
    if (u <= ratio) {
      x(i) = y;
      i += 1;
    }
  }
  return x;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
temp <- seq(0, 1, length.out = 100)
grid <- expand.grid(temp, temp)
test <- reject_samp(50000,
                    dep_par = 0.1,
                    grid_x = as.matrix(grid),
                    dim = 2,
                    ang_gauge_type = "asym_log")
*/
