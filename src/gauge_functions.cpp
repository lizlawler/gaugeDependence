#include "gauge_functions.h"
#include <RcppArmadillo.h>
#include <string>

arma::vec gauss_gauge(arma::vec const& W, const double& dep) {
  arma::vec w2 = 1 - W;
  return (W + w2 - 2 * dep * sqrt(W % w2)) / (1 - pow(dep, 2));
}

arma::vec logistic_gauge(arma::vec const& W, const double& dep) {
  arma::mat w1_and_2(W.size(), 2);
  w1_and_2.col(0) = W;
  w1_and_2.col(1) = 1 - W;
  double r_inv = 1 / dep;
  return r_inv * max(w1_and_2, 1) + (1 - r_inv) * min(w1_and_2, 1);
}

// Function to map strings to gauge functions
gauge_function get_gauge_function(const std::string& type_str) {
  if (type_str == "gauss") {
    return &gauss_gauge;
  } else if (type_str == "logistic") {
    return &logistic_gauge;
  } else {
    Rcpp::stop("Unknown gauge type: " + type_str);
  }
}