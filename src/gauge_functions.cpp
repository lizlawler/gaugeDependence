#include "gauge_functions.h"
#include <RcppArmadillo.h>
#include <string>

// Gauge function definitions
// Asymptotic independent gauges

// [[Rcpp::export]]
arma::vec gauss_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep) {
  double dep_value = dep(0);
  return (w1 + w2 - 2 * dep_value * sqrt(w1 % w2)) / (1 - pow(dep_value, 2));
}

// [[Rcpp::export]]
arma::vec inv_log_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep) {
  double r_inv = 1 / dep(0);
  return pow((pow(w1, r_inv) + pow(w2, r_inv)), dep(0));
}

// [[Rcpp::export]]
arma::vec rectangular_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep) {
  double dep_value = dep(0);
  arma::mat w_max(w1.size(), 3);
  w_max.col(0) = (w1 - w2) / dep_value;
  w_max.col(1) = (w2 - w1) / dep_value;
  w_max.col(2) = ((w1 + w2) / (2 - dep_value));
  return max(w_max, 1);
}

// [[Rcpp::export]]
arma::vec logistic_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep) {
  double r_inv = 1 / dep(0);
  arma::mat w1_and_2(w1.size(), 2);
  w1_and_2.col(0) = w1;
  w1_and_2.col(1) = w2;
  return r_inv * max(w1_and_2, 1) + (1 - r_inv) * min(w1_and_2, 1);
}

// [[Rcpp::export]]
arma::vec asym_log_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep) {
  double r_inv = 1 / dep(0);
  arma::mat w1_and_2(w1.size(), 2);
  w1_and_2.col(0) = w1;
  w1_and_2.col(1) = w2;
  return min(w1 + w2, (r_inv * max(w1_and_2, 1) + (1 - r_inv) * min(w1_and_2, 1)));
}

// [[Rcpp::export]]
arma::vec dirichlet_gauge(const arma::vec& w1, const arma::vec& w2, const arma::vec& dep) {
  double theta1 = dep(0), theta2 = dep(1);
  arma::mat w1_and_2(w1.size(), 2);
  w1_and_2.col(0) = w1;
  w1_and_2.col(1) = w2;
  return (1 + theta1 + theta2) * max(w1_and_2, 1) - (theta1 * w1 + theta2 * w2);  
}

// Function to map strings to gauge functions
gauge_function get_gauge_function(const std::string& type_str) {
  if (type_str == "gauss") {
    return &gauss_gauge;
  } else if (type_str == "inv_log") {
    return &inv_log_gauge;
  } else if (type_str == "rectangular") {
    return &rectangular_gauge;
  } else if (type_str == "logistic") {
    return &logistic_gauge;
  } else if (type_str == "asym_log") {
    return &asym_log_gauge;
  } else if (type_str == "dirichlet") {
    return &dirichlet_gauge;
  } else {
    Rcpp::stop("Unknown gauge type: " + type_str);
  }
}

/*** R
# w <- seq(0, 1, length.out = 200)
# gw <- rectangular_gauge(w, 1-w, 0.5)
# plot(w/gw, (1-w)/gw)
*/