// [[Rcpp::depends(BH)]]
#include <boost/math/distributions/beta.hpp>
#include <algorithm>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]
arma::uword rmultinom_cluster(const arma::vec& arma_probs) {
  unsigned int N = arma_probs.size();
  Rcpp::NumericVector probs = Rcpp::NumericVector(arma_probs.begin(), arma_probs.end());
  Rcpp::IntegerVector outcome(N);
  rmultinom(1, probs.begin(), N, outcome.begin());
  
  // Extract cluster number based on position of '1' from on 'outcome'
  arma::uword cluster = arma::index_max(Rcpp::as<arma::vec>(Rcpp::wrap(outcome)));
  return cluster;
}

// [[Rcpp::export]]
arma::vec beta_lpdf_vectorized(const arma::vec& x, double shape1, double shape2) {
  arma::vec result(x.n_elem);
  for(int i = 0; i < x.n_elem; ++i) {
    result(i) = R::dbeta(x(i), shape1, shape2, true);
  }  
  return result;
}  

// [[Rcpp::export]]
double beta_pdf(const double& x, double shape1, double shape2) {
  return R::dbeta(x, shape1, shape2, false);
}  

// [[Rcpp::export]]
double MH_loglik_fn(arma::vec const& W, double const& mu, double const& tau) { 
  // W is all the W_i associated with this specific cluster
  double alpha = mu * tau;
  double beta = (1-mu) * tau;
  arma::vec loglik = beta_lpdf_vectorized(W, alpha, beta);
  return(arma::accu(loglik));
}

// [[Rcpp::export]]
arma::vec stick_breaking_probs(const arma::vec& sticks) {
  int L = sticks.n_elem + 1;
  arma::vec probs(L);
  probs(0) = sticks(0);
  double prod_probs = 1 - probs(0);
  for(int i = 1; i < L - 1; ++i) {
    probs(i) = sticks(i) * prod_probs;
    prod_probs *= (1-sticks(i));
  }  
  probs(L - 1) = prod_probs;
  return(probs);
}  

// [[Rcpp::export]]
arma::uvec set_diff(const arma::uvec& full_set, const arma::uvec& subset) {
  arma::uvec diff;
  for (arma::uword i = 0; i < full_set.n_elem; ++i) {
    if (!arma::any(subset == full_set(i))) {
      diff.insert_rows(diff.n_elem, arma::uvec({full_set(i)})); 
    }  
  }  
  return diff;
}  

// [[Rcpp::export]]
Rcpp::List stick_angular_mcmc(const arma::vec& angles,
                              const int& L = 10,
                              const int& n_updates = 10000, const int& update_freq = 100,
                              const int& n_burnin = 0, const int& n_thin = 1,
                              const double& r_opt = 0.234, 
                              const double& c0 = 10, const double& c1 = 0.8, 
                              const int& K = 10) {
  
  // Initialize adaptive parameters
  int k = 3; // iteration offest
  vec sigma_m_mu(L, fill::value(pow(2.4, 2)));
  vec sigma_m_tau(L, fill::value(pow(2.4, 2)));
  
  // Initialize cluster variables
  double alpha = R::rgamma(1, 1);
  vec sticks(L - 1);
  vec probs = ones(L);
  probs /= accu(probs); // Give all clusters equal probability initially
  uvec clusters(angles.size());
  for(int i = 0; i < angles.size(); ++i) {
    clusters(i) = rmultinom_cluster(probs);
  }  
  
  // Initialization for MH steps
  vec mu = randu(L);
  double beta = randg(distr_param(1.0, 8.0));
  vec tau = 1 / randg(L, distr_param(2.0, 1 / beta));
  
  // Initialize trace objects
  mat trace_mu = mat(n_updates, L, fill::zeros), trace_tau = mat(n_updates, L, fill::zeros), trace_probs = mat(n_updates, L, fill::zeros), 
    jump_trace_mu = mat(n_updates, L, fill::zeros), jump_trace_tau = mat(n_updates, L, fill::zeros),
    sigma_m_trace_mu = mat(n_updates, L, fill::zeros), sigma_m_trace_tau = mat(n_updates, L, fill::zeros);
  vec trace_beta = vec(n_updates, fill::zeros), trace_alpha = vec(n_updates, fill::zeros);
  
  trace_mu.row(0) = mu.t(), trace_tau.row(0) = tau.t(), trace_probs.row(0) = probs.t(), 
    trace_beta(0) = beta, trace_alpha(0) = alpha,
    sigma_m_trace_mu.row(0) = sigma_m_mu.t(), sigma_m_trace_tau.row(0) = sigma_m_tau.t();
  
  // Create full list of clusters for use in set difference 
  uvec all_clusters = regspace<uvec>(0, L - 1);
  
  // Main loop of MCMC (MH within Gibbs)
  for(int iter = 1; iter < n_updates; ++iter) {
    uvec unique_clusters = unique(clusters);
    uvec diff_clusters = set_diff(all_clusters, unique_clusters);
    
    // Start with Metropolis-Hastings updates of mu and tau for the m unique clusters (may not be all L)
    for(uword j = 0; j < unique_clusters.n_elem; ++j) {
      uvec mask = (clusters == unique_clusters(j));
      double prior, prior_star, likelihood, likelihood_star, mu_star, tau_star;
      
      // update mu first
      prior = R::dunif(mu(unique_clusters(j)), 0, 1, true);
      mu_star = mu(unique_clusters(j)) + sigma_m_mu(unique_clusters(j)) * randn();
      prior_star = R::dunif(mu_star, 0, 1, true);
      if(prior_star != -datum::inf) {
        likelihood = MH_loglik_fn(angles.elem(mask), mu(unique_clusters(j)), tau(unique_clusters(j)));
        likelihood_star = MH_loglik_fn(angles.elem(mask), mu_star, tau(unique_clusters(j)));
        if ((prior_star + likelihood_star - prior - likelihood) > log(randu())) {
          mu(unique_clusters(j)) = mu_star;
          likelihood = likelihood_star; // use this for the "original" likelihood in the tau update
          jump_trace_mu(iter, unique_clusters(j)) = 1;
        } 
      } 
      
      // update tau next
      prior = R::dgamma(1/tau(unique_clusters(j)), 2, 1/beta, true) - 2 * log(tau(unique_clusters(j)));
      tau_star = tau(unique_clusters(j)) + sigma_m_tau(unique_clusters(j)) * randn();
      prior_star = R::dgamma(1/tau_star, 2, 1/beta, true) - 2 * log(tau_star);
      if(prior_star != -datum::inf) {
        likelihood_star = MH_loglik_fn(angles.elem(mask), mu(unique_clusters(j)), tau_star);
        if ((prior_star + likelihood_star - prior - likelihood) > log(randu())) {
          tau(unique_clusters(j)) = tau_star;
          jump_trace_tau(iter, unique_clusters(j)) = 1;
        }
      }
    }
    
    // Sample independently from priors of mu and tau for the clusters not represented above
    if(diff_clusters.n_elem >= 1) {
      for(uword j = 0; j < diff_clusters.n_elem; ++j) {
        mu(diff_clusters(j)) = randu();
        tau(diff_clusters(j)) = 1 / R::rgamma(2, 1 / beta);
      }    
    } 
    
    // Adaptation portion of MH update
    if(((iter + 1) % K) == 0) {
      double gamma1 = c0 / pow((((iter + 1) / K) + k), c1),
        gamma2 = 1 / pow((((iter + 1) / K) + k), c1);
      for(int j = 0; j < L; ++j) {
        sigma_m_mu(j) = exp(log(sigma_m_mu(j) + 
          gamma1 * (mean(jump_trace_mu.col(j).subvec((iter - K + 1), iter)) - r_opt)));
        sigma_m_tau(j) = exp(log(sigma_m_tau(j) + 
          gamma1 * (mean(jump_trace_tau.col(j).subvec((iter - K + 1), iter)) - r_opt)));
      }
      sigma_m_trace_mu.row(iter) = sigma_m_mu.t(), sigma_m_trace_tau.row(iter) = sigma_m_tau.t();
    } // End adaptation
    
    
    // Update cluster assignments
    vec updated_probs(L);
    for(uword i = 0; i < clusters.n_elem; ++i) {
      double w_i = angles(i);
      for(uword l = 0; l < L; ++l) {
        double a = mu(l) * tau(l);
        double b = (1 - mu(l)) * tau(l);
        updated_probs(l) = probs(l) * beta_pdf(w_i, a, b);
      } 
      updated_probs /= accu(updated_probs);
      clusters(i) = rmultinom_cluster(updated_probs);
    } 
    
    // Update cluster probabilities
    uvec cluster_counts = hist(clusters, all_clusters);
    for(int l = 0; l < L - 2; ++l) {
      double a = 1 + cluster_counts(l);
      double b = alpha + accu(cluster_counts.subvec(l + 1, L - 1));
      sticks(l) = R::rbeta(a, b);
    } 
    probs = stick_breaking_probs(sticks);
    
    // Update remaining hyperparameters
    alpha = R::rgamma(L + 1, 1 / (2 - accu(log(1.0 - sticks))));
    beta = R::rgamma(2 * L + 1, 1 / (1/8 + accu(1/tau)));
    
    // Update trace/sample objects
    trace_mu.row(iter) = mu.t(), trace_tau.row(iter) = tau.t(), trace_probs.row(iter) = probs.t(), trace_alpha(iter) = alpha, trace_beta(iter) = beta;
    
    // Report progress to console
    if ((iter + 1) % update_freq == 0) {
      int progress = (100 * (iter + 1)) / n_updates;
      Rcpp::Rcout << "\rProgress: |"
                  << std::string(progress / 2, '-')
                  << std::string(50 - (progress / 2), ' ')
                  << "| " << progress << "% completed." << std::flush;
      Rcpp::checkUserInterrupt();
    }
  }
  
  // Return list with the samples and accept ance probabilities
  return Rcpp::List::create(
    Rcpp::Named("mu") = trace_mu,
    Rcpp::Named("tau") = trace_tau,
    Rcpp::Named("probs") = trace_probs,
    Rcpp::Named("alpha") = trace_alpha,
    Rcpp::Named("beta") = trace_beta,
    Rcpp::Named("acc_prob_mu") = mean(jump_trace_mu, 0),
    Rcpp::Named("acc_prob_tau") = mean(jump_trace_tau, 0)
  ); 
}

/*** R
# data <- RcppSimdJson::fload("data/gauss/low_1.json")
w <- data$W
results <- stick_angular_mcmc(w, L = 6, n_updates = 10000, update_freq = 250)
mix_beta <- function(angles, weights, mus, taus) {
  alphas <- mus * taus
  betas <- (1-mus) * taus
  n <- length(weights)
  dens <- 0.0
  for(i in 1:n) {
    dens <- dens + weights[i] * dbeta(angles, alphas[i], betas[i])
  }
  return(dens)
}
# 
post_dens <- mix_beta(w, colMeans(results$probs[5000:10000,]), colMeans(results$mu[5000:10000,]), colMeans(results$tau[5000:10000,]))
hist(w, freq= FALSE, breaks = 50)
points(w, post_dens)
round(colMeans(results$probs[5000:15000,]), 3)
round(colMeans(results$mu[5000:15000,]) , 3)
round(colMeans(results$tau[5000:15000,]), 3)
*/
