#' Run adaptive MH algorithm for radial-angular parameters
#' 

radial_top5_mcmc <- function(radii, angles,
                             starting_theta,
                             gauge_type,
                             prop_Sigma_ = NULL,
                             n_updates = 10000, update_freq = 100,
                             n_burnin = 0, n_thin = 1,
                             adapt_cov = FALSE,
                             r_opt = 0.234, 
                             c0 = 10, c1 = 0.8, 
                             K = 10) {
  if(length(radii) != length(angles)) {
    stop("Radii, angles, and/or threshold values are not of the same length.")
  } else {
    return(adaptive_mh(radii, angles,
                       starting_theta, gauge_type,
                       prop_Sigma_, 
                       n_updates, update_freq, n_burnin, n_thin, 
                       adapt_cov, r_opt, c0, c1, K))
  }
}