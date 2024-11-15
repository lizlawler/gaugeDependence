#' Run adaptive MH algorithm for radial-angular parameters
#' 

angular_mcmc <- function(angles, 
                         sum_term, sqrt_term, dim,
                         starting_theta_ang,
                         ang_gauge_type,
                         prop_Sigma_angles_ = NULL,
                         n_updates = 10000, update_freq = 100,
                         n_burnin = 0, n_thin = 1,
                         adapt_cov = FALSE,
                         r_opt = 0.234, 
                         c0 = 10, c1 = 0.8, 
                         K = 10) {
  if(length(sum_term) != length(sqrt_term)) {
    stop("Sum term and sqrt term must have the same number of values.")
  } else {
    return(adaptive_mh_ang_vol(angles, 
                               sum_term, sqrt_term, dim, 
                               starting_theta_ang,  
                               ang_gauge_type,
                               prop_Sigma_angles_, 
                               n_updates, update_freq, n_burnin, n_thin, 
                               adapt_cov, r_opt, c0, c1, K))
  }
}