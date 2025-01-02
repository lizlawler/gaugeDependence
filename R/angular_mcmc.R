#' Run adaptive MH algorithm for angular parameters
#' 

angular_mcmc <- function(angles,
                         dim,
                         starting_theta,
                         gauge_type,
                         prop_Sigma_ = NULL,
                         n_updates = 10000, update_freq = 100,
                         n_burnin = 0, n_thin = 1,
                         adapt_cov = FALSE,
                         r_opt = 0.234, 
                         c0 = 10, c1 = 0.8, 
                         K = 10) {
  x1 <- seq(0, 1, length.out = 100)
  grid <- as.matrix(expand.grid(x1, x1))
  return(angular_vol_adaptive_mh(angles, 
                                 grid, dim, 
                                 starting_theta,  
                                 gauge_type,
                                 prop_Sigma_, 
                                 n_updates, update_freq, n_burnin, n_thin, 
                                 adapt_cov, r_opt, c0, c1, K))
}