#' Run adaptive MH algorithm for radial parameters
#' 

radial_mcmc <- function(radii, threshold, angles,
                        starting_theta,
                        likelihood_type,
                        gauge_type,
                        prop_Sigma_ = NULL,
                        n_updates = 10000, update_freq = 100,
                        n_burnin = 0, n_thin = 1,
                        adapt_cov = FALSE,
                        r_opt = 0.234, 
                        c0 = 10, c1 = 0.8, 
                        K = 10) {
  if((length(radii) != length(angles)) || (length(radii) != length(threshold))) {
    stop("Radii, angles, and/or threshold values are not of the same length.")
  } else if(gauge_type != "dirichlet") {
    if(starting_theta != 2) {
      warning(paste0("Length of initial parameter values is ", length(starting_theta), ", but should be 2. Initializing from priors."))
      starting_theta <- c(rgamma(1, 4, 2), runif(1))
    }
  } else if(gauge_type == "dirichlet") {
    if(starting_theta != 3) {
      warning(paste0("Length of initial parameter values is ", length(starting_theta), ", but should be 3. Initializing from priors."))
      starting_theta <- c(rgamma(1, 4, 2), abs(rt(1, 4, ncp = 0))*4, abs(rt(1, 4, ncp = 0))*2)
    }
  } else {
    return(radial_adaptive_mh(radii, threshold, angles,
                              starting_theta, 
                              likelihood_type, gauge_type,
                              prop_Sigma_, 
                              n_updates, update_freq, n_burnin, n_thin, 
                              adapt_cov, r_opt, c0, c1, K))
  }
}