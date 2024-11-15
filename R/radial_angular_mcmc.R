#' Run adaptive MH algorithm for radial-angular parameters
#' 

radial_angular_mcmc <- function(radii, angles, threshold, 
                                sum_term, sqrt_term, dim,
                                starting_theta_rad,
                                starting_theta_ang,
                                rad_gauge_type,
                                ang_gauge_type,
                                prop_Sigma_radii_ = NULL,
                                prop_Sigma_angles_ = NULL,
                                n_updates = 10000, update_freq = 100,
                                n_burnin = 0, n_thin = 1,
                                adapt_cov = FALSE,
                                r_opt = 0.234, 
                                c0 = 10, c1 = 0.8, 
                                K = 10) {
  if(length(radii) != length(angles) | length(angles) != length(threshold)) {
    stop("Radii, angles, and/or threshold values are not of the same length.")
  } else if(length(sum_term) != length(sqrt_term)) {
    stop("Sum term and sqrt term must have the same number of values.")
  } else {
    return(adaptive_mh_rad_ang_vol(radii, angles, threshold, sum_term, sqrt_term, dim,
                            starting_theta_rad, starting_theta_ang, rad_gauge_type, ang_gauge_type,
                            prop_Sigma_radii_, prop_Sigma_angles_, 
                            n_updates, update_freq, n_burnin, n_thin, 
                            adapt_cov, r_opt, c0, c1, K))
  }
}