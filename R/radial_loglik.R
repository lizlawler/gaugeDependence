#' Calculate pointwise loglikelihood using posterior parameters from radial MCMC
#' 

radial_loglik <- function(radii, threshold, angles,
                        posterior_params,
                        likelihood_type,
                        gauge_type) {
  
  return(calc_loglik_radial(radii, threshold, angles,
                            posterior_params, 
                            likelihood_type, gauge_type))
}