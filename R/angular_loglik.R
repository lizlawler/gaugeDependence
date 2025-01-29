#' Calculate pointwise loglikelihood using posterior parameters from angular volume MCMC
#' 

angular_loglik <- function(angles,
                           dim,
                           posterior_params,
                           gauge_type) {
  x1 <- seq(0, 1, length.out = 100)
  grid <- as.matrix(expand.grid(x1, x1))
  return(calc_loglik_angular(angles, 
                             grid, dim, 
                             as.matrix(posterior_params),  
                             gauge_type))
}