#' Calculate pointwise loglikelihood using posterior parameters from angular volume MCMC
#' 

reject_sampling <- function(n,
                            posterior_params,
                            dim,
                            gauge_type) {
  x1 <- seq(0, 1, length.out = 100)
  grid <- as.matrix(expand.grid(x1, x1))
  return(reject_samp(n, 
                     as.matrix(posterior_params),
                     grid, 
                     dim,
                     gauge_type))
}