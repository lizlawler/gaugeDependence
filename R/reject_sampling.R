#' Sample angles from the star-shaped angular density via rejection sampling
#'
#' Generates independent draws from the star-shaped angular density
#' f_W(w) proportional to g(w, 1-w; theta)^(-dim) using a uniform
#' envelope on (0, 1). Used for simulation and importance-sampling
#' prediction tasks.
#'
#' @param n Integer. Number of angular samples to generate.
#' @param posterior_params Numeric vector of dependence parameter(s) for the
#'   gauge function. For single-parameter gauges (all except Dirichlet),
#'   a scalar or length-1 vector. For the Dirichlet gauge, a length-2 vector
#'   c(theta1, theta2).
#' @param dim Integer. Dimension of the data (currently only dim = 2 supported).
#' @param gauge_type Character string specifying the gauge function family.
#'   One of \code{"gauss"}, \code{"logistic"}, \code{"inv_log"},
#'   \code{"asym_log"}, \code{"dirichlet"}, \code{"rectangular"}.
#'
#' @return Numeric vector of length \code{n} containing sampled angles
#'   in (0, 1).
#'
#' @examples
#' \dontrun{
#' # Sample 5000 angles from a Gaussian gauge star density with dep = 0.7
#' w_samples <- reject_sampling(n = 5000, posterior_params = 0.7,
#'                              dim = 2, gauge_type = "gauss")
#' hist(w_samples, breaks = 40, freq = FALSE,
#'      main = "Sampled angles", xlab = "W")
#' }
#'
#' @seealso \code{\link{angular_mcmc}}
#' @export
reject_sampling <- function(n,
                            posterior_params,
                            dim,
                            gauge_type) {
  # Build the 100 x 100 grid over [0,1]^2 for star-volume Monte Carlo estimation
  x1   <- seq(0, 1, length.out = 100)
  grid <- as.matrix(expand.grid(x1, x1))

  return(reject_samp(n,
                     as.numeric(posterior_params),
                     grid,
                     dim,
                     gauge_type))
}
