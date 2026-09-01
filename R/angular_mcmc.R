#' Run adaptive Metropolis-Hastings MCMC for angular dependence parameters
#'
#' Fits the star-shaped angular density model by running an adaptive MH sampler
#' for the gauge function dependence parameter(s). The angular log-likelihood is
#' evaluated on a 100 x 100 grid over [0,1]^2 used to estimate the star-body
#' volume via Monte Carlo.
#'
#' @param angles Numeric vector of observed angles W in (0, 1).
#' @param dim Integer. Dimension of the data (currently only dim = 2 is supported).
#' @param starting_theta Numeric vector of initial parameter values. Length 1
#'   for all gauge types except \code{"dirichlet"}, which requires length 2
#'   (theta1, theta2).
#' @param gauge_type Character string specifying the gauge function family.
#'   One of \code{"gauss"}, \code{"logistic"}, \code{"inv_log"},
#'   \code{"asym_log"}, \code{"dirichlet"}, \code{"rectangular"}.
#' @param prop_Sigma_ Optional numeric matrix. Initial proposal covariance
#'   matrix for the MH step. Defaults to an identity matrix if not provided or
#'   invalid.
#' @param n_updates Integer. Total number of MCMC iterations (default 10000).
#' @param update_freq Integer. Frequency at which progress is reported to the
#'   console (default 100).
#' @param n_burnin Integer. Number of iterations to discard as burn-in
#'   (default 0).
#' @param n_thin Integer. Thinning interval; every \code{n_thin}-th post-burnin
#'   sample is retained (default 1).
#' @param adapt_cov Logical. If TRUE, the proposal covariance matrix is
#'   updated adaptively during sampling (default FALSE).
#' @param r_opt Numeric. Target acceptance rate for the adaptive scaling step
#'   (default 0.234, the optimal rate for a single-parameter Gaussian target).
#' @param c0 Numeric. Adaptation rate constant controlling the speed of
#'   sigma_m adaptation (default 10).
#' @param c1 Numeric. Adaptation decay exponent in (0.5, 1]; controls how
#'   quickly the adaptation rate decreases (default 0.8).
#' @param K Integer. Number of iterations between each adaptive update step
#'   (default 10).
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{samples}{Numeric matrix of post-burnin, thinned MCMC samples.
#'       Columns are named \code{"dep"} and \code{"sigma_m"} for single-parameter
#'       gauges, or \code{"theta1"}, \code{"theta2"}, and \code{"sigma_m"} for
#'       the Dirichlet gauge. \code{sigma_m} is the adaptive proposal scale.}
#'     \item{acc_prob}{Numeric scalar. Overall acceptance rate of the sampler.}
#'   }
#'
#' @examples
#' \dontrun{
#' data <- RcppSimdJson::fload("data/gauss/high_1.json")
#' fit <- angular_mcmc(
#'   angles        = data$W,
#'   dim           = 2,
#'   starting_theta = 0.5,
#'   gauge_type    = "gauss",
#'   n_updates     = 15000,
#'   update_freq   = 250,
#'   n_burnin      = 5000,
#'   n_thin        = 5,
#'   adapt_cov     = TRUE
#' )
#' plot(fit$samples[, "dep"], type = "l", ylab = "dep", main = "Traceplot")
#' }
#'
#' @seealso \code{\link{radial_mcmc}}, \code{\link{angular_loglik}}
#' @export
angular_mcmc <- function(angles,
                         dim,
                         starting_theta,
                         gauge_type,
                         prop_Sigma_  = NULL,
                         n_updates    = 10000,
                         update_freq  = 100,
                         n_burnin     = 0,
                         n_thin       = 1,
                         adapt_cov    = FALSE,
                         r_opt        = 0.234,
                         c0           = 10,
                         c1           = 0.8,
                         K            = 10) {
  # Build the 100 x 100 grid over [0,1]^2 for star-volume Monte Carlo estimation
  x1   <- seq(0, 1, length.out = 100)
  grid <- as.matrix(expand.grid(x1, x1))

  return(angular_vol_adaptive_mh(angles,
                                 grid, dim,
                                 starting_theta,
                                 gauge_type,
                                 prop_Sigma_,
                                 n_updates, update_freq, n_burnin, n_thin,
                                 adapt_cov, r_opt, c0, c1, K))
}
