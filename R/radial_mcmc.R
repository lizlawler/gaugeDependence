#' Run adaptive Metropolis-Hastings MCMC for radial dependence parameters
#'
#' Fits the gauge-function radial model by running an adaptive MH sampler for
#' the Gamma shape parameter (alpha) and gauge dependence parameter(s). Supports
#' both truncated (exceedances only) and censored (full-data) likelihoods.
#'
#' @param radii Numeric vector of observed radii R. For the truncated
#'   likelihood, pass only exceedances (R > r0(W)); for the censored likelihood,
#'   pass all observations.
#' @param threshold Numeric vector of gauge threshold values r0(W), evaluated
#'   at the corresponding angles. Must be the same length as \code{radii}.
#' @param angles Numeric vector of observed angles W in (0, 1). Must be the
#'   same length as \code{radii}.
#' @param starting_theta Numeric vector of initial parameter values.
#'   For all gauge types except \code{"dirichlet"}: length 2,
#'   c(alpha, dep). For \code{"dirichlet"}: length 3, c(alpha, theta1, theta2).
#'   If the wrong length is supplied, a warning is issued and values are drawn
#'   from the prior.
#' @param likelihood_type Character string. \code{"cens"} for the censored
#'   likelihood (all observations used, sub-threshold observations enter as
#'   CDF contributions) or \code{"trunc"} for the truncated likelihood
#'   (renormalized to the exceedance region).
#' @param gauge_type Character string specifying the gauge function family.
#'   One of \code{"gauss"}, \code{"logistic"}, \code{"inv_log"},
#'   \code{"asym_log"}, \code{"dirichlet"}, \code{"rectangular"}.
#' @param prop_Sigma_ Optional numeric matrix. Initial proposal covariance
#'   for the MH step. Defaults to an identity matrix if not provided or invalid.
#' @param n_updates Integer. Total number of MCMC iterations (default 10000).
#' @param update_freq Integer. Progress reporting frequency (default 100).
#' @param n_burnin Integer. Burn-in iterations to discard (default 0).
#' @param n_thin Integer. Thinning interval (default 1).
#' @param adapt_cov Logical. Adaptively update the proposal covariance
#'   (default FALSE).
#' @param r_opt Numeric. Target acceptance rate (default 0.234).
#' @param c0 Numeric. Adaptation rate constant (default 10).
#' @param c1 Numeric. Adaptation decay exponent (default 0.8).
#' @param K Integer. Iterations between adaptive updates (default 10).
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{samples}{Numeric matrix of post-burnin, thinned MCMC samples.
#'       Columns are named \code{"alpha"}, \code{"dep"}, and \code{"sigma_m"}
#'       for single-parameter gauges; or \code{"alpha"}, \code{"theta1"},
#'       \code{"theta2"}, and \code{"sigma_m"} for the Dirichlet gauge.}
#'     \item{acc_prob}{Numeric scalar. Overall sampler acceptance rate.}
#'   }
#'
#' @examples
#' \dontrun{
#' data <- RcppSimdJson::fload("data/gauss/high_1.json")
#' idx  <- data$idx  # threshold exceedance indices
#' fit  <- radial_mcmc(
#'   radii          = data$R[idx],
#'   threshold      = data$r0_w[idx],
#'   angles         = data$W[idx],
#'   starting_theta = c(rgamma(1, 4, 2), runif(1)),
#'   likelihood_type = "trunc",
#'   gauge_type      = "gauss",
#'   n_updates       = 15000,
#'   update_freq     = 250,
#'   n_burnin        = 5000,
#'   n_thin          = 5,
#'   adapt_cov       = TRUE
#' )
#' plot(fit$samples[, "dep"], type = "l", ylab = "dep", main = "Traceplot")
#' }
#'
#' @seealso \code{\link{angular_mcmc}}, \code{\link{radial_loglik}}
#' @export
radial_mcmc <- function(radii,
                        threshold,
                        angles,
                        starting_theta,
                        likelihood_type,
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

  # Input validation: radii, angles, and threshold must be the same length
  if ((length(radii) != length(angles)) || (length(radii) != length(threshold))) {
    stop("radii, angles, and threshold must all be the same length.")
  }

  # Validate and reinitialise starting values if wrong length is supplied
  expected_len <- if (gauge_type == "dirichlet") 3L else 2L
  if (length(starting_theta) != expected_len) {
    warning(sprintf(
      "starting_theta has length %d but should be %d for gauge_type = '%s'. Reinitialising from prior.",
      length(starting_theta), expected_len, gauge_type
    ))
    starting_theta <- if (gauge_type == "dirichlet") {
      c(rgamma(1, 4, 2), abs(rt(1, 4, ncp = 0)) * 4, abs(rt(1, 4, ncp = 0)) * 2)
    } else {
      c(rgamma(1, 4, 2), runif(1))
    }
  }

  return(radial_adaptive_mh(radii, threshold, angles,
                            starting_theta,
                            likelihood_type, gauge_type,
                            prop_Sigma_,
                            n_updates, update_freq, n_burnin, n_thin,
                            adapt_cov, r_opt, c0, c1, K))
}
