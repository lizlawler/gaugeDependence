#' Compute posterior pointwise log-likelihood for the radial model
#'
#' For each row of \code{posterior_params} (one MCMC iteration), evaluates the
#' gauge-function radial log-likelihood at all observations. Returns an
#' (n_iter x n_obs) matrix of pointwise log-likelihoods for use with
#' \code{loo::loo()} or \code{loo::stacking_weights()}.
#'
#' @param radii Numeric vector of observed radii R. Should match the data
#'   passed to \code{radial_mcmc()}: exceedances only for \code{"trunc"},
#'   all observations for \code{"cens"}.
#' @param threshold Numeric vector of gauge threshold values r0(W), same
#'   length as \code{radii}.
#' @param angles Numeric vector of observed angles W in (0, 1), same length
#'   as \code{radii}.
#' @param posterior_params Numeric matrix of posterior MCMC samples
#'   (n_iter x n_params). Columns should correspond to alpha (shape), then the
#'   dependence parameter(s). The \code{sigma_m} column from
#'   \code{radial_mcmc()} output should be dropped before passing.
#' @param likelihood_type Character string. \code{"cens"} for the censored
#'   likelihood or \code{"trunc"} for the truncated likelihood. Must match the
#'   likelihood used to fit the model in \code{radial_mcmc()}.
#' @param gauge_type Character string specifying the gauge function family.
#'   One of \code{"gauss"}, \code{"logistic"}, \code{"inv_log"},
#'   \code{"asym_log"}, \code{"dirichlet"}, \code{"rectangular"}.
#'
#' @return Numeric matrix of dimension (n_iter x n_obs) containing the
#'   pointwise log-likelihood at each observation for each MCMC iteration.
#'
#' @examples
#' \dontrun{
#' data   <- RcppSimdJson::fload("data/gauss/high_1.json")
#' idx    <- data$idx
#' fit    <- radial_mcmc(radii = data$R[idx], threshold = data$r0_w[idx],
#'                       angles = data$W[idx], starting_theta = c(2, 0.5),
#'                       likelihood_type = "trunc", gauge_type = "gauss",
#'                       n_updates = 15000, n_burnin = 5000, n_thin = 5)
#' # Drop sigma_m before computing log-likelihoods
#' params    <- fit$samples[, c("alpha", "dep")]
#' pw_loglik <- radial_loglik(radii = data$R[idx], threshold = data$r0_w[idx],
#'                            angles = data$W[idx],
#'                            posterior_params = params,
#'                            likelihood_type = "trunc",
#'                            gauge_type = "gauss")
#' loo_result <- loo::loo(pw_loglik)
#' }
#'
#' @seealso \code{\link{radial_mcmc}}, \code{\link{angular_loglik}}
#' @export
radial_loglik <- function(radii,
                          threshold,
                          angles,
                          posterior_params,
                          likelihood_type,
                          gauge_type) {
  return(calc_loglik_radial(radii, threshold, angles,
                            as.matrix(posterior_params),
                            likelihood_type, gauge_type))
}
