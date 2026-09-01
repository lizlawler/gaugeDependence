## gaugeDependence-package.R
## Package-level documentation and namespace directives.
## roxygen2 uses the tags below to generate the NAMESPACE file, so NAMESPACE
## should NOT be hand-edited -- run devtools::document() to regenerate it.

#' gaugeDependence: Bayesian MCMC for Bivariate Extremes via Gauge Functions
#'
#' Adaptive Metropolis-Hastings MCMC samplers for fitting bivariate extreme
#' value models using gauge function dependence structures. The joint
#' distribution of a bivariate extreme is decomposed into radial (R) and
#' angular (W) components: radii above a threshold are modelled as
#' Gamma(2, g(W)), where g is a gauge function; angles are modelled via a
#' star-shaped angular density. Six parametric gauge families are implemented,
#' with truncated and censored likelihood options for the radial component.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{radial_mcmc}}}{Fit the radial model via adaptive MH MCMC.}
#'   \item{\code{\link{angular_mcmc}}}{Fit the star-shaped angular model via adaptive MH MCMC.}
#'   \item{\code{\link{radial_loglik}}}{Posterior pointwise log-likelihood for the radial model.}
#'   \item{\code{\link{angular_loglik}}}{Posterior pointwise log-likelihood for the angular model.}
#'   \item{\code{\link{reject_sampling}}}{Sample angles from the star-shaped density.}
#' }
#'
#' @section Gauge functions:
#' The six gauge function families are also exported and can be evaluated
#' directly: \code{\link{gauss_gauge}}, \code{\link{logistic_gauge}},
#' \code{\link{inv_log_gauge}}, \code{\link{asym_log_gauge}},
#' \code{\link{dirichlet_gauge}}, and \code{\link{rectangular_gauge}}.
#'
#' @keywords internal
#' @aliases gaugeDependence-package
#' @useDynLib gaugeDependence, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats rgamma rt runif
"_PACKAGE"

## ---------------------------------------------------------------------------
## Explicit exports for the Rcpp-defined gauge functions.
## Their R bindings live in the auto-generated RcppExports.R (which cannot
## carry roxygen tags), so the @export directives are declared here instead.
## Only the gauge functions are exported; the low-level MCMC / likelihood
## internals (radial_adaptive_mh, angular_vol_adaptive_mh, calc_loglik_*,
## reject_samp, est_star_vol) remain unexported and are reached through the
## user-facing R wrappers.
## ---------------------------------------------------------------------------

#' Gaussian gauge function
#'
#' Evaluates the Gaussian (asymptotic independence) gauge function, which
#' corresponds to a Gaussian copula with correlation \code{dep}.
#'
#' @param w1,w2 Numeric vectors of angular coordinates (with \code{w1 + w2 = 1}).
#' @param dep Numeric dependence parameter (correlation) in (0, 1).
#' @return Numeric vector of gauge function values.
#' @name gauss_gauge
#' @export
NULL

#' Logistic (Gumbel) gauge function
#'
#' @param w1,w2 Numeric vectors of angular coordinates (with \code{w1 + w2 = 1}).
#' @param dep Numeric dependence parameter in (0, 1].
#' @return Numeric vector of gauge function values.
#' @name logistic_gauge
#' @export
NULL

#' Inverted logistic gauge function
#'
#' @param w1,w2 Numeric vectors of angular coordinates (with \code{w1 + w2 = 1}).
#' @param dep Numeric dependence parameter in (0, 1].
#' @return Numeric vector of gauge function values.
#' @name inv_log_gauge
#' @export
NULL

#' Asymmetric logistic gauge function
#'
#' @param w1,w2 Numeric vectors of angular coordinates (with \code{w1 + w2 = 1}).
#' @param dep Numeric dependence parameter in (0, 1].
#' @return Numeric vector of gauge function values.
#' @name asym_log_gauge
#' @export
NULL

#' Dirichlet gauge function
#'
#' @param w1,w2 Numeric vectors of angular coordinates (with \code{w1 + w2 = 1}).
#' @param dep Numeric vector of two positive dependence parameters,
#'   \code{c(theta1, theta2)}.
#' @return Numeric vector of gauge function values.
#' @name dirichlet_gauge
#' @export
NULL

#' Rectangular gauge function
#'
#' @param w1,w2 Numeric vectors of angular coordinates (with \code{w1 + w2 = 1}).
#' @param dep Numeric dependence parameter in (0, 1].
#' @return Numeric vector of gauge function values.
#' @name rectangular_gauge
#' @export
NULL
