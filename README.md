# gaugeDependence

An R package providing Bayesian MCMC samplers for modeling bivariate extremes 
through a geometric framework. The core samplers and gauge functions 
are implemented in C++ via Rcpp and RcppArmadillo for performance.

## Background

The package implements a limit-set representation of the bivariate extremes framework. 
A bivariate observation (X1, X2) with standard exponential margins is decomposed into 
a radial component R = X1 + X2 and an angular component W = X1 / R. Above a threshold 
r0(W), the radial component is modelled as:

```
R | W ~ Gamma(2, g(W; theta))
```

where g(·; theta) is a parametric gauge function. The angular component W is 
modelled via a star-shaped density:

```
f_W(w) = g(w, 1-w; theta)^{-2} / (2 * vol(L))
```

where vol(L) is the volume of the unit star body {x : g(x) <= 1}.

## Installation

```r
# Install from GitHub (requires remotes)
remotes::install_github("lizlawler/gaugeDependence")
```

The package requires a C++ compiler and the RcppArmadillo headers, which are 
installed automatically with the package.

## Gauge families

Six parametric gauge functions are implemented:

| Name | `gauge_type` | Parameter | Dependence class |
|------|-------------|-----------|-----------------|
| Gaussian | `"gauss"` | rho in (0, 1) | Asymptotic independence (AI) |
| Logistic | `"logistic"` | dep in (0, 1] | Asymptotic dependence (AD) |
| Inverted logistic | `"inv_log"` | dep in (0, 1] | AI |
| Asymmetric logistic | `"asym_log"` | dep in (0, 1] | AD |
| Rectangular | `"rectangular"` | dep in (0, 1] | AI |
| Dirichlet | `"dirichlet"` | c(theta1, theta2) > 0 | AD |

All six gauge functions are exported to R and can be evaluated directly:

```r
w <- seq(0, 1, length.out = 200)
plot(w / gauss_gauge(w, 1-w, 0.5), (1-w) / gauss_gauge(w, 1-w, 0.5),
     type = "l", xlab = "x1", ylab = "x2",
     main = "Unit level set: Gaussian gauge, rho = 0.5")
```

## Main functions

### MCMC samplers

```r
# Fit the angular model (star-shaped density)
ang_fit <- angular_mcmc(
  angles        = w,             # observed angles in (0,1)
  dim           = 2,
  starting_theta = 0.5,          # initial dep value
  gauge_type    = "gauss",
  n_updates     = 15000,
  update_freq   = 250,
  n_burnin      = 5000,
  n_thin        = 5,
  adapt_cov     = TRUE
)

# Fit the radial model
rad_fit <- radial_mcmc(
  radii          = r[idx],       # threshold exceedances
  threshold      = r0w[idx],     # gauge threshold at each angle
  angles         = w[idx],
  starting_theta = c(rgamma(1, 4, 2), runif(1)),  # c(alpha, dep)
  likelihood_type = "trunc",     # "trunc" or "cens"
  gauge_type     = "gauss",
  n_updates      = 15000,
  update_freq    = 250,
  n_burnin       = 5000,
  n_thin         = 5,
  adapt_cov      = TRUE
)
```

Both functions return a list with:
- `$samples`: matrix of post-burnin thinned MCMC samples
- `$acc_prob`: overall acceptance rate

### Posterior log-likelihoods (for model comparison via `loo`)

```r
# Angular pointwise log-likelihood (n_iter x n_obs matrix)
ang_pw <- angular_loglik(
  angles           = w,
  dim              = 2,
  posterior_params = ang_fit$samples[, "dep", drop = FALSE],
  gauge_type       = "gauss"
)

# Radial pointwise log-likelihood
rad_pw <- radial_loglik(
  radii            = r[idx],
  threshold        = r0w[idx],
  angles           = w[idx],
  posterior_params = rad_fit$samples[, c("alpha", "dep")],
  likelihood_type  = "trunc",
  gauge_type       = "gauss"
)

# LOO-CV and model comparison
library(loo)
joint_pw   <- rad_pw + ang_pw  # joint pointwise log-likelihood
loo_result <- loo(joint_pw)
```

### Rejection sampling from the angular density

```r
# Sample angles from the posterior predictive angular density
w_sim <- reject_sampling(
  n               = 5000,
  posterior_params = median(ang_fit$samples[, "dep"]),
  dim             = 2,
  gauge_type      = "gauss"
)
```

## Likelihood options

The radial model supports two likelihoods:

- **`"cens"` (censored):** uses all observations; sub-threshold values enter as 
  CDF contributions, exceedances as density contributions.
- **`"trunc"` (truncated):** uses exceedances only; the density is renormalised 
  to the exceedance region P(R > r0(W) | W).

## Citation

If you use this package, please cite the associated manuscript:

> E. Lawler and B. Shaby (2026). *Bayesian model averaging of risk set probabilities
> using a geometric representation of multivariate extremes*.
> arXiv:2608.28888. https://arxiv.org/abs/2608.28888

The manuscript is currently under review; this citation will be updated
if and when it is published.

## License

GPL (>= 2)
