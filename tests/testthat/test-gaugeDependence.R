## Tests for gaugeDependence gauge functions and samplers.
## Run via devtools::test() or during R CMD check.

test_that("gauge functions return known values", {
  # gauss_gauge(0.5, 0.5, 0.5) = (1 - 2*0.5*0.5) / (1 - 0.25) = 0.5 / 0.75 = 2/3
  expect_equal(as.numeric(gauss_gauge(0.5, 0.5, 0.5)), 2 / 3, tolerance = 1e-10)

  # inv_log_gauge(0.5, 0.5, 0.5) = (0.5^2 + 0.5^2)^0.5 = sqrt(0.5)
  expect_equal(as.numeric(inv_log_gauge(0.5, 0.5, 0.5)), sqrt(0.5), tolerance = 1e-10)

  # logistic_gauge at w1 = w2: max and min are equal, so dep drops out -> 0.5
  expect_equal(as.numeric(logistic_gauge(0.5, 0.5, 0.7)), 0.5, tolerance = 1e-10)
  expect_equal(as.numeric(logistic_gauge(0.5, 0.5, 0.2)), 0.5, tolerance = 1e-10)
})

test_that("gauge functions preserve input length", {
  w1 <- c(0.2, 0.4, 0.6, 0.8)
  w2 <- 1 - w1
  expect_equal(length(gauss_gauge(w1, w2, 0.5)), length(w1))
  expect_equal(length(logistic_gauge(w1, w2, 0.5)), length(w1))
  expect_equal(length(dirichlet_gauge(w1, w2, c(1, 1))), length(w1))
})

test_that("gauge functions are homogeneous of order 1", {
  # g(t*w1, t*w2) = t * g(w1, w2)
  expect_equal(as.numeric(gauss_gauge(c(0.2), c(0.8), 0.3)),
    as.numeric(2 * gauss_gauge(c(0.1), c(0.4), 0.3)),
    tolerance = 1e-10
  )
  expect_equal(as.numeric(logistic_gauge(c(0.6), c(0.4), 0.5)),
    as.numeric(3 * logistic_gauge(c(0.2), c(0.4 / 3), 0.5)),
    tolerance = 1e-10
  )
})

test_that("gauge functions are strictly positive on the simplex interior", {
  w1 <- c(0.2, 0.4, 0.6, 0.8)
  w2 <- 1 - w1
  expect_true(all(gauss_gauge(w1, w2, 0.5) > 0))
  expect_true(all(inv_log_gauge(w1, w2, 0.5) > 0))
  expect_true(all(asym_log_gauge(w1, w2, 0.5) > 0))
  expect_true(all(rectangular_gauge(w1, w2, 0.5) > 0))
  expect_true(all(dirichlet_gauge(w1, w2, c(1, 1)) > 0))
})

test_that("radial_mcmc returns the documented output structure", {
  set.seed(1)
  n <- 50
  w <- runif(n, 0.05, 0.95)
  r0w <- rep(1, n)
  r <- r0w + rexp(n) # all above threshold

  fit <- radial_mcmc(
    radii = r,
    threshold = r0w,
    angles = w,
    starting_theta = c(2, 0.5),
    likelihood_type = "trunc",
    gauge_type = "gauss",
    n_updates = 200,
    update_freq = 1000, # large value suppresses progress printing
    n_burnin = 50,
    n_thin = 1
  )

  expect_true(is.list(fit))
  expect_true(all(c("samples", "acc_prob") %in% names(fit)))
  expect_true(all(c("alpha", "dep", "sigma_m") %in% colnames(fit$samples)))
  expect_true(fit$acc_prob >= 0 && fit$acc_prob <= 1)
})

test_that("radial_mcmc validates its inputs", {
  set.seed(1)
  n <- 50
  w <- runif(n, 0.05, 0.95)
  r0w <- rep(1, n)
  r <- r0w + rexp(n)

  # Mismatched input lengths raise an error
  expect_error(
    radial_mcmc(
      radii = r[1:10], threshold = r0w, angles = w,
      starting_theta = c(2, 0.5),
      likelihood_type = "trunc", gauge_type = "gauss",
      n_updates = 10
    ),
    regexp = "same length"
  )

  # Wrong-length starting_theta triggers a warning (and reinitialises)
  expect_warning(
    radial_mcmc(
      radii = r, threshold = r0w, angles = w,
      starting_theta = c(2), # should be length 2
      likelihood_type = "trunc", gauge_type = "gauss",
      n_updates = 50, update_freq = 1000
    ),
    regexp = "starting_theta"
  )
})

test_that("angular_mcmc returns the documented output structure", {
  set.seed(2)
  n <- 50
  w <- runif(n, 0.05, 0.95)

  ang_fit <- angular_mcmc(
    angles         = w,
    dim            = 2,
    starting_theta = 0.5,
    gauge_type     = "gauss",
    n_updates      = 200,
    update_freq    = 1000,
    n_burnin       = 50,
    n_thin         = 1
  )

  expect_true(is.list(ang_fit))
  expect_true(all(c("samples", "acc_prob") %in% names(ang_fit)))
  expect_true(all(c("dep", "sigma_m") %in% colnames(ang_fit$samples)))
})

test_that("log-likelihood helpers return (n_iter x n_obs) matrices", {
  set.seed(3)
  n <- 50
  w <- runif(n, 0.05, 0.95)
  r0w <- rep(1, n)
  r <- r0w + rexp(n)

  fit <- radial_mcmc(
    radii = r, threshold = r0w, angles = w,
    starting_theta = c(2, 0.5),
    likelihood_type = "trunc", gauge_type = "gauss",
    n_updates = 200, update_freq = 1000, n_burnin = 50, n_thin = 1
  )
  ang_fit <- angular_mcmc(
    angles = w, dim = 2, starting_theta = 0.5, gauge_type = "gauss",
    n_updates = 200, update_freq = 1000, n_burnin = 50, n_thin = 1
  )

  ang_params <- ang_fit$samples[, "dep", drop = FALSE]
  ang_ll <- angular_loglik(
    angles = w, dim = 2,
    posterior_params = ang_params,
    gauge_type = "gauss"
  )
  expect_equal(nrow(ang_ll), nrow(ang_params))
  expect_equal(ncol(ang_ll), length(w))

  rad_params <- fit$samples[, c("alpha", "dep")]
  rad_ll <- radial_loglik(
    radii = r, threshold = r0w, angles = w,
    posterior_params = rad_params,
    likelihood_type = "trunc",
    gauge_type = "gauss"
  )
  expect_equal(nrow(rad_ll), nrow(rad_params))
  expect_equal(ncol(rad_ll), length(r))
})

test_that("reject_sampling returns angles of the right length and range", {
  set.seed(4)
  samp <- reject_sampling(
    n = 500, posterior_params = 0.5,
    dim = 2, gauge_type = "gauss"
  )
  expect_equal(length(samp), 500)
  expect_true(all(samp > 0 & samp < 1))
})
