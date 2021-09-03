test_that("likelihood_shapes", {
  set.seed(45)
  n <- 2500
  d <- 4
  adj_test <- diag(d)
  adj_test[2, 1] <- 0.5
  adj_test[3, 1] <- 0.1

  beta_value <- 0.4999
  delta_time <- 0.1
  y_init <- rep(0, d)

  times <- seq(0, by = delta_time, length.out = n)
  noise <- matrix(rnorm(n * d, sd = sqrt(delta_time)), ncol = d)
  path <- construct_path(
    nw_topo = adj_test, noise = noise, y_init = y_init, delta_time = delta_time
  )
  thresholds <- rep(delta_time^beta_value, d)
  # without scaling
  loglik <- likelihood_fn(times, path, thresholds)
  loglik_val <- loglik(adj_test)
  testthat::expect_equal(length(loglik_val), 1)

  loglik <- likelihood_fn(times, path, thresholds, log = FALSE)
  loglik_val <- loglik(adj_test)
  testthat::expect_true(loglik_val < 0)

  # with scaling
  loglik <- likelihood_fn(times, path, thresholds, use_scaling = TRUE)
  loglik_val <- loglik(adj_test)
  testthat::expect_equal(length(loglik_val), 1)

  # small lambda
  loglik <- likelihood_fn(times, path, thresholds, lambda = 1e-4)
  loglik_val_small <- loglik(adj_test)

  # large lambda
  loglik <- likelihood_fn(times, path, thresholds, lambda = 1e4)
  loglik_val_large <- loglik(adj_test)
  testthat::expect_true(loglik_val_large > loglik_val_small)
})

test_that("grad_likelihood_shapes", {
  set.seed(45)
  n <- 2500
  d <- 4
  adj_test <- diag(d)
  adj_test[2, 1] <- 0.5
  adj_test[3, 1] <- 0.1

  beta_value <- 0.4999
  delta_time <- 0.1
  y_init <- rep(0, d)

  times <- seq(0, by = delta_time, length.out = n)
  noise <- matrix(rnorm(n * d, sd = sqrt(delta_time)), ncol = d)
  path <- construct_path(
    nw_topo = adj_test, noise = noise, y_init = y_init, delta_time = delta_time
  )
  thresholds <- rep(delta_time^beta_value, d)

  loglik <- likelihood_fn(times, path, thresholds)
  grad_loglik_numerical_val <- numDeriv::grad(function(x) -loglik(x), adj_test)
  # without scaling
  grad_loglik <- grad_likelihood_fn(times, path, thresholds)
  grad_loglik_val <- grad_loglik(adj_test)
  testthat::expect_equal(
    as.vector(grad_loglik_val), as.vector(grad_loglik_numerical_val),
    tolerance = 1e-3
  )

  # with scaling
  loglik <- likelihood_fn(times, path, thresholds, use_scaling = TRUE)
  grad_loglik_numerical_val <- numDeriv::grad(function(x) -loglik(x), adj_test)

  grad_loglik <- grad_likelihood_fn(times, path, thresholds, use_scaling = TRUE)
  grad_loglik_val <- grad_loglik(adj_test)
  testthat::expect_equal(
    as.vector(grad_loglik_val), as.vector(grad_loglik_numerical_val),
    tolerance = 1e-3
  )
})
