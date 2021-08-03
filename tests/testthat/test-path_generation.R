
test_that("construct_path__shape", {
  set.seed(42)

  corr_mat <- matrix(c(1.0, 0.2, 0.2, 0.2, 1.0, 0.2, 0.2, 0.2, 1.0), 3, 3)
  d <- dim(corr_mat)[1]
  n <- 1000
  delta_time <- 0.01
  noise <- matrix(rnorm(n * d), ncol = d)
  y_init <- rep(0, d)
  path <- construct_path(corr_mat, noise, y_init, delta_time)
  testthat::expect_equal(dim(path), c(n, d))
  lapply(
    seq(-10, 10, by = 5),
    function(x) {
      y_init <- rep(0, d)
      path <- construct_path(corr_mat, noise, y_init, delta_time)
      testthat::expect_equal(dim(path), c(n, d))
    }
  )
})

test_that("correlated_brownian_noise__shape", {
  set.seed(42)

  corr_mat <- matrix(c(1.0, 0.2, 0.2, 0.2, 1.0, 0.2, 0.2, 0.2, 1.0), 3, 3)
  d <- dim(corr_mat)[1]
  n <- 1000
  delta_time <- 0.01
  corr_inc <- correlated_brownian_noise(
    sigma_matrix = corr_mat, n = n, delta_time = delta_time
  )
  testthat::expect_equal(dim(corr_inc), c(n, d))
})

test_that("correlated_brownian_noise__corr", {
  set.seed(42)

  corr_mat <- matrix(c(1.0, 0.2, 0.2, 0.2, 1.0, 0.2, 0.2, 0.2, 1.0), 3, 3)
  d <- dim(corr_mat)[1]
  n <- 2000
  delta_time <- 0.01
  corr_inc <- correlated_brownian_noise(
    sigma_matrix = corr_mat, n = n, delta_time = delta_time
  )
  testthat::expect_equal(corr_mat, cor(corr_inc), tolerance = .5)
})

test_that("correlated_brownian_noise__shape", {
  set.seed(42)
  n <- 1000
  d <- 3
  delta_time <- 0.001
  n_jumps <- 1000
  jump_vals <- matrix(
    rep(rnorm(n_jumps, 0, sd = 0.1), d),
    nrow = n_jumps, ncol = d
  )
  cmpnd_poisson <- compound_poisson_jumps(
    d = d, n = n, delta_time = delta_time, jump_values = jump_vals
  )
  testthat::expect_equal(dim(cmpnd_poisson$noise), c(n, d))
  testthat::expect_equal(dim(cmpnd_poisson$jump_times), c(n, d))
})

test_that("correlated_brownian_noise__n_jumps_higher", {
  set.seed(42)
  n <- 1000
  d <- 3
  delta_time <- 0.001
  n_jumps <- 10000
  jump_vals <- matrix(
    rep(rnorm(n_jumps, 0, sd = 0.1), d),
    nrow = n_jumps, ncol = d
  )
  testthat::expect_warning(
    cmpnd_poisson <- compound_poisson_jumps(
      d = d, n = n, delta_time = delta_time, jump_values = jump_vals
    )
  )
  testthat::expect_equal(dim(cmpnd_poisson$noise), c(n, d))
  testthat::expect_equal(dim(cmpnd_poisson$jump_times), c(n, d))
})

test_that("correlated_brownian_noise__unidimensional", {
  set.seed(42)
  n <- 1000
  d <- 1
  delta_time <- 0.001
  n_jumps <- 500
  jump_vals <- rnorm(n_jumps, 0, sd = 0.1)

  cmpnd_poisson <- compound_poisson_jumps(
    d = d, n = n, delta_time = delta_time, jump_values = jump_vals
  )

  testthat::expect_equal(dim(cmpnd_poisson$noise), c(n, d))
  testthat::expect_equal(dim(cmpnd_poisson$jump_times), c(n, d))
})

test_that("addition_bm_and_compound_poisson__shape", {
  set.seed(42)
  n <- 1000
  d <- 3
  delta_time <- 0.001
  n_jumps <- 1000
  jump_vals <- matrix(
    rep(rnorm(n_jumps, 0, sd = 0.1), d),
    nrow = n_jumps, ncol = d
  )
  cmpnd_poisson <- compound_poisson_jumps(
    d = d, n = n, delta_time = delta_time, jump_values = jump_vals
  )
  cmpnd_poisson_noise <- apply(cmpnd_poisson$noise, 2, cumsum)
  sig_matrix <- matrix(-.2, 3, 3)
  diag(sig_matrix) <- 1
  corr_bw <- correlated_brownian_noise(
    sigma_matrix = sig_matrix, n = n, delta_time = delta_time
  )
  corr_bw <- apply(corr_bw, 2, cumsum)
  noise_add <- cmpnd_poisson_noise + corr_bw
  testthat::expect_equal(dim(noise_add), c(n, d))
})

test_that("correlated_jumps__shape", {
  set.seed(42)
  n <- 1000
  d <- 3
  delta_time <- 0.001
  n_jumps <- 1000
  corr_mat <- matrix(c(1.0, 0.2, 0.2, 0.2, 1.0, 0.2, 0.2, 0.2, 1.0), 3, 3)
  corr_jumps <- correlated_jumps(n, corr_mat)

  testthat::expect_equal(dim(corr_jumps), c(n, d))
})

test_that("bm_compound_poisson__shape", {
  set.seed(42)
  n <- 1000
  d <- 3
  delta_time <- 0.001
  n_jumps <- 1000
  corr_mat <- matrix(c(1.0, 0.2, 0.2, 0.2, 1.0, 0.2, 0.2, 0.2, 1.0), 3, 3)
  corr_jumps <- correlated_jumps(n, corr_mat)
  lapply(
    c(0, 1, 10, 100, 1000),
    function(n_j) {
      bm_cp <- bm_compound_poisson(
        n = n, corr_mat, diag(d), n_jumps = n_j, delta_time = delta_time
      )
      testthat::expect_equal(dim(bm_cp), c(n, d))
    }
  )
})
