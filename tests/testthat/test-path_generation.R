
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

test_that("construct_path__errors", {
  set.seed(42)

  corr_mat <- matrix(c(1.0, 0.2, 0.2, 0.2, 1.0, 0.2, 0.2, 0.2, 1.0), 3, 3)
  d <- dim(corr_mat)[1]
  n <- 1000
  delta_t <- 0.01
  noise <- matrix(rnorm(n * d), ncol = d)
  y_init <- rep(0, d)
  path <- construct_path(corr_mat, noise, y_init, delta_t)
  testthat::expect_error(construct_path(corr_mat, noise, y_init, -1))
  wg_corr_mat <- cbind(corr_mat, rep(0, d))
  testthat::expect_error(construct_path(wg_corr_mat, noise, y_init, delta_t))
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

test_that("correlated_brownian_noise__errors", {
  set.seed(42)
  n <- 1000
  delta_time <- 0.01
  testthat::expect_error(correlated_brownian_noise(
    sigma_matrix = matrix(1, 1, 1), n = n, delta_time = delta_time
  ))
  testthat::expect_error(correlated_brownian_noise(
    sigma_matrix = 2, n = n, delta_time = delta_time
  ))
  testthat::expect_error(correlated_brownian_noise(
    sigma_matrix = NA, n = n, delta_time = delta_time
  ))
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

test_that("correlated_brownian_noise__corr_structure", {
  set.seed(42)
  n <- 10000
  d <- 10
  dd <- 2 * diag(d)
  dd[2, 1] <- -1
  dd[1, 2] <- -1
  dd[3, 1] <- -1
  dd[1, 3] <- -1
  delta_time <- 0.01
  corr_bm <- correlated_brownian_noise(
    sigma_matrix = dd, n = n, delta_time = delta_time
  )
  testthat::expect_equal(
    cov(corr_bm)[1:5, 1:5] / delta_time,
    cov(mvtnorm::rmvnorm(n, sigma = dd))[1:5, 1:5],
    tolerance = 0.05
  )
})

test_that("compound_poisson_jumps__shape", {
  set.seed(42)
  n <- 100
  d <- 3
  delta_time <- 0.001
  n_jumps <- 100
  jump_vals <- matrix(
    rep(rnorm(n_jumps, 0, sd = 0.1), d),
    nrow = n_jumps, ncol = d
  )
  cmpnd_poisson <- compound_poisson_jumps(
    d = d, n = n, delta_time = delta_time, jump_values = jump_vals
  )
  testthat::expect_equal(dim(cmpnd_poisson$noise), c(n, d))
  testthat::expect_equal(dim(cmpnd_poisson$jump_times), c(n_jumps, d))

  # synchro
  cmpnd_poisson <- compound_poisson_jumps(
    d = d, n = n, delta_time = delta_time, jump_values = jump_vals,
    synchronised = T
  )
  testthat::expect_equal(dim(cmpnd_poisson$noise), c(n, d))
  testthat::expect_true(is.vector(cmpnd_poisson$jump_times))
  testthat::expect_equal(length(cmpnd_poisson$jump_times), n_jumps)
})

test_that("compound_poisson_jumps__NAs", {
  set.seed(42)
  n <- 1000
  d <- 3
  delta_time <- 0.001
  n_jumps <- 10000
  jump_vals <- NA
  cmpnd_poisson <- compound_poisson_jumps(
    d = d, n = n, delta_time = delta_time, jump_values = jump_vals
  )

  testthat::expect_equal(dim(cmpnd_poisson$noise), c(n, d))
  testthat::expect_equal(dim(cmpnd_poisson$jump_times), c(1, d))
})

test_that("compound_poisson_jumps__vector", {
  set.seed(42)
  n <- 1000
  d <- 1
  delta_time <- 0.001
  n_jumps <- 100
  jump_vals <- rnorm(n_jumps, 0, sd = 0.1)
  cmpnd_poisson <- compound_poisson_jumps(
    d = d, n = n, delta_time = delta_time, jump_values = jump_vals
  )

  testthat::expect_true(is.vector(cmpnd_poisson$noise))
  testthat::expect_true(is.vector(cmpnd_poisson$jump_times))
  testthat::expect_equal(length(cmpnd_poisson$jump_times), n_jumps)
})

test_that("compound_poisson_jumps__n_jumps_higher", {
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
  corr_jumps <- correlated_jumps(
    n = n, sigma = corr_mat, delta_time = delta_time
  )

  testthat::expect_equal(dim(corr_jumps), c(n, d))
})

test_that("bm_compound_poisson__shape", {
  set.seed(42)
  n <- 1000
  d <- 3
  delta_time <- 0.001
  corr_mat <- matrix(c(1.0, 0.2, 0.2, 0.2, 1.0, 0.2, 0.2, 0.2, 1.0), 3, 3)
  lapply(
    c(-10, 0, 1, 10, 100, 1000),
    function(n_j) {
      bm_cp <- bm_compound_poisson(
        n = n, sigma = corr_mat, jump_sigma = diag(d),
        n_jumps = n_j, delta_time = delta_time
      )
      testthat::expect_equal(dim(bm_cp), c(n, d))
    }
  )
})

test_that("bm_compound_poisson_ghyp__shape", {
  set.seed(42)
  n <- 1000
  d <- 3
  delta_time <- 0.001
  corr_mat <- matrix(c(1.0, 0.2, 0.2, 0.2, 1.0, 0.2, 0.2, 0.2, 1.0), 3, 3)
  jump_distr <- ghyp::ghyp(mu = rep(0, d))
  lapply(
    c(-10, 0, 1, 10, 100, 1000),
    function(n_j) {
      bm_cp <- bm_compound_poisson_ghyp(
        n = n, sigma = corr_mat, ghyp_distr = jump_distr,
        n_jumps = n_j, delta_time = delta_time
      )
      testthat::expect_equal(dim(bm_cp), c(n, d))
    }
  )
})
