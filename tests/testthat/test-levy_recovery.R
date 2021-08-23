testthat::test_that("levy_recovery_shapes", {
  n <- 1000
  d <- 2
  times <- seq(n)
  delta_time <- 0.01
  noise <- ghyp::rghyp(n, ghyp::ghyp(mu = rep(0, d)))
  data <- construct_path(
    diag(d),
    noise = noise, y_init = rep(0, d), delta_time = delta_time
  )
  levy_recov <- levy_recovery(adj = diag(d), data = data, times = times)
  dim_increments <- c(dim(data)[1] - 1, dim(data)[2])
  testthat::expect_equal(dim(levy_recov$increments), dim_increments)
  testthat::expect_equal(dim(levy_recov$cumsum), dim_increments)
  testthat::expect_equal(c(d, d), dim(levy_recov$adj))
})

testthat::test_that("levy_recovery_look_ahead", {
  # TODO(val) change this if look_ahead > 1 is implemented
  n <- 1000
  d <- 2
  look_ahead_max <- 10
  times <- seq(n)
  delta_time <- 0.01
  noise <- ghyp::rghyp(n, ghyp::ghyp(mu = rep(0, d)))
  data <- construct_path(
    diag(d),
    noise = noise, y_init = rep(0, d), delta_time = delta_time
  )
  for (i in seq_len(look_ahead_max)) {
    if (i > 1) {
      testthat::expect_error(
        levy_recovery(
          adj = diag(d), data = data, times = times, look_ahead = i
        )
      )
    } else {
      levy_recovery(adj = diag(d), data = data, times = times, look_ahead = i)
    }
  }
})

testthat::test_that("levy_recovery_sparse_adj", {
  n <- 5000
  d <- 2
  times <- seq(n)
  delta_time <- 0.01
  noise <- ghyp::rghyp(n, ghyp::ghyp(mu = rep(0, d)))
  data <- construct_path(
    diag(d),
    noise = noise, y_init = rep(0, d), delta_time = delta_time
  )
  sparse_adj <- Matrix::Matrix(diag(d), sparse = T)
  levy_recov <- levy_recovery(adj = sparse_adj, data = data, times = times)
  testthat::expect_equal(class(sparse_adj)[1], "ddiMatrix")
  testthat::expect_equal(
    class(levy_recov$increments), class(matrix(seq_len(d), 1, d))
  )
})

testthat::test_that("fit_bm_compound_poisson_output", {
  n <- 5000
  d <- 10
  times <- seq(n)
  delta_time <- 0.01

  noise <- ghyp::rghyp(n, ghyp::ghyp(mu = rep(0, d)))
  data <- noise

  bm_cp <- fit_bm_compound_poisson(data, delta_time, jump_quantile = 0.5)
  testthat::expect_equal(dim(bm_cp$sigma), c(d, d))
  testthat::expect_true(is.na(bm_cp$jump_sigma))
  testthat::expect_true(is.na(bm_cp$n_jumps))
  testthat::expect_true(is.na(bm_cp$poisson))
  testthat::expect_true(is.na(bm_cp$poisson_unique))
  testthat::expect_true(is.na(bm_cp$ghyp))

  bm_cp <- fit_bm_compound_poisson(
    data, delta_time,
    thresholds = rep(0.1, d), jump_quantile = 0.05
  )
  testthat::expect_equal(dim(bm_cp$sigma), c(d, d))
  testthat::expect_equal(dim(bm_cp$jump_sigma), c(d, d))
  testthat::expect_equal(length(bm_cp$poisson), d)
  testthat::expect_equal(length(bm_cp$poisson_unique), 1)
  testthat::expect_true(bm_cp$poisson_unique > 0)
  testthat::expect_equal(class(bm_cp$ghyp$FULL)[1], "mle.ghyp")
})

testthat::test_that("fit_ghyp_diffusion_output", {
  n <- 500
  d <- 2
  times <- seq(n)
  delta_time <- 0.01
  ghyp_model_list <- list(
    "NIG" = ghyp::fit.NIGmv,
    "GAUSS" = ghyp::fit.gaussmv,
    "VG" = ghyp::fit.VGmv,
    "T" = ghyp::fit.tmv,
    "FULL" = ghyp::fit.ghypmv
  )

  ghyp_names <- names(ghyp_model_list)
  set.seed(42)

  lapply(
    seq_len(length(ghyp_names)),
    function(n_pick) {
      # noise generation
      noise <- ghyp::rghyp(n, ghyp::ghyp(mu = rep(0, d)))

      # pick n model families and fit
      fit_ghyp <- fit_ghyp_diffusion(
        noise,
        ghyp_names = sample(ghyp_names, size = n_pick, replace = F),
        silent = T
      )
      testthat::expect_equal(length(fit_ghyp), n_pick)
    }
  )
})

testthat::test_that("fit_ghyp_diffusion_wrong_name", {
  n <- 1000
  d <- 2
  data <- matrix(rnorm(n * d), ncol = d)
  testthat::expect_error(fit_ghyp_diffusion(data, "Normal", silent = T))
  testthat::expect_error(
    fit_ghyp_diffusion(data, c("NIG", "Normal"), silent = T)
  )
})

testthat::test_that("bi_power_variation_shapes", {
  n <- 100
  d <- 2
  delta_time <- 0.01
  data <- matrix(rnorm(n * d), ncol = d)
  dt <- bi_power_variation(data, delta_time)
  testthat::expect_equal(dim(dt), c(d, d))
})
