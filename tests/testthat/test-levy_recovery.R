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
  noise <- ghyp::rghyp(n, ghyp::ghyp(mu = rep(0, d)))

  lapply(
    seq_len(length(ghyp_names)),
    function(n_pick) {
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
  n <- 5000
  d <- 2
  data <- matrix(rnorm(n * d), ncol = d)
  testthat::expect_error(fit_ghyp_diffusion(data, "Normal", silent = T))
  testthat::expect_error(
    fit_ghyp_diffusion(data, c("NIG", "Normal"), silent = T)
  )
})
