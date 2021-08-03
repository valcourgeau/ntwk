testthat::test_that("levy_recovery_shapes", {
  n <- 1000
  d <- 2
  times <- seq(n)
  delta_time <- 0.01
  noise <- cbind(ghyp::rghyp(n, ghyp::ghyp()), ghyp::rghyp(n, ghyp::ghyp()))
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

testthat::test_that("levy_recovery_sdistr", {
  n <- 5000
  d <- 2
  times <- seq(n)
  delta_time <- 0.01
  noise <- cbind(ghyp::rghyp(n, ghyp::ghyp()), ghyp::rghyp(n, ghyp::ghyp()))
  data <- construct_path(
    diag(d),
    noise = noise, y_init = rep(0, d), delta_time = delta_time
  )
  levy_recov <- levy_recovery(adj = diag(d), data = data, times = times)
  print(summary(levy_recov$increments))
  print(summary(noise))
})
