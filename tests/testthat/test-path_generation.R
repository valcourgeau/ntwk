
test_that("CorrelatedBrownianNoise_shape", {
  set.seed(42)

  corr_mat <- matrix(c(1.0, 0.2, 0.2, 0.2, 1.0, 0.2, 0.2, 0.2, 1.0), 3, 3)
  d <- dim(corr_mat)[1]
  n <- 1000
  delta_time <- 0.01
  corr_inc <- CorrelatedBrownianNoise(sigma_matrix = corr_mat, n = n, delta_time = delta_time)
  testthat::expect_equal(dim(corr_inc), c(n, d))
})

test_that("CorrelatedBrownianNoise_corr", {
  set.seed(42)

  corr_mat <- matrix(c(1.0, 0.2, 0.2, 0.2, 1.0, 0.2, 0.2, 0.2, 1.0), 3, 3)
  d <- dim(corr_mat)[1]
  n <- 2000
  delta_time <- 0.01
  corr_inc <- CorrelatedBrownianNoise(sigma_matrix = corr_mat, n = n, delta_time = delta_time)
  testthat::expect_equal(corr_mat, cor(corr_inc), tolerance = .5)
})

