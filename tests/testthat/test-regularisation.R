
testthat::test_that("get_reg_fn__type", {
  set.seed(45)
  n <- 10000
  d <- 4
  adj_test <- diag(d)
  adj_test[2, 1] <- 0.5
  adj_test[3, 1] <- 0.1

  reg_fn <- get_reg_fn("l1")
  x <- sample(seq(-100, 100, by = 1), size = 30, replace = FALSE)
  testthat::expect_equal(reg_fn(x), sum(abs(x)))

  reg_fn <- get_reg_fn("l2")
  x <- sample(seq(-100, 100, by = 1), size = 30, replace = FALSE)
  testthat::expect_equal(reg_fn(x), sqrt(sum(x^2)))

  beta_value <- 0.4999
  delta_time <- 0.1
  y_init <- rep(0, d)

  times <- seq(0, by = delta_time, length.out = n)
  noise <- matrix(rnorm(n * d, sd = sqrt(delta_time)), ncol = d)
  path <- construct_path(
    nw_topo = adj_test, noise = noise, y_init = y_init, delta_time = delta_time
  )

  grou_mle_fit <- grou_mle(
    times = times, data = path, thresholds = NA,
    mode = "node", output = "vector"
  )

  reg_fn <- get_reg_fn("adaptive", mle = grou_mle_fit, gamma = 1)
  x <- sample(seq(-100, 100, by = 1), size = 30, replace = FALSE)
  testthat::expect_equal(
    reg_fn(grou_mle_fit),
    sum(abs(rep(1, length(grou_mle_fit))))
  )
  reg_fn <- get_reg_fn("adaptive", mle = grou_mle_fit, gamma = 3)
  testthat::expect_equal(
    reg_fn(grou_mle_fit),
    sum(abs(grou_mle_fit^{
      1 - 3
    }))
  )
})

testthat::test_that("get_reg_fn__errors", {
  set.seed(45)
  n <- 500
  d <- 4
  adj_test <- diag(d)
  adj_test[2, 1] <- 0.5
  adj_test[3, 1] <- 0.1
  delta_time <- 0.01
  y_init <- rep(0, d)
  times <- seq(0, by = delta_time, length.out = n)
  noise <- matrix(rnorm(n * d, sd = sqrt(delta_time)), ncol = d)
  path <- construct_path(
    nw_topo = adj_test, noise = noise, y_init = y_init, delta_time = delta_time
  )

  grou_mle_fit <- grou_mle(
    times = times, data = path, thresholds = NA,
    mode = "node", output = "vector"
  )

  testthat::expect_error(
    get_reg_fn("adaptive", mle = grou_mle_fit, gamma = NA),
    regexp = "gamma"
  )
  testthat::expect_error(
    get_reg_fn("adaptive", mle = NA, gamma = NA),
    regexp = "mle"
  )
})


testthat::test_that("grou_regularisation__shape", {
  set.seed(23)
  n <- 10000
  d <- 4
  adj_test <- diag(d)
  adj_test[2, 1] <- 0.5
  adj_test[3, 1] <- 0.1

  corr_mat <- matrix(c(1.0, 0.2, 0.2, 0.2, 1.0, 0.2, 0.2, 0.2, 1.0), 3, 3)
  beta_value <- 0.4999
  delta_time <- 0.1
  y_init <- rep(0, d)

  times <- seq(0, by = delta_time, length.out = n)
  noise <- matrix(rnorm(n * d, sd = sqrt(delta_time)), ncol = d)
  path <- construct_path(
    nw_topo = adj_test, noise = noise, y_init = y_init, delta_time = delta_time
  )

  reg_adj <- grou_regularisation(
    times = times,
    data = path,
    thresholds = NA,
    lambda = 10,
    reg = "l1",
    output = "matrix",
    gamma = 3, use_scaling = TRUE, cut_off = .65
  )
  testthat::expect_equal(dim(reg_adj), c(d, d))

  grou_mle_fit <- grou_mle(
    times = times, data = path, thresholds = NA,
    mode = "node", output = "matrix"
  )
  testthat::expect_equal(
    mean(reg_adj[adj_test > 0] / adj_test[adj_test > 0]),
    1.0,
    tolerance = .11
  )

  testthat::expect_equal(reg_adj, adj_test, tolerance = .2)
})


testthat::test_that("grou_regularisation__errors", {
  set.seed(23)
  n <- 100
  d <- 4
  adj_test <- diag(d)
  adj_test[2, 1] <- 0.5
  adj_test[3, 1] <- 0.1

  beta_value <- 0.4999
  delta_time <- 0.1
  y_init <- rep(0, d)

  times <- seq(0, by = delta_time, length.out = n)
  noise <- matrix(rnorm(n * d, sd = sqrt(delta_time)), ncol = d)
  path <- noise

  testthat::expect_error(grou_regularisation(
    times = times, data = path, thresholds = NA, lambda = 10, reg = "l1",
    output = "qwe", gamma = 3, use_scaling = TRUE, cut_off = .65
  ), regexp = "output")
  testthat::expect_error(grou_regularisation(
    times = times, data = path, thresholds = NA, lambda = -1, reg = "l1",
    output = "matrix", gamma = 3, use_scaling = TRUE, cut_off = .65
  ), regexp = "lambda")
  testthat::expect_error(grou_regularisation(
    times = times, data = path, thresholds = NA, lambda = 10, reg = "l1",
    output = "cut_off", gamma = 3, use_scaling = TRUE, cut_off = -1
  ))
})
