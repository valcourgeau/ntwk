
test_that("utils_row_normalised__shape", {
  n <- 10000
  d <- 10

  adj <- diag(5, d, d)
  row_adj <- row_normalised(adj)
  testthat::expect_equal(dim(row_adj), c(d, d))
})

test_that("utils_row_normalised__zero_diag", {
  d <- 10
  adj <- diag(5, d, d)
  row_adj <- row_normalised(adj, keep_values = FALSE)
  testthat::expect_equal(sum(abs(diag(row_adj))), 0.0)
})

test_that("utils_row_normalised__off_diag", {
  d <- 10
  adj <- matrix(1, d, d)
  row_adj <- row_normalised(adj)
  testthat::expect_equal(sum(abs(row_adj)), d * (d - 1) / (d - 1))
  testthat::expect_equal(row_adj[1, d - 1], 1 / (d - 1))
  testthat::expect_true(all(rowSums(row_adj) == 1.))
})

test_that("utils_row_normalised__keep_value", {
  d <- 3
  theta_2 <- 2.0
  thetas_1 <- seq(-1, 1, length.out = 10)

  for (theta_1 in thetas_1) {
    adj <- matrix(1, d, d)
    diag(adj) <- 0
    adj <- theta_1 * adj + theta_2 * diag(d)
    row_adj <- row_normalised(adj, keep_value = TRUE)
    testthat::expect_equal(sum(abs(diag(row_adj))), 0)
    testthat::expect_equal(
      apply(row_adj, 1, function(x) sum(x)),
      rep(theta_1, d)
    )
  }

  for (theta_1 in thetas_1) {
    adj <- matrix(1, d, d)
    diag(adj) <- 0
    adj <- theta_1 * adj + theta_2 * diag(d)
    row_adj <- row_normalised(adj, keep_value = FALSE)
    testthat::expect_equal(sum(abs(diag(row_adj))), 0)
    testthat::expect_equal(
      apply(row_adj, 1, function(x) sum(x)),
      rep(1, d)
    )
  }
})

test_that("utils_concatenate_col__shape_and_values", {
  d <- 10
  col_to_concat <- seq(d * d)
  for (i in 1:d) {
    row_adj <- concatenate_col(col_to_concat, n = i)
    testthat::expect_equal(dim(row_adj), c(d * d, i))
    testthat::expect_true(all(
      apply(row_adj, 2, function(x) {
        x == col_to_concat
      })
    ))
  }
})

test_that("utils_clean_data__shape_and_values", {
  n <- 200
  d <- 10
  freq <- 24
  delta_time <- 0.01
  noise <- matrix(rnorm(n * d, 0, sd = sqrt(delta_time)), ncol = d)
  seasonality <- 2 * sin(2 * pi * seq_len(n) / freq)
  seasonality <- seasonality - 1 * cos(2 * pi * seq_len(n) / freq)
  data <- noise + seasonality
  cln_dt <- clean_data(data)
  lapply(
    cln_dt$stl_obj,
    function(x) {
      testthat::expect_equal(dim(x$time.series), c(n, 3))
    }
  )
  testthat::expect_equal(dim(cln_dt$remainders), c(n, d))
  testthat::expect_equal(length(cln_dt$std.dev), d)
})



test_that("utils_augmented_diag__shape_values", {
  d <- 10
  col_to_concat <- seq(d * d)
  for (i in seq(-2 * d, 2 * d)) {
    if (-d < i & i < d) {
      aug_diag <- augmented_diag(d = d, offset = i)
      testthat::expect_equal(dim(aug_diag), c(d, d))
      testthat::expect_equal(sum(abs(aug_diag)), d - abs(i))
      testthat::expect_true(all(aug_diag[aug_diag != 0] == 1))
    } else {
      testthat::expect_error(augmented_diag(d = d, offset = i))
    }
  }
})

test_that("utils_data_filtering_diag__shape_values", {
  n <- 10000
  d <- 5
  data <- matrix(rnorm(n * d), ncol = d)

  # No filtering
  testthat::expect_equal(data_filtering(data)$data, data)

  # Filtering everything
  expected_vals <- matrix(0, nrow = nrow(data), ncol = ncol(data))
  testthat::expect_equal(data_filtering(data, 0.0)$data, expected_vals)
  testthat::expect_equal(data_filtering(data, rep(0.0, d))$data, expected_vals)
  testthat::expect_equal(
    dim(data_filtering(data, rep(0.0, d))$filter), c(n, d)
  )
  testthat::expect_error(data_filtering(data, rep(0.0, d * 2)))

  # Filtering only the first col
  tmp <- data_filtering(data, c(0.0, rep(NA, d - 1)))$data
  testthat::expect_equal(tmp[, 2:d], data[, 2:d])
  testthat::expect_equal(tmp[, 1], rep(0, n))

  # One filter
  tmp <- data_filtering(data, NA)$data
  testthat::expect_equal(tmp, data)
})

test_that("utils_data_filtering_vector", {
  n <- 10000
  d <- 1
  data <- rnorm(n * d)

  # No filtering
  testthat::expect_equal(data_filtering(data)$data, data)

  # Filtering everything
  expected_vals <- rep(0, length(data))
  testthat::expect_equal(data_filtering(data, 0.0)$data, expected_vals)
  testthat::expect_equal(data_filtering(data, rep(0.0, d))$data, expected_vals)
  testthat::expect_equal(
    length(data_filtering(data, rep(0.0, d))$filter), n
  )
  testthat::expect_true(
    is.null(dim(data_filtering(data, rep(0.0, d))$filter))
  )
  testthat::expect_error(data_filtering(data, rep(0.0, 3)))
  testthat::expect_error(data_filtering(data, c()))
})
