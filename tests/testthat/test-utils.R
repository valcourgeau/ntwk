
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
  row_adj <- row_normalised(adj)
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
    row_adj <- row_normalised(adj, keep_value = T)
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
    row_adj <- row_normalised(adj, keep_value = F)
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

test_that("utils_augmented_diag__shape_values", {
  d <- 10
  col_to_concat <- seq(d * d)
  for (i in seq(2 * d)) {
    if (i < d) {
      aug_diag <- augmented_diag(d = d, offset = i)
      testthat::expect_equal(dim(aug_diag), c(d, d))
      testthat::expect_equal(sum(abs(aug_diag)), d - i)
      testthat::expect_true(all(aug_diag[aug_diag != 0] == 1))
    } else {
      testthat::expect_error(augmented_diag(d = d, offset = i))
    }
  }
})
