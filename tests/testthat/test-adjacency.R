testthat::test_that("Adjacency_IsolatedNetwork__shapes", {
  d <- 10
  theta_2 <- 2
  adj <- IsolatedNetwork(d = d, theta_2 = theta_2)
  testthat::expect_equal(dim(adj), c(d, d))
})

testthat::test_that("Adjacency_IsolatedNetwork__values", {
  d <- 10
  theta_2 <- 2
  adj <- IsolatedNetwork(d = d, theta_2 = theta_2)
  testthat::expect_equal(sum(abs(diag(adj))), d * theta_2)
  diag(adj) <- 0
  testthat::expect_equal(sum(abs(adj)), 0)
})

testthat::test_that("Adjacency_PolymerNetwork__shapes", {
  d <- 10
  theta_1 <- 1
  theta_2 <- 2
  adj <- PolymerNetwork(d = d, theta_1 = theta_1, theta_2 = theta_2)
  testthat::expect_equal(dim(adj), c(d, d))
})

testthat::test_that("Adjacency_PolymerNetwork__values", {
  d <- 10
  theta_1 <- 1.5
  theta_2 <- 2
  adj <- PolymerNetwork(d = d, theta_1 = theta_1, theta_2 = theta_2)
  testthat::expect_equal(sum(abs(diag(adj))), d * theta_2)
  diag(adj) <- 0
  testthat::expect_equal(sum(abs(adj)), d * theta_1)
})

testthat::test_that("Adjacency_LatticeNetwork__shapes_values", {
  d <- 10
  theta_1 <- 1.2
  theta_2 <- 2
  adj <- LatticeNetwork(d = d, theta_1 = theta_1, theta_2 = theta_2)
  testthat::expect_equal(dim(adj), c(d, d))
})

testthat::test_that("Adjacency_LatticeNetwork__values", {
  d <- 10
  theta_1 <- 1.2
  theta_2 <- 2
  adj <- LatticeNetwork(d = d, theta_1 = theta_1, theta_2 = theta_2)
  testthat::expect_equal(sum(abs(diag(adj))), d * theta_2)
  diag(adj) <- 0
  testthat::expect_equal(
    sum(abs(adj) * rowSums(adj > 0)) / sum(abs(adj) > 0),
    theta_1
  )
  testthat::expect_true(all(rowSums(adj / theta_1)[1:(d - 1)] == 1))
})
