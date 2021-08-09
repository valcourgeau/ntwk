testthat::test_that("Adjacency_check_thetas__shapes", {
  d <- 10
  theta_2 <- 2
  testthat::expect_error(
    check_thetas(theta_1 = 2, theta_2 = 1)
  )
  testthat::expect_true(check_thetas(theta_1 = 1, theta_2 = 2))
  testthat::expect_true(check_thetas(theta_1 = -1, theta_2 = 2))
})

testthat::test_that("Adjacency_isolated_network__shapes", {
  d <- 10
  theta_2 <- 2
  adj <- isolated_network(d = d, theta_2 = theta_2)
  testthat::expect_equal(dim(adj), c(d, d))
})

testthat::test_that("Adjacency_isolated_network__values", {
  d <- 10
  theta_2 <- 2
  adj <- isolated_network(d = d, theta_2 = theta_2)
  testthat::expect_equal(sum(abs(diag(adj))), d * theta_2)
  diag(adj) <- 0
  testthat::expect_equal(sum(abs(adj)), 0)
})

testthat::test_that("Adjacency_polymer_network__shapes", {
  d <- 10
  theta_1 <- 1
  theta_2 <- 2
  adj <- polymer_network(d = d, theta_1 = theta_1, theta_2 = theta_2)
  testthat::expect_equal(dim(adj), c(d, d))
})

testthat::test_that("Adjacency_polymer_network__values", {
  d <- 10
  theta_1 <- 1.5
  theta_2 <- 2
  adj <- polymer_network(d = d, theta_1 = theta_1, theta_2 = theta_2)
  testthat::expect_equal(sum(abs(diag(adj))), d * theta_2)
  diag(adj) <- 0
  testthat::expect_equal(sum(abs(adj)), d * theta_1)
})

testthat::test_that("Adjacency_lattice_network__shapes_values", {
  d <- 10
  theta_1 <- 1.2
  theta_2 <- 2
  adj <- lattice_network(d = d, theta_1 = theta_1, theta_2 = theta_2)
  testthat::expect_equal(dim(adj), c(d, d))
})

testthat::test_that("Adjacency_lattice_network__values", {
  d <- 10
  theta_1 <- 1.2
  theta_2 <- 2
  adj <- lattice_network(d = d, theta_1 = theta_1, theta_2 = theta_2)
  testthat::expect_equal(sum(abs(diag(adj))), d * theta_2)
  diag(adj) <- 0
  testthat::expect_equal(
    sum(abs(adj) * rowSums(adj > 0)) / sum(abs(adj) > 0),
    theta_1
  )
  testthat::expect_true(all(rowSums(adj / theta_1)[1:(d - 1)] == 1))
})

testthat::test_that("Adjacency_lattice_network__sqrt_d", {
  d <- 25
  theta_1 <- 1.2
  theta_2 <- 2
  testthat::expect_error(
    lattice_network(d = d, theta_1 = theta_1, theta_2 = theta_2)
  )
})

testthat::test_that("Adjacency_fc_network__shapes_values", {
  d <- 10
  theta_1 <- 1.2
  theta_2 <- 2
  adj <- fully_connected_network(d = d, theta_1 = theta_1, theta_2 = theta_2)
  testthat::expect_equal(dim(adj), c(d, d))
})

testthat::test_that("Adjacency_fc_network__values", {
  d <- 10
  theta_1 <- 1.2
  theta_2 <- 2
  adj <- fully_connected_network(d = d, theta_1 = theta_1, theta_2 = theta_2)
  testthat::expect_equal(sum(abs(diag(adj))), d * theta_2)
  diag(adj) <- 0
  testthat::expect_equal(
    sum(abs(adj) * rowSums(adj > 0)) / sum(abs(adj) > 0),
    theta_1
  )
  testthat::expect_true(all(rowSums(adj / theta_1)[1:(d - 1)] == 1))
})


testthat::test_that("Adjacency_fc_network__directed", {
  d <- 10
  theta_1 <- 1.2
  theta_2 <- 2
  adj <- fully_connected_network(
    d = d, theta_1 = theta_1, theta_2 = theta_2, directed = T
  )
  testthat::expect_equal(sum(adj[upper.tri(adj)]), d / 2 * theta_1)
  testthat::expect_equal(sum(adj[upper.tri(adj)]), -sum(adj[lower.tri(adj)]))
})
