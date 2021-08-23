testthat::test_that("core_node_mle_kron_numerator_unif_times", {
  n <- 5
  d <- 5
  times <- seq_len(n)
  data <- matrix(rep(seq_len(n), d), nrow = n, ncol = d, byrow = F)
  cn_mle <- core_node_mle(times = times, data = data, thresholds = NA)
  perfect_denominator <- matrix(sum(seq_len(n - 1)), d, d)
  testthat::expect_equal(cn_mle$numerator, perfect_denominator)
})

testthat::test_that("core_node_mle_kron_numerator_doubling_unif_times", {
  n <- 5
  d <- 5
  times <- seq_len(n)
  # component j have increments of size j
  data <- do.call(cbind, lapply(1:d, function(j) {
    j * seq_len(n)
  }))
  cn_mle <- core_node_mle(times = times, data = data, thresholds = NA)
  # i is for the increment, j for the component
  component_wise <- vapply(1:d, function(j) {
    sum(j * seq_len(n - 1))
  }, 1.0) # rows
  perfect_denominator <- do.call(cbind, lapply(1:d, function(i) {
    i * component_wise
  })) # increments-wise, building cols
  testthat::expect_equal(cn_mle$numerator, perfect_denominator)
})


testthat::test_that("grou_mle__shapes", {
  n_nodes <- 50
  n_sample <- 10000
  set.seed(42)
  adj_test <- diag(n_nodes)
  adj_test[2, 1] <- 0.5

  mesh_size <- 0.01
  sample_path <- construct_path(
    adj_test,
    matrix(rnorm(n_sample * n_nodes, 0, 1 * mesh_size^{
      1 / 2
    }), ncol = n_nodes),
    rep(0, n_nodes),
    mesh_size
  )
  node_mat <- grou_mle(
    times = seq(0, by = mesh_size, length.out = n_sample), data = sample_path,
    adj = adj_test, div = 1e3, mode = "node", output = "matrix"
  )
  node_vec <- grou_mle(
    times = seq(0, by = mesh_size, length.out = n_sample), data = sample_path,
    adj = adj_test, div = 1e3, mode = "node", output = "vector"
  )
  net_vec <- grou_mle(
    times = seq(0, by = mesh_size, length.out = n_sample), data = sample_path,
    adj = adj_test, div = 1e3, mode = "network", output = "vector"
  )
  net_mat <- grou_mle(
    times = seq(0, by = mesh_size, length.out = n_sample), data = sample_path,
    adj = adj_test, div = 1e3, mode = "network", output = "matrix"
  )
  testthat::expect_equal(dim(node_mat), dim(net_mat))
  testthat::expect_equal(dim(node_mat), c(n_nodes, n_nodes))

  testthat::expect_equal(length(node_vec), n_nodes * n_nodes)
  testthat::expect_equal(length(net_vec), 2)
})

testthat::test_that("grou_mle__values", {
  n_nodes <- 50
  n_sample <- 15000
  set.seed(42)
  adj_test <- diag(n_nodes)
  adj_test[2, 1] <- 0.5

  mesh_size <- 0.01
  sample_path <- construct_path(
    adj_test,
    noise = matrix(
      rnorm(n_sample * n_nodes, 0, 1 * mesh_size^{
        1 / 2
      }),
      ncol = n_nodes
    ), rep(0, n_nodes), mesh_size
  )
  times <- seq(0, by = mesh_size, length.out = n_sample)
  node_mat <- grou_mle(
    times = times, data = sample_path,
    adj = adj_test, div = 1e3, mode = "node", output = "matrix"
  )
  net_vec <- grou_mle(
    times = times, data = sample_path,
    adj = adj_test, div = 1e3, mode = "network", output = "vector"
  )
  node_mat_long <- node_mle_long(
    times = times, data = sample_path, thresholds = NA, output = "matrix"
  )
  testthat::expect_equal(net_vec, c(.5, 1.), tolerance = .4)
  testthat::expect_equal(sum(diag(node_mat)) / n_nodes, 1., tolerance = .4)
  testthat::expect_equal(sum(diag(node_mat_long)) / n_nodes, 1., tolerance = .4)
  testthat::expect_equal(node_mat_long, node_mat, tolerance = .2)

  adj_test_na <- adj_test
  adj_test_na[1, 1] <- NA
  net_vec_nas <- grou_mle(
    times = times, data = sample_path,
    adj = adj_test, div = 1e3, mode = "network", output = "vector"
  )
  testthat::expect_equal(net_vec_nas, c(.5, 1.), tolerance = .4)
})

testthat::test_that("grou_mle__zero_l2", {
  n_nodes <- 5
  n_sample <- 5000
  set.seed(42)
  adj_test <- diag(n_nodes)
  mesh_size <- 0.01
  beta_value <- 0.499
  sample_path <- construct_path(
    adj_test,
    noise = matrix(
      rnorm(n_sample * n_nodes, 0, 1 * mesh_size^beta_value),
      ncol = n_nodes
    ), rep(0, n_nodes), mesh_size
  )
  times <- seq(0, by = mesh_size, length.out = n_sample)
  net_vec <- grou_mle(
    times = times, data = sample_path,
    adj = adj_test, div = 1e3, mode = "network", output = "vector"
  )

  testthat::expect_equal(net_vec, c(0, 1.), tolerance = .4)
})

testthat::test_that("grou_mle__random_graphs", {
  n_nodes <- 5
  n_sample <- 15000
  p_link <- 0.3
  set.seed(42)

  adj_test <- as.matrix(igraph::as_adj(igraph::sample_gnp(n_nodes, p = p_link)))
  theta_1 <- 0.3
  theta_2 <- 2.0
  diag(adj_test) <- 0
  adj_test <- theta_1 * adj_test
  adj_test <- row_normalised(adj_test, keep_value = T) + theta_2 * diag(n_nodes)

  mesh_size <- 0.01
  noise <- matrix(
    rnorm(n_sample * n_nodes, 0, sqrt(mesh_size)),
    ncol = n_nodes
  )
  sample_path <- construct_path(
    adj_test,
    noise = noise, y_init = rep(0, n_nodes), delta_time = mesh_size
  )
  times <- seq(0, by = mesh_size, length.out = n_sample)

  net_vec <- grou_mle(
    times = times, data = sample_path,
    adj = adj_test, div = 1e3, mode = "network", output = "vector"
  )

  testthat::expect_equal(net_vec, c(theta_1, theta_2), tolerance = .2)
})

testthat::test_that("grou_mle__div_larger_than_sample_size", {
  n_nodes <- 50
  n_sample <- 1000
  set.seed(42)
  adj_test <- diag(n_nodes)
  adj_test[2, 1] <- 0.5

  div <- n_sample * 10

  mesh_size <- 0.01
  sample_path <- construct_path(
    adj_test,
    noise = matrix(
      rnorm(n_sample * n_nodes, 0, 1 * mesh_size^{
        1 / 2
      }),
      ncol = n_nodes
    ), rep(0, n_nodes), mesh_size
  )
  node_mat <- grou_mle(
    times = seq(0, by = mesh_size, length.out = n_sample), data = sample_path,
    adj = adj_test, div = div, mode = "node", output = "matrix"
  )
  testthat::expect_equal(dim(node_mat), c(n_nodes, n_nodes))
})

testthat::test_that("grou_mle__div_smaller_than_sample_size", {
  n_nodes <- 50
  n_sample <- 1000
  set.seed(42)
  adj_test <- diag(n_nodes)
  adj_test[2, 1] <- 0.5

  div <- n_sample / 10

  mesh_size <- 0.01
  sample_path <- construct_path(
    adj_test,
    noise = matrix(
      rnorm(n_sample * n_nodes, 0, 1 * mesh_size^{
        1 / 2
      }),
      ncol = n_nodes
    ), rep(0, n_nodes), mesh_size
  )
  node_mat <- grou_mle(
    times = seq(0, by = mesh_size, length.out = n_sample), data = sample_path,
    adj = adj_test, div = 1e4, mode = "node", output = "matrix"
  )
  testthat::expect_equal(dim(node_mat), c(n_nodes, n_nodes))
})
