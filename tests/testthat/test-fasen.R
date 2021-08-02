
test_that("fasen__shape", {
  n <- 10000
  d <- 10
  set.seed(42)
  data <- matrix(
    data = runif(n * d, -1, 10),
    nrow = n,
    ncol = d
  )
  fasen_reg <- ntwk::FasenRegression(data)
  testthat::expect_equal(dim(fasen_reg), c(d, d))
})
