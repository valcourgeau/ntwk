
test_that("utils_RowNormalised__shape", {
  n <- 10000
  d <- 10

  adj <- diag(5, d, d)
  row_adj <- RowNormalised(adj)
  testthat::expect_equal(dim(row_adj), c(d, d))
})

test_that("utils_RowNormalised__zero_diag", {
  d <- 10
  adj <- diag(5, d, d)
  row_adj <- RowNormalised(adj)
  testthat::expect_equal(sum(abs(diag(row_adj))), 0.0)
})

test_that("utils_RowNormalised__off_diag", {
  d <- 10
  adj <- matrix(1, d, d)
  row_adj <- RowNormalised(adj)
  testthat::expect_equal(sum(abs(row_adj)), d*(d-1)/(d-1))
  testthat::expect_equal(row_adj[1, d-1], 1/(d-1))
  testthat::expect_true(all(rowSums(row_adj) == 1.))
})

test_that("utils_ConcatenateCol__shape_and_values", {
  d <- 10
  col_to_concat <- seq(d * d)
  for(i in 1:d){
    row_adj <- ConcatenateCol(col_to_concat, n = i)
    testthat::expect_equal(dim(row_adj), c(d * d, i))
    testthat::expect_true(all(
        apply(row_adj, 2, function(x){x == col_to_concat})
      )
    )
  }
})

test_that("utils_AugmentedDiag__shape_values", {
  d <- 10
  col_to_concat <- seq(d * d)
  for(i in seq(2 * d)){
    if(i < d){
      aug_diag <- AugmentedDiag(d = d, offset = i)
      testthat::expect_equal(dim(aug_diag), c(d, d))
      testthat::expect_equal(sum(abs(aug_diag)), d-i)
      testthat::expect_true(all(aug_diag[aug_diag!=0] == 1))
    }else{
      testthat::expect_error(AugmentedDiag(d = d, offset = i))
    }
  }
})
