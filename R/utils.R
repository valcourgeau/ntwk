#' Normalise a square matrix by dividing elements by the sum of
#' off-diagonal entries row-wise. Note: sets diagonal to zero.
#' @param adj Adjacency matrix to normalise
#' @return A matrix based on `adj` with zero diagonals
#'     and off-diagonal elements sum to one.
#' @example row_normalised(matrix(1, 5, 5))
row_normalised <- function(adj) {
  # zero on diag, row sum to 1 with positive weights

  diag(adj) <- 0.0
  adj[adj != 0.0] <- 1.0
  divisors <- pmax(rowSums(adj), 1.0)

  return(diag(1.0 / divisors) %*% adj)
}

#' Concatenate `col_vec` columnwise `n` times
#' @param col_vec a vector to copy/concatenate
#' @param n number of times to concatenate (number of columns)
#' @return A matrix where `col_vec` is concatenated
#'     `n` times along the column axis.
#' @example concatenate_col(seq(10), 5) # (10, 5) matrix created
concatenate_col <- function(col_vec, n) {
  return(matrix(rep(col_vec, n), ncol = n))
}

#' Clean data using STL / Loess function.
#' @param data Data to clean.
#' @param frequency Time series frequency of observations
#' @param s_window Seasonal window.
#' @param t_window Trend window.
#' @param ... additional inputs given to stats::stl
#' @return List with the cleaned time series, the remainders
#'     of the cleaning and their standard deviations.
#' @importFrom stats stl
#' @export
clean_data <- function(data, frequency = 24,
                       s_window = 24, t_window = 24, ...) {
  stl_cleaned <- lapply(seq_len(ncol(data)),
    FUN = function(i) {
      stl(
        ts(data[, i], frequency = frequency),
        s.window = s_window, t.window = t_window,
        ...
      )
    }
  )
  names(stl_cleaned) <- colnames(data)
  remainders <- vapply(stl_cleaned, function(x) {
    as.numeric(x$time.series[, 3])
  }, rep(0, length(data[, 1])))
  colnames(remainders) <- colnames(data)
  std_dev <- apply(remainders, MARGIN = 2, sd)

  return(
    list(
      stl_obj = stl_cleaned,
      remainders = remainders,
      std.dev = std_dev
    )
  )
}

#' Returns an offset diagonal matrix.
#' @param d Matrix dimensional (square matrix)
#' @param offset Offset from the diagonal
#' @return A matrix with an offset diagonal of ones.
#' @example augmented_diag(10, 1) # zeroes everywhere except 1s on the +1 diag
#' @export
augmented_diag <- function(d, offset) {
  assertthat::assert_that(d > offset)

  mat_temp <- cbind(
    matrix(0, nrow = d, ncol = abs(offset)),
    rbind(
      diag(d - abs(offset)),
      matrix(0, ncol = d - abs(offset), nrow = abs(offset))
    )
  )
  if (offset < 0) {
    return(t(mat_temp))
  } else {
    return(mat_temp)
  }
}
