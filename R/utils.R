#' Normalise a square matrix by dividing elements by the sum of
#' off-diagonal entries row-wise. Note: sets diagonal to zero.
#' @param adj Adjacency matrix to normalise
#' @param keep_values Whether to keep off-diagonal values
#'     or normalise fully (defaults to `FALSE`).
#' @return A matrix based on `adj` with zero diagonals
#'     and off-diagonal elements sum to one.
#' @examples
#' d <- 5
#' adj <- matrix(1, d, d)
#' row_normalised(adj)
#' @export
row_normalised <- function(adj, keep_values = F) {
  # zero on diag, row sum to 1 with positive weights

  diag(adj) <- 0.0
  if (!keep_values) adj[abs(adj) != 0.0] <- 1.0
  divisors <- pmax(rowSums(adj != 0.0), 1.0)

  return(diag(1.0 / divisors) %*% adj)
}

#' Concatenate `col_vec` columnwise `n` times
#' @param col_vec a vector to copy/concatenate
#' @param n number of times to concatenate (number of columns)
#' @return A matrix where `col_vec` is concatenated
#'     `n` times along the column axis.
#' @examples
#' vec_to_concat <- seq(10)
#' n <- 5
#' concatenate_col(vec_to_concat, n) # (10, 5) matrix created
#' @export
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
#' @examples
#' n <- 1000
#' d <- 3
#' delta_time <- 0.01
#' data <- matrix(
#'   rnorm(n = n * d, mean = 0, sd = sqrt(delta_time)),
#'   ncol = d
#' )
#' data <- apply(data, 2, cumsum)
#' clean_data(data)
#' @importFrom stats stl rnorm
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
#' @examples
#' d <- 10
#' offset <- 1
#' augmented_diag(d, offset) # zeroes everywhere except 1s on the +1 diag
#' @export
augmented_diag <- function(d, offset) {
  assertthat::assert_that(d > abs(offset))

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


#' Filters the data above the thresholds (set to zero).
#' Thresholds can be a single value, or equal to the number of rows.
#' If `thresholds` is set to `NA`, no thresholding happens.
#' @param data Data to filter. Can be vector or matrix.
#' @param thresholds Filtering thresholds. Defaults to `NA` but can be
#'     a mixture of `NA`s and numerical values
#'     (with length equal to `nrow(data)`).
#' @return List of filtered data vector or matrix (`$data`)
#'     and binary matrix of filtered times (`filter`) (`T` means filtered).
#' @importFrom stats rnorm
#' @examples
#' n <- 10000
#' d <- 5
#' data <- matrix(rnorm(n * d), ncol = d)
#'
#' # No filtering
#' data_filtering(data)
#'
#' # Filtering everything
#' data_filtering(data, 0.0)
#' data_filtering(data, rep(0.0, d))
#'
#' # Filtering only the first col
#' data_filtering(data, c(0.0, rep(NA, d - 1)))
#' @export
data_filtering <- function(data, thresholds = NA) {
  assertthat::assert_that(any(c(is.na(thresholds), is.numeric(thresholds))))
  if (is.vector(data)) {
    if (all(is.na(thresholds))) {
      filtered_data <- data
      binary_data <- rep(F, length(data))
    } else {
      assertthat::assert_that(
        assertthat::are_equal(length(thresholds), 1)
      )
      binary_data <- abs(data) > thresholds
      binary_data[is.na(binary_data)] <- F
      filtered_data <- data * !binary_data
    }
  } else {
    if (is.matrix(data)) {
      if (length(thresholds) == 1) thresholds <- rep(thresholds, ncol(data))
      assertthat::assert_that(
        assertthat::are_equal(length(thresholds), ncol(data))
      )
      binary_data <- t(apply(abs(data), 1, ">", thresholds))
      binary_data[is.na(binary_data)] <- F
      filtered_data <- data * !binary_data
    } else {
      stop("`data` should be a matrix or a vector.")
    }
  }

  return(list(data = filtered_data, filter = binary_data))
}
