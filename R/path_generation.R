
#' Generates a single Graph Ornstein-Uhlenbeck path with the given
#' network topology (`nw_topo`), noise increments (`noise`),
#' initial values (`y_init`) and unique time difference (`delta_time`).
#' @param nw_topo Graph topology or adjacency matrix.
#' @param noise Noise increments to apply.
#' @param y_init Start value.
#' @param delta_time Time step.
#' @examples
#' n <- 100
#' d <- 10
#' delta_time <- 0.01
#' noise <- matrix(rnorm(n * d), ncol = d)
#' construct_path(diag(10), noise, d, delta_time)
#' @export
construct_path <- function(nw_topo, noise, y_init, delta_time) {
  # generates path given noise increments and fixed delta_time.
  assertthat::assert_that(delta_time > 0)
  assertthat::are_equal(ncol(nw_topo), nrow(nw_topo))

  n <- nrow(noise)
  d <- ncol(nw_topo)
  d_y <- matrix(0, nrow = n, ncol = d)
  d_y[1, ] <- y_init
  for (i in 2:n) {
    increment <- -as.vector(nw_topo %*% d_y[i - 1, ]) * delta_time + noise[i, ]
    d_y[i, ] <- d_y[i - 1, ] + increment
  }

  return(d_y)
}

#' Generates multidimensional correlated Brownian motion increments
#' @param sigma_matrix Correlation/Covariance matrix
#' @param n Sample length
#' @param delta_time Time difference between steps
#' @examples
#' set.seed(42)
#' corr_mat <- matrix(c(1.0, 0.2, 0.2, 0.2, 1.0, 0.2, 0.2, 0.2, 1.0), 3, 3)
#' n <- 1000
#' delta_time <- 0.01
#' correlated_brownian_noise(
#'   sigma_matrix = corr_mat, n = n, delta_time = delta_time
#' )
#' @importFrom mvtnorm rmvnorm
#' @export
correlated_brownian_noise <- function(n, sigma_matrix, delta_time) {
  assertthat::assert_that(ncol(sigma_matrix) > 1)
  return(mvtnorm::rmvnorm(n = n, sigma = sigma_matrix) * sqrt(delta_time))
}

compound_poisson_jumps <- function(d, n, delta_time,
                                   jump_values = NA, synchronised = F) {
  if (any(is.na(jump_values))) {
    # one jump is enough
    jump_values <- matrix(0, nrow = 1, ncol = d)
  }
  if (d == 1) {
    n_jumps <- length(jump_values)
    jump_values <- as.matrix(jump_values)
  } else {
    n_jumps <- nrow(jump_values)
  }
  if (n < n_jumps) {
    warning("`n` < `n_jumps` in compound_poisson_jumps:
            only `n` steps will be created.")
  }
  results <- matrix(0, nrow = n, ncol = d)
  assertthat::assert_that(assertthat::are_equal(ncol(jump_values), d))
  horizon <- n * delta_time

  if (synchronised) {
    # synchronise all jumps at the same times
    jump_times_tmp <- runif(min = 0, max = horizon, n = min(n, n_jumps))
    jump_times_tmp <- sort(jump_times_tmp)
    time_grid <- seq(0, to = delta_time * n, by = delta_time)
    idx_jump_injection <- vapply(
      X = jump_times_tmp,
      FUN = function(jump_time) {
        pmin(which.max(jump_time < time_grid), n)
      },
      FUN.VALUE = 3
    )
    jump_value_idx <- 1
    assertthat::assert_that(
      assertthat::are_equal(nrow(jump_values), length(idx_jump_injection))
    )

    jump_value_idx <- 1
    for (idx in idx_jump_injection) {
      results[idx, ] <- results[idx, ] + jump_values[jump_value_idx, ]
      jump_value_idx <- jump_value_idx + 1
    }
    jump_times <- jump_times_tmp
  } else {
    jump_times_collection <- lapply(seq_len(d), function(x) {
      jump_times_tmp <- runif(min = 0, max = horizon, n = min(n, n_jumps))
      return(sort(jump_times_tmp))
    })
    time_grid <- seq(0, to = delta_time * n, by = delta_time)

    idx_jump_injection <- lapply(
      jump_times_collection,
      function(times) {
        vapply(
          X = times,
          FUN = function(jump_time) {
            pmin(which.max(jump_time < time_grid), n)
          },
          FUN.VALUE = 3
        )
      }
    )
    assertthat::assert_that(
      assertthat::are_equal(length(idx_jump_injection), d)
    )

    for (i in seq_len(d)) {
      jump_value_idx <- 1
      for (idx in idx_jump_injection[[i]]) {
        results[idx, i] <- results[idx, i] + jump_values[jump_value_idx, i]
        jump_value_idx <- jump_value_idx + 1
      }
    }
    assertthat::assert_that(
      assertthat::are_equal(dim(results), c(n, d))
    )
    jump_times <- do.call(cbind, jump_times_collection)
  }

  return(
    list(
      "noise" = results,
      "jump_times" = jump_times
    )
  )
}

correlated_jumps <- function(n, sigma, delta_time) {
  return(
    correlated_brownian_noise(
      sigma_matrix = sigma, n = n, delta_time = delta_time
    )
  )
}


#' Generates a (correlated) Brownian motion path with
#' correlated but unsynchronised Gaussian jumps.
#' @param n Length of the path
#' @param sigma Correlation matrix for the Brownian part.
#' @param jump_sigma Correlation matrix for the jump part.
#' @param n_jumps Number of jumps.
#' @param delta_time Time step.
#' @param synchronised Boolean to synchronise all jumps or not.
#' @return A BM path with Correlated Gaussian jumps
#' @examples
#' n <- 1000
#' sigma <- matrix(c(1.0, 0.2, 0.2, 0.2, 1.0, 0.2, 0.2, 0.2, 1.0), 3, 3)
#' jump_sigma <- sigma
#' n_jumps <- 50
#' delta_time <- 0.5
#' bm_compound_poisson(n, sigma, jump_sigma, n_jumps, delta_time)
#' @export
bm_compound_poisson <- function(n, sigma, jump_sigma,
                                n_jumps, delta_time, synchronised = F) {
  d <- ncol(sigma)
  if (n_jumps > 0) {
    jump_vals <- correlated_jumps(n = n_jumps, sigma = jump_sigma, delta_time)
  } else {
    jump_vals <- NA
  }
  cmpnd_poisson <- compound_poisson_jumps(
    d = d, n = n, delta_time = delta_time,
    jump_values = jump_vals, synchronised = synchronised
  )
  cmpnd_poisson_noise <- cmpnd_poisson$noise
  corr_bm_noise <- correlated_brownian_noise(
    sigma_matrix = sigma, n = n, delta_time = delta_time
  )

  return(cmpnd_poisson_noise + corr_bm_noise)
}


#' Generates a (correlated) Brownian motion path with
#' correlated but unsynchronised GHYP jumps.
#' @param n Length of the path
#' @param sigma Correlation matrix for the Brownian part.
#' @param ghyp_distr GHYP distribution (`ghyp::ghyp`).
#' @param n_jumps Number of jumps.
#' @param delta_time Time step.
#' @return A BM path with Correlated Gaussian jumps
#' @importFrom ghyp rghyp
#' @examples
#' n <- 1000
#' d <- 3
#' sigma <- matrix(c(1.0, 0.2, 0.2, 0.2, 1.0, 0.2, 0.2, 0.2, 1.0), 3, 3)
#' jump_sigma <- sigma
#' n_jumps <- 50
#' delta_time <- 0.5
#' ghyp_distr <- ghyp::ghyp(mu = rep(0, d))
#' bm_compound_poisson_ghyp(n, sigma, ghyp_distr, n_jumps, delta_time)
#' @export
bm_compound_poisson_ghyp <- function(n, sigma, ghyp_distr,
                                     n_jumps, delta_time) {
  d <- ncol(sigma)
  jump_vals <- ghyp::rghyp(n = n_jumps, object = ghyp_distr)
  cmpnd_poisson <- compound_poisson_jumps(
    d = d, n = n, delta_time = delta_time,
    jump_values = jump_vals, synchronised = T
  )
  cmpnd_poisson_noise <- cmpnd_poisson$noise
  corr_bm_noise <- correlated_brownian_noise(
    sigma_matrix = sigma, n = n, delta_time = delta_time
  )

  return(cmpnd_poisson_noise + corr_bm_noise)
}
