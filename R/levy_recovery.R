#' Apply the Levy recovery on data with (fitted) adjacency/dynamics matrix.
#' @param adj Adjacency or dynamics matrix.
#' @param data Data to apply the Levy recovery (vector).
#' @param times Observation times (vector).
#' @param look_ahead Step increment.
#' @return List with recovered increments, adjacency matrix and cumulative sum.
#' @examples
#' n <- 1000
#' d <- 2
#' times <- seq(n)
#' delta_time <- 0.01
#' noise <- cbind(ghyp::rghyp(n, ghyp::ghyp()), ghyp::rghyp(n, ghyp::ghyp()))
#' data <- construct_path(
#'   diag(d),
#'   noise = noise, y_init = rep(0, d), delta_time = delta_time
#' )
#' levy_recovery(adj = diag(d), data = data, times = times)
#' @export
levy_recovery <- function(adj, data, times, look_ahead = 1) {
  # only for OU processes; t_{n} -> t_{n+look_ahead}
  if (look_ahead > 1) stop("NotImplementedError")

  assertthat::are_equal(length(times), nrow(data))

  diff_x <- apply(data, MARGIN = 2, FUN = function(x) {
    diff(x, lag = look_ahead)
  })
  integrated_x <- zoo::rollapply(
    cbind(data, times),
    width = look_ahead + 1,
    by.column = FALSE,
    FUN = function(sub_integ) {
      tt <- sub_integ[, ncol(sub_integ)]
      y_col_bind <- sub_integ[, -ncol(sub_integ)]
      return(
        apply(
          y_col_bind,
          MARGIN = 2,
          FUN = function(y) {
            pracma::trapz(x = tt, y = y)
          }
        )
      )
    }
  ) # create integral component-wise

  # few dimension checks
  assertthat::assert_that(nrow(integrated_x) == nrow(diff_x))
  assertthat::assert_that((nrow(data) - look_ahead) == nrow(diff_x))
  assertthat::assert_that(dim(adj)[1] == ncol(integrated_x))

  q_integrated_x <- apply(
    integrated_x,
    MARGIN = 1, FUN = function(x) {
      return(as.matrix(adj %*% x))
    }
  )
  q_integrated_x <- t(q_integrated_x)

  # TODO implement ghyp
  # TODO implement with m > 1
  assertthat::assert_that(
    assertthat::are_equal(dim(diff_x), dim(q_integrated_x))
  )

  recover <- diff_x + q_integrated_x

  recovery_results <- list()
  recovery_results[["increments"]] <- recover
  recovery_results[["cumsum"]] <- apply(recover, MARGIN = 2, FUN = cumsum)
  recovery_results[["adj"]] <- adj
  return(recovery_results)
}

#' Fit a Generalised Hyperbolic distribution.
#' @param data Data to fit.
#' @param ghyp_names Ghyp distribution names
#'     (in `c('NIG', 'GAUSS', 'VG', 'ST', 'FULL')`).
#' @param silent Boolean to silence the noise inference.
#' @param ... Extra options given to the ghyp fitting.
#' @note `ghyp` fitting is forced to be silent.
#' @return List of fitted Generalised Hyperbolic distribution.
#' @examples
#' n <- 500
#' data <- cbind(ghyp::rghyp(n, ghyp::ghyp()), ghyp::rghyp(n, ghyp::ghyp()))
#' fit_ghyp_diffusion(data = data, silent = T)
#' @export
fit_ghyp_diffusion <- function(data, ghyp_names = "FULL", silent = TRUE, ...) {
  ghyp_model_list <- list(
    "NIG" = ghyp::fit.NIGmv,
    "GAUSS" = ghyp::fit.gaussmv,
    "VG" = ghyp::fit.VGmv,
    "ST" = ghyp::fit.tmv,
    "FULL" = ghyp::fit.ghypmv
  )

  assertthat::assert_that(
    all(
      vapply(
        ghyp_names, function(x) x %in% names(ghyp_model_list),
        FUN.VALUE = TRUE
      )
    )
  )

  ghyp_models <- ghyp_model_list[ghyp_names]
  i <- 1
  res <- list()
  for (model in ghyp_models) {
    if (ghyp_names[i] == "GAUSS") {
      res[[ghyp_names[i]]] <- model(data = data, ...)
    } else {
      res[[ghyp_names[i]]] <- model(data = data, silent = silent, ...)
    }

    i <- i + 1
  }
  return(res)
}

#' Fit Brownian motion mixture with Gaussian jumps.
#' @param data Data to fit (increments).
#' @param mesh_size Time difference of data.
#' @param thresholds Jump thresholds (Defaults to `NA`, no filtering).
#' @param jump_quantile Quantile above which a data vector is a jump.
#' @return List of fitted Generalised Hyperbolic distribution.
#' @examples
#' n <- 10000
#' delta_time <- 0.01
#' data <- cbind(ghyp::rghyp(n, ghyp::ghyp()), ghyp::rghyp(n, ghyp::ghyp()))
#' fit_bm_compound_poisson(data = data, mesh_size = delta_time)
#' @importFrom stats optim
#' @export
fit_bm_compound_poisson <- function(data, mesh_size,
                                    thresholds = NA, jump_quantile = 0.9) {
  # data are increments
  n_data <- nrow(data)
  data_filtered_whole <- data_filtering(data, thresholds)
  filter_jumps <- data_filtered_whole$filter
  filter_data <- data_filtered_whole$data
  data_l2 <- apply(data, 1, function(x) {
    sqrt(sum(x * x))
  })
  jump_cutoff <- quantile(data_l2, jump_quantile)

  # Realised variance - cts
  rm_filter_data <- filter_data
  rm_filter_data[data_l2 > jump_cutoff, ] <- NA
  rv_c <- cov(rm_filter_data, use = "complete.obs") / mesh_size

  if (sum(filter_jumps) > 0) {
    jumps_n <- apply(filter_jumps, 2, sum)

    likelihood_poisson <- function(n_jump) {
      function(log_lambda) {
        loglik <- n_jump * log_lambda - exp(log_lambda) * mesh_size * n_data
        return(-loglik)
      }
    }
    jumps_fn_optim <- lapply(jumps_n, likelihood_poisson)
    poisson_intensities <- lapply(jumps_fn_optim, function(fn_optim) {
      exp(optim(par = c(-1), fn = fn_optim, method = "BFGS")$par)
    })
    poisson_intensities <- do.call(c, poisson_intensities)

    # synchronised / unique jumps
    n_unique_jumps <- sum(data_l2 > jump_cutoff)
    poisson_unique <- exp(
      optim(
        par = c(-1),
        fn = likelihood_poisson(n_unique_jumps),
        method = "BFGS"
      )$par
    )

    # Realised variance - total
    rv_total <- cov(data) / mesh_size

    # Realised variance - discontinuous
    rv_d <- as.matrix(rv_total - rv_c)
    jump_sigma <- t(t(rv_d) / poisson_unique)
    jump_sigma <- as.matrix(Matrix::nearPD(jump_sigma, corr = FALSE)$mat)

    # Fit ghyp only on jumps
    ghyp <- fit_ghyp_diffusion(data[abs(data_l2) > jump_cutoff, ])
  } else {
    poisson_intensities <- NA
    poisson_unique <- NA
    jump_sigma <- NA
    jumps_n <- NA
    ghyp <- NA
  }
  return(
    list(
      "sigma" = rv_c, "jump_sigma" = jump_sigma, "n_jumps" = jumps_n,
      "poisson" = poisson_intensities, "poisson_unique" = poisson_unique,
      "ghyp" = ghyp
    )
  )
}

bi_power_variation <- function(data, mesh_size) {
  abs_diff_data <- abs(apply(data, 2, diff))
  n_data <- nrow(abs_diff_data)
  abs_diff_data_wo_first <- abs_diff_data[2:n_data, ]
  abs_diff_data_wo_last <- abs_diff_data[1:(n_data - 1), ]
  cross_prod <- t(abs_diff_data_wo_first) %*% abs_diff_data_wo_last
  cross_prod <- cross_prod / (mesh_size * n_data)
  return(pi / 2 * cross_prod)
}
