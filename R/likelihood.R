
#' Likelihood function for the GrOU process with penalty.
#' @param times Times at which data is given
#' @param data Values to compute the MLE with.
#' @param thresholds Jump threshold values.
#' @param lambda Penalty parameter (defaults to `NA` with no penalty).
#' @param reg Type of penalty (`l1`, `l2` or `adaptive`).
#' @param div Batch size/divisor to avoid large memory allocation.
#' @param gamma Adaptive MLE scaling parameter.
#' @param use_scaling Brownian motion covariance matrix scaling
#'     in the likelihood.
#' @param mle MLE for `adaptive` regularisation.
#' @param log Log-scale for the likelihood or not (defaults to `FALSE`).
#' @return (Log)likelihood of the GrOU process with penalty.
#' @examples
#' n <- 1000
#' d <- 10
#' times <- seq(n)
#' delta_time <- 0.01
#' beta_value <- 0.499
#' noise <- matrix(rnorm(n * d, sd = sqrt(delta_time)), ncol = d)
#' data <- construct_path(
#'   diag(d),
#'   noise = noise, y_init = rep(0, d), delta_time = delta_time
#' )
#' thresholds <- rep(delta_time^beta_value, d)
#' loglik <- likelihood_fn(
#'   times = times, data = data,
#'   thresholds = thresholds,
#'   lambda = 1, div = 1e2
#' )
#' loglik(diag(d))
#' @export
likelihood_fn <- function(times, data, thresholds, lambda = NA,
                          reg = "l1", div = 1e5, gamma = NA,
                          use_scaling = FALSE, mle = NA, log = TRUE) {
  # TODO(val) add warning for lambda
  # TODO(val) add scaling
  # TODO(val) add filtering on levy incrementss
  n_nodes <- ncol(data)

  if (use_scaling) {
    mle_value <- grou_mle(
      times = times, data = data, thresholds = thresholds,
      div = div, mode = "node", output = "vector"
    )
    levy_increments <- levy_recovery(
      adj = matrix(mle_value, n_nodes, n_nodes),
      times = times,
      data = data,
      look_ahead = 1
    )$increments
    scaling_matrix <- cov(levy_increments) / max(diff(times))
    # covariance matrix to use in scaling
    inv_scaling_matrix <- solve(scaling_matrix)
  } else {
    inv_scaling_matrix <- Matrix::Diagonal(n_nodes)
  }

  reg_fn <- get_reg_fn(reg, mle, gamma)
  components <- node_mle_components(
    times = times, data = data, thresholds = thresholds,
    div = div, output = "vector"
  )
  fn_optim <- function(adj_vector) {
    # Compute the penalised log-likelihood
    adj_vector <- c(adj_vector)
    adj <- matrix(adj_vector, n_nodes, n_nodes)
    to_maximise <- c(
      sum(
        adj_vector
        %*% kronecker(Matrix::Diagonal(n_nodes), inv_scaling_matrix)
          %*% as.vector(t(components$numerator))
      ),
      sum(
        adj_vector
        %*% kronecker(t(components$denominator), inv_scaling_matrix)
          %*% adj_vector
      )
    )

    # scaling of variance possible:  2 / n_nodes
    to_maximise[2] <- 0.5 * to_maximise[2]

    to_maximise <- -sum(to_maximise)
    if (!is.na(lambda) & abs(lambda) > .Machine$double.eps) {
      to_maximise <- to_maximise - lambda * reg_fn(adj)
    }
    return(-to_maximise)
  }
  return(fn_optim)
}

#' Gradient of the GrOU likelihood function with penalty.
#' @param times Times at which data is given
#' @param data Values to compute the MLE with.
#' @param thresholds Jump threshold values.
#' @param div Batch size/divisor to avoid large memory allocation.
#' @param use_scaling Brownian motion covariance matrix scaling
#'     in the likelihood.
#' @param log Log-scale for the likelihood or not (defaults to `FALSE`).
#' @return (Log)likelihood of the GrOU process with penalty.
#' @examples
#' n <- 1000
#' d <- 10
#' times <- seq(n)
#' delta_time <- 0.01
#' beta_value <- 0.499
#' noise <- matrix(rnorm(n * d, sd = sqrt(delta_time)), ncol = d)
#' data <- construct_path(
#'   diag(d),
#'   noise = noise, y_init = rep(0, d), delta_time = delta_time
#' )
#' thresholds <- rep(delta_time^beta_value, d)
#' grad_loglik <- grad_likelihood_fn(
#'   times = times, data = data, thresholds = thresholds, div = 1e2
#' )
#' grad_loglik(diag(d))
#' @export
grad_likelihood_fn <- function(times, data, thresholds, div = 1e5,
                               use_scaling = FALSE, log = TRUE) {
  n_nodes <- ncol(data)

  if (use_scaling) {
    mle_value <- grou_mle(
      times = times, data = data, thresholds = thresholds,
      div = div, mode = "node", output = "vector"
    )
    levy_increments <- levy_recovery(
      adj = matrix(mle_value, n_nodes, n_nodes),
      times = times,
      data = data,
      look_ahead = 1
    )$increments
    scaling_matrix <- cov(levy_increments) / max(diff(times))
    # covariance matrix to use in scaling
    inv_scaling_matrix <- solve(scaling_matrix)
  } else {
    inv_scaling_matrix <- Matrix::Diagonal(n_nodes)
  }

  components <- node_mle_components(
    times = times, data = data, thresholds = thresholds,
    div = div, output = "vector"
  )
  fn_optim <- function(adj_vector) {
    # Compute the gradient of the log-likelihood
    adj_vector <- as.vector(adj_vector)
    grad_val <- kronecker(
      Matrix::Diagonal(n_nodes), inv_scaling_matrix
    ) %*% as.vector(t(components$numerator))
    grad_val <- grad_val + kronecker(
      t(components$denominator), inv_scaling_matrix
    ) %*% adj_vector
    return(-grad_val)
  }
  return(fn_optim)
}
