
get_reg_fn <- function(reg, mle = NA, gamma = NA) {
  assertthat::assert_that(reg %in% c("l1", "l2", "adaptive"))
  if (reg == "l1") {
    fn <- function(x) sum(abs(x))
  }
  if (reg == "l2") {
    fn <- function(x) sqrt(sum(x * x))
  }
  if (reg == "adaptive") {
    if (any(is.na(mle))) {
      stop("`mle` is required to compute the adaptive lasso.")
    }
    if (any(is.na(gamma))) {
      stop("`gamma` is required to compute the adaptive lasso.")
    }
    denominator <- (abs(mle) + .Machine$double.eps)^gamma
    fn <- function(x) sum(abs(x / denominator))
  }
  return(fn)
}

#' Finding `lambda` for regularisation.
#' @param times Times at which data is given
#' @param data Values to compute the MLE with.
#' @param thresholds Jump threshold values.
#' @param reg Type of penalty (`l1`, `l2` or `adaptive`).
#' @param div Batch size/divisor to avoid large memory allocation.
#' @param gamma Adaptive MLE scaling parameter.
#' @param use_scaling Brownian motion covariance matrix scaling
#'     in the likelihood.
#' @param cross Cross-validation or not (defaults to `FALSE`).
#' @return List of top-choice for lambda, the collection of lambdas test
#'     and corresponding log-likelihood values
#' @export
finding_lambda <- function(times, data, thresholds,
                           reg = "l1", div = 1e5, gamma = NA,
                           use_scaling = F, cross = F) {
  assertthat::assert_that(isFALSE(cross) | (cross <= 1 & cross >= 0))
  if (!isFALSE(cross)) {
    n_train <- round((1 - cross) * nrow(data))
    data_test <- data[-seq_len(n_train), ]
    times_test <- times[-seq_len(n_train)]
    data <- data[seq_len(n_train), ]
    times <- times[seq_len(n_train)]
  } else {
    data_test <- data
    times_test <- times
  }

  mle_value <- grou_mle(
    times = times, data = data, thresholds = thresholds,
    div = div, mode = "node", output = "vector"
  )
  lambdas <- 10^{
    seq(-4, 2, length.out = 8)
  }
  loglik_vals <- vapply(
    lambdas, function(lmbd) {
      reg_lmbd <- grou_regularisation(
        times = times, data = data, thresholds = thresholds, lambda = lmbd,
        reg = "l1", div = div, output = "vector",
        gamma = gamma, use_scaling = use_scaling
      )
      loglike_lmbd <- likelihood_fn(
        times = times_test, data = data_test, thresholds = thresholds,
        lambda = 0, reg = "l1", div = div, log = F
      )
      return(-loglike_lmbd(reg_lmbd))
    }, 1
  )
  return(list(
    "lambda_picked" = lambdas[which.max(loglik_vals)],
    "lambdas" = lambdas,
    "loglik" = loglik_vals
  ))
}

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
#' @param MLE MLE for `adaptive` regularisation.
#' @param log Log-scale for the likelihood or not (defaults to `FALSE`).
#' @return (Log)likelihood of the GrOU process with penalty.
#' @examples
#' n <- 1000
#' d <- 10
#' times <- seq(n)
#' delta_time <- 0.01
#' noise <- matrix(rnorm(n * d, sd = sqrt(delta_time)), ncol = d)
#' data <- construct_path(
#'   diag(d),
#'   noise = noise, y_init = rep(0, d), delta_time = delta_time
#' )
#' loglik <- likelihood_fn(
#'   times = times, data = data, lambda = 1, div = 1e2
#' )
#' loglik(diag(d))
#' @export
likelihood_fn <- function(times, data, thresholds, lambda = NA,
                          reg = "l1", div = 1e5, gamma = NA,
                          use_scaling = F, mle = NA, log = T) {
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
#' @param MLE MLE for `adaptive` regularisation.
#' @param log Log-scale for the likelihood or not (defaults to `FALSE`).
#' @return (Log)likelihood of the GrOU process with penalty.
#' @examples
#' n <- 1000
#' d <- 10
#' times <- seq(n)
#' delta_time <- 0.01
#' noise <- matrix(rnorm(n * d, sd = sqrt(delta_time)), ncol = d)
#' data <- construct_path(
#'   diag(d),
#'   noise = noise, y_init = rep(0, d), delta_time = delta_time
#' )
#' loglik <- grad_likelihood_fn(
#'   times = times, data = data, div = 1e2
#' )
#' loglik(diag(d))
#' @export
grad_likelihood_fn <- function(times, data, thresholds, div = 1e5,
                               use_scaling = F, log = T) {
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



#' Regularisation schemes for the GrOU process that implements
#' a Lasso, Ridge or Adaptive Lasso penalty.
#' @param times Times at which data is given
#' @param data Values to compute the MLE with.
#' @param thresholds Jump threshold values.
#' @param lambda Penalty parameter.
#' @param reg Type of penalty (`l1`, `l2` or `adaptive`).
#' @param div Batch size/divisor to avoid large memory allocation.
#' @param output Output type: either "vector"or "matrix".
#' @param gamma Adaptive MLE scaling parameter.
#' @param cut_off Sparsity proportion, defaults to `NA`.
#' @param use_scaling Brownian motion covariance matrix scaling
#'     in the likelihood.
#' @return Regularised dynamics matrix.
#' @examples
#' n <- 1000
#' d <- 10
#' times <- seq(n)
#' delta_time <- 0.01
#' noise <- matrix(rnorm(n * d, sd = sqrt(delta_time)), ncol = d)
#' data <- construct_path(
#'   diag(d),
#'   noise = noise, y_init = rep(0, d), delta_time = delta_time
#' )
#' grou_regularisation(times = times, data = data, lambda = 1, div = 1e2)
#' @importFrom stats cov quantile
#' @export
grou_regularisation <- function(times, data, thresholds = NA, lambda = NA,
                                reg = "l1", div = 1e5, output = "vector",
                                gamma = NA, cut_off = NA, use_scaling = F) {
  assertthat::assert_that(output %in% c("vector", "matrix"))
  assertthat::assert_that(any(c(is.na(lambda), lambda >= 0.0)))
  assertthat::assert_that(any(c(is.na(cut_off), cut_off >= 0.0, cut_off <= 1.)))
  n_nodes <- ncol(data)
  mle_value <- grou_mle(
    times = times, data = data, thresholds = thresholds,
    div = div, mode = "node", output = "vector"
  )

  loglik_fn <- likelihood_fn(
    times = times, data = data, thresholds = thresholds,
    lambda = lambda, reg = reg, use_scaling = use_scaling,
    gamma = gamma, mle = mle_value, log = F
  )

  reg_optim <- optim(par = mle_value, fn = loglik_fn, method = "BFGS")
  reg_vector <- reg_optim$par
  if (!is.na(cut_off)) {
    assertthat::assert_that(0 <= cut_off & cut_off <= 1)
    reg_vector[abs(reg_vector) < quantile(abs(reg_vector), cut_off)] <- 0
  }

  if (output == "vector") {
    return(reg_vector)
  } else {
    return(matrix(reg_vector, n_nodes, n_nodes))
  }
}
