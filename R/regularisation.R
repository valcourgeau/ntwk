
get_reg_fn <- function(reg, mle = NA, gamma = NA) {
  assertthat::assert_that(reg %in% c("l1", "l2", "adaptive"))
  if (reg == "l1") {
    fn <- function(x) sum(abs(x))
  }
  if (reg == "l2") {
    fn <- function(x) sum(x * x)
  }
  if (reg == "adaptive") {
    if (any(is.na(mle))) {
      stop("`mle` is required to compute the adaptive lasso.")
    }
    if (any(is.na(gamma))) {
      stop("`mle` is required to compute the adaptive lasso.")
    }
    denominator <- (abs(mle) + .Machine$double.eps)^gamma
    fn <- function(x) sum(abs(x / denominator))
  }
  return(fn)
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
#' @importFrom stats cov
#' @export
grou_regularisation <- function(times, data, thresholds = NA, lambda = NA,
                                reg = "l1", div = 1e5, output = "vector",
                                gamma = NA, cut_off = NA, use_scaling = F) {
  # TODO(val) add warning for lambda
  # TODO(val) add scaling
  # TODO(val) add filtering on levy increments
  assertthat::assert_that(output %in% c("vector", "matrix"))
  assertthat::assert_that(any(c(is.na(lambda), lambda >= 0.0)))
  assertthat::assert_that(any(c(is.na(cut_off), cut_off >= 0.0, cut_off <= 1.)))


  n_nodes <- ncol(data)
  mle_value <- grou_mle(
    times = times, data = data, thresholds = thresholds,
    div = div, mode = "node", output = "vector"
  )

  if (use_scaling) {
    levy_increments <- levy_recovery(
      adj = matrix(mle_value, n_nodes, n_nodes),
      times = times,
      data = data,
      look_ahead = 1
    )$increments
    scaling_matrix <- cov(levy_increments)
    # covariance matrix to use in scaling
    inv_scaling_matrix <- solve(scaling_matrix)
  } else {
    inv_scaling_matrix <- Matrix::Diagonal(n_nodes)
  }

  reg_fn <- get_reg_fn(reg, mle_value, gamma)
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
        %*% kronecker(inv_scaling_matrix, t(components$denominator))
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

  reg_optim <- optim(par = mle_value, fn = fn_optim, method = "BFGS")
  reg_vector <- reg_optim$par
  if (!is.na(cut_off)) {
    reg_vector[abs(reg_vector) < quantile(abs(reg_vector), cut_off)] <- 0
  }

  if (output == "vector") {
    return(reg_vector)
  } else {
    return(matrix(reg_vector, n_nodes, n_nodes))
  }
}
