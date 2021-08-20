
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
