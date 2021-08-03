#' Fasen's multivariate OU least squares regression
#' \doi{/10.1016/j.jeconom.2012.08.019}.
#' @param data Multivariate ata to perform the estimation on.
#' @return A \code{(d, d)} least square matrix.
#' @examples
#' n <- 100
#' d <- 10
#' delta_time <- 0.01
#' noise <- matrix(rnorm(n * d), ncol = d)
#' path <- construct_path(diag(10), noise, d, delta_time)
#' fasen_regression(path)
#' @export
fasen_regression <- function(data) {
  data <- as.matrix(data)
  data_without_first <- data[-1, ]
  data_without_last <- data[-nrow(data), ]
  sub <- t(data_without_last) %*% data_without_last
  return(t(data_without_first) %*% data_without_last %*% solve(sub))
}
