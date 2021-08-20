#' Computes the psi-GrOU MLE of a batch of data
#' @param times Times at which data is given
#' @param data Values to compute the MLE with.
#' @param output Output type: either "vector"or "matrix".
#' @return the psi-GrOU (node-wise) MLE
#' @importFrom stats rnorm runif sd stl ts
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
#' node_mle(times, data)
#' @export
node_mle <- function(times, data, output = "vector") {
  assertthat::assert_that(
    output %in% c("vector", "matrix"),
    msg = paste('output should be "node" or "network", given', output)
  )
  assertthat::assert_that(
    assertthat::are_equal(length(times), nrow(data))
  )

  n_nodes <- ncol(data)
  n_data <- nrow(data)
  delta_t <- diff(times)
  delta_data <- apply(data, MARGIN = 2, diff)
  wo_last <- data[-n_data, ]
  expanded_delta_data <- t(rep(1, n_nodes)) %x% delta_data
  expanded_wo_last <- wo_last %x% t(rep(1, n_nodes))

  numerator <- colSums(expanded_wo_last * expanded_delta_data)
  denominator <- t(wo_last * delta_t) %*% wo_last
  sd_denom <- sd(denominator)
  inv_denominator <- solve(denominator / sd_denom)
  inv_denominator <- inv_denominator / sd_denom
  mle <- -inv_denominator %x% diag(n_nodes) %*% numerator
  if (output == "vector") {
    return(mle)
  } else {
    return(matrix(mle, nrow = n_nodes, ncol = n_nodes, byrow = FALSE))
  }
}

#' Computes the MLE numerator and denominator of a batch of data
#' @param times Times at which data is given
#' @param data Values to compute the MLE with.
#' @param thresholds Jump threshold values. Default is `NA` (no filtering).
#' @return A list of two matrices (numerator and denominator).
core_node_mle <- function(times, data, thresholds = NA) {
  if (is.character(thresholds)) {
    stop(paste(
      "`thresholds` should be NA or numeric. Received:", class(thresholds)
    ))
  }

  n_nodes <- ncol(data)
  n_data <- nrow(data)
  delta_t <- diff(times)
  delta_data <- apply(data, MARGIN = 2, diff)
  wo_last <- data[-n_data, ]

  kron_delta_data <- lapply(
    seq_len(ncol(wo_last)),
    function(i) {
      # transpose for the filtering
      tmp <- t(delta_data * as.vector(wo_last[, i]))
      # filtering
      if (!any(is.na(thresholds))) {
        filters <- apply(abs(delta_data), 1, "<=", rep(thresholds[i], n_nodes))
        filters <- t(filters)
        tmp <- tmp * as.vector(filters)
      }
      return(colSums(t(tmp)))
    }
  )
  numerator <- t(do.call(cbind, kron_delta_data))
  denominator <- t(wo_last * as.vector(delta_t)) %*% wo_last

  return(list(numerator = numerator, denominator = denominator))
}

#' Returns the components (i.e. numerator and denominator) of the MLE
#' @param times Times at which data is given
#' @param data Values to compute the MLE with.
#' @param thresholds Jump threshold values.
#' @param div Batch size/divisor to avoid large memory allocation.
#' @param output String to indicate
#'     the form of the output: `vector` or `matrix`.
#' @return Computes the GrOU MLE components.
node_mle_components <- function(times, data, thresholds, div = 1e5,
                                output = "vector") {
  assertthat::assert_that(
    output %in% c("vector", "matrix"),
    msg = paste('output should be "node" or "network", given', output)
  )
  assertthat::assert_that(
    assertthat::are_equal(length(times), nrow(data))
  )

  n_nodes <- ncol(data)
  n_data <- nrow(data)

  if (div < n_data) {
    idx <- seq(1, n_data + div - 1, by = div)
    idx <- idx[idx <= n_data]
    if (!(n_data %in% idx)) {
      idx <- c(idx, n_data)
    }
  } else {
    idx <- c(1, n_data + 1)
  }
  collection_num_denom <- lapply(
    seq_len(length(idx) - 1),
    FUN = function(i) {
      idx_couple <- idx[i]:(idx[i + 1] - 1)
      core_node_mle(
        times = times[idx_couple], data = data[idx_couple, ],
        thresholds = thresholds
      )
    }
  )

  numerator <- Reduce("+", lapply(collection_num_denom, function(coll) {
    coll$numerator
  }))
  denominator <- Reduce("+", lapply(collection_num_denom, function(coll) {
    coll$denominator
  }))

  # post-processing
  numerator <- matrix(numerator, nrow = n_nodes, byrow = FALSE)
  sd_numerator <- sd(numerator)
  sd_denominator <- sd(denominator)

  return(list(
    numerator = numerator, denominator = denominator,
    sd_numerator = sd_numerator, sd_denominator = sd_denominator
  ))
}

#' Constructs the MLE from its components
#' @param components List of components necessary to implement GrOU MLE.
#' @param output String to indicate
#'     the form of the output: `vector` or `matrix`.
#' @return GrOU MLE made from `components`.
make_node_mle <- function(components, output) {
  # unpacking components
  numerator <- components$numerator
  sd_numerator <- components$sd_numerator
  denominator <- components$denominator
  sd_denominator <- components$sd_denominator

  # Numerator
  numerator <- numerator / sd_numerator

  # Denominator
  inv_denominator <- solve(denominator / sd_denominator)
  inv_denominator <- inv_denominator / sd_denominator

  mle <- apply(numerator, MARGIN = 2, function(x) {
    -inv_denominator %*% x
  }) * sd_numerator

  if (output == "vector") {
    return(c(mle))
  } else {
    return(mle)
  }
}

#' Constructs the Node MLE with jump thresholding
#' @param times Times at which data is given
#' @param data Values to compute the MLE with.
#' @param thresholds Jump threshold values.
#' @param div Batch size/divisor to avoid large memory allocation.
#' @param output String to indicate
#'     the form of the output: `vector` or `matrix`.
#' @return The GrOU MLE in matrix/vector form in psi parametrisation.
node_mle_long <- function(times, data, thresholds, div = 1e5,
                          output = "vector") {
  components <- node_mle_components(times, data, thresholds, div, output)
  return(make_node_mle(components, output))
}

#' Computes the GrOUs MLE of a batch of data
#' @param times Times at which data is given
#' @param data Values to compute the MLE with.
#' @param adj Adjacency matrix of the underlying network.
#' @param thresholds Jump threshold values.
#' @param div Batch size/divisor to avoid large memory allocation.
#' @param mode GrOU mode: either "node" (psi-GrOU) or "network" (theta-GrOU).
#' @param output Output type: either "vector"or "matrix".
#' @return The GrOU MLE in matrix/vector form in theta or psi parametrisation.
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
#' grou_mle(times, data, adj = diag(d), div = 1e2)
#' @export
grou_mle <- function(times, data, adj = NA, thresholds = NA,
                     div = 1e3, mode = "node", output = "vector") {
  assertthat::assert_that(
    mode %in% c("node", "network"),
    msg = paste('mode should be "node" or "network", given', mode)
  )
  assertthat::assert_that(
    length(times) == nrow(data),
    msg = "`times` and `data` should have compatible length and shape.s"
  )

  n_nodes <- ncol(data)
  if (any(is.na(adj))) {
    adj_normalised <- row_normalised(matrix(1, n_nodes, n_nodes))
  } else {
    adj_normalised <- row_normalised(adj)
  }
  diag(adj_normalised) <- 0.0 # to be sure

  node_mle_value <- node_mle_long(
    times, data,
    thresholds = thresholds, div = div, output = "vector"
  )
  node_mle_value <- as.vector(t(matrix(node_mle_value, nrow = n_nodes)))

  # post-processing with mode, output and adj
  if (mode == "node") {
    adj_normalised[adj_normalised != 0] <- 1
    if (output == "vector") {
      d_a <- c(diag(n_nodes) + adj_normalised)
      return(d_a * node_mle_value)
    } else {
      if (output == "matrix") {
        adj_full <- diag(n_nodes) + adj_normalised
        return(adj_full * matrix(node_mle_value, n_nodes, n_nodes, byrow = F))
      }
    }
  } else {
    a_l2 <- sum(adj_normalised)
    adj_ones <- adj_normalised
    adj_ones[adj_ones != 0] <- 1
    diag(adj_ones) <- 0

    if (a_l2 < .Machine$double.eps) {
      theta_1 <- 0.0
    } else {
      theta_1 <- (as.vector(adj_ones) %*% node_mle_value) / a_l2
    }

    theta_2 <- (c(diag(n_nodes)) %*% node_mle_value) / n_nodes
    if (output == "vector") {
      return(c(theta_1, theta_2))
    } else {
      if (output == "matrix") {
        return(c(theta_1) * adj_normalised + c(theta_2) * diag(n_nodes))
      }
    }
  }
}
