#' Apply the Levy recovery on data with (fitted) adjacency/dynamics matrix.
#' @param adj Adjacency or dynamics matrix.
#' @param data Data to apply the Levy recovery (vector).
#' @param times Observation times (vector).
#' @param look_ahead Step increment.
#' @return List with recovered increments, the adjacency matrix and cumulative sum.
#' @examples
#' n <- 1000
#' d <- 2
#' times <- seq(n)
#' delta_time <- 0.01
#' noise <- cbind(ghyp::rghyp(n, ghyp::ghyp()), ghyp::rghyp(n, ghyp::ghyp()))
#' data <- ConstructPath(diag(d), noise=noise, y_init=rep(0, d), delta_time=delta_time)
#' LevyRecovery(adj=diag(d), data=data, times=times)
#' @export
LevyRecovery <- function(adj, data, times, look_ahead=1){
  # only for OU processes; t_{n} -> t_{n+look_ahead}
  assertthat::are_equal(length(times), nrow(data))

  diff_x <- apply(data, MARGIN = 2, FUN = function(x){diff(x, lag = look_ahead)})
  diff_times <- diff(times, lag =  1)
  integrated_x <- zoo::rollapply(cbind(data, times),
                                 width = look_ahead+1,
                                 by.column = FALSE,
                                 FUN = function(sub_integ){
                                   tt <- sub_integ[,ncol(sub_integ)]
                                   y_col_bind <- sub_integ[,-ncol(sub_integ)]
                                   return(apply(y_col_bind,
                                                MARGIN = 2,
                                                FUN = function(y){
                                                  pracma::trapz(x = tt,
                                                                y = y)}
                                   ))}
  ) # create integral component-wise

  # few dimension checks
  assertthat::are_equal(nrow(integrated_x), nrow(diff_x))
  assertthat::are_equal(nrow(data)-look_ahead, nrow(diff_x))
  assertthat::are_equal(dim(adj)[1], ncol(integrated_x))

  q_integrated_x <- apply(integrated_x, MARGIN = 1, function(x){return(adj %*% x)})
  if(any(vapply(class(adj), function(x){x == 'dgCMatrix'}, T))){
    q_integrated_x <- lapply(q_integrated_x, as.matrix) # getting rid of potential sparse matrix
    q_integrated_x <- t(as.matrix(as.data.frame(q_integrated_x)))
  }else{
    q_integrated_x <- t(q_integrated_x)
  }

  # TODO implement ghyp
  # TODO implement with m > 1
  assertthat::are_equal(dim(diff_x)[1], dim(q_integrated_x)[1])
  assertthat::are_equal(dim(diff_x)[2], dim(q_integrated_x)[2])

  recover <- diff_x + q_integrated_x

  recovery_results <- list()
  recovery_results[['increments']] <- recover
  recovery_results[['cumsum']] <- apply(recover, MARGIN = 2, FUN = cumsum)
  recovery_results[['adj']] <- adj
  return(recovery_results)
}

#' Fit a Generalised Hyperbolic distribution.
#' @param data Data to fit.
#' @param ghyp_names Ghyp distribution names (in `c('NIG', 'GAUSS', 'VG', 'T', 'FULL')`).
#' @return List of fitted Generalised Hyperbolic distribution.
#' @examples
#' n <- 500
#' data <- cbind(ghyp::rghyp(n, ghyp::ghyp()), ghyp::rghyp(n, ghyp::ghyp()))
#' FitGhypDiffusion(data=data)
#' @export
FitGhypDiffusion <- function(data, ghyp_names='FULL'){
  # ghyp_models <- c(ghyp::fit.NIGmv, ghyp::fit.gaussmv, ghyp::fit.VGmv, ghyp::fit.tmv, ghyp::fit.ghypmv)
  ghyp_model_list <- list(
    'NIG'=ghyp::fit.NIGmv,
    'GAUSS'=ghyp::fit.gaussmv,
    'VG'=ghyp::fit.VGmv,
    'T'=ghyp::fit.tmv,
    'FULL'=ghyp::fit.ghypmv
  )
  assertthat::assert_that(all(lapply(ghyp_names, function(x){x %in% c('NIG', 'GAUSS', 'VG', 'T', 'FULL')})))

  ghyp_models <- ghyp_model_list[ghyp_names]
  i <- 1
  res <- list()
  for(model in ghyp_models){
    res[[ghyp_names[i]]] <- model(data = data)
    i <- i + 1
  }
  return(res)
}

#' Fit Brownian motion mixture with Gaussian jumps.
#' @param data Data to fit.
#' @param mesh_size Time difference of data.
#' @param thresholds Jump thresholds (Defaults to `NA`, no filtering).
#' @return List of fitted Generalised Hyperbolic distribution.
#' @examples
#' n <- 10000
#' delta_time <- 0.01
#' data <- cbind(ghyp::rghyp(n, ghyp::ghyp()), ghyp::rghyp(n, ghyp::ghyp()))
#' FitBrownianMotionCompoundPoisson(data=data, mesh_size=delta_time)
#' @importFrom stats optim
#' @export
FitBrownianMotionCompoundPoisson <- function(data, mesh_size, thresholds=NA){
  data <- apply(data, 2, cumsum)
  # data are increments
  n_data <- nrow(data)
  delta_data <- apply(data, 2, diff)
  abs_delta_data <- abs(delta_data)

  if(any(is.na(thresholds))){
    # stop('update FitBrownianMotionCompoundPoisson with data_original')
    # No jump filtering
    filter_jumps <- c(T)
  }else{
    assertthat::are_equal(length(thresholds), ncol(data))
    filter_jumps <- t(apply(abs_delta_data, 1, '<=', thresholds))
  }

  if(sum(!filter_jumps) > 0){
    jumps_dt <- apply(filter_jumps, 2, function(x){diff(which(as.numeric(x) == 0))*mesh_size})
    jumps_n <- lapply(jumps_dt, length)
    jumps_fn_optim <- lapply(
      jumps_dt,
      function(times){
        function(lambda){
          prob_no_jump <- exp(-exp(lambda))
          return(-(sum(log(1.0-exp(-exp(lambda)*times)))+log(prob_no_jump)*(n_data*mesh_size-sum(times))))
        }
      }
    )
    poisson_intensities <- lapply(jumps_fn_optim, function(fn_optim){exp(optim(par = c(-1), fn = fn_optim, method = 'BFGS')$par)})
    lapply(jumps_fn_optim, function(fn_optim){print(optim(par = c(-1), fn = fn_optim, method = 'BFGS'))})
    poisson_intensities <- do.call(c, poisson_intensities)
    filtered_delta_data <- abs_delta_data * filter_jumps
    rv_c <- (t(filtered_delta_data) %*% filtered_delta_data) / (mesh_size * n_data)
    rv_total <- (t(abs_delta_data) %*% abs_delta_data) / (mesh_size * n_data)
    rv_d <- as.matrix(rv_total - rv_c)
    jump_sigma <- t(t(rv_d)/(poisson_intensities))

    jump_sigma <- as.matrix(Matrix::nearPD(jump_sigma, corr=FALSE)$mat)
    return(list(sigma=rv_c, poisson=poisson_intensities, jump_sigma=jump_sigma, n_jumps=jumps_n))
  }else{
    rv_c <- (t(abs_delta_data) %*% abs_delta_data) / (mesh_size * n_data)
    return(list(sigma=rv_c, poisson=NA, jump_sigma=NA, n_jumps=NA))
  }
}

BiPowerVariation <- function(data, mesh_size){
  abs_diff_data <- abs(apply(data, 2, diff))
  n_data <- nrow(abs_diff_data)
  abs_diff_data_wo_first <- abs_diff_data[2:n_data,]
  abs_diff_data_wo_last <- abs_diff_data[1:(n_data-1),]
  cross_prod <- t(abs_diff_data_wo_first) %*% abs_diff_data_wo_last / (mesh_size * n_data)
  return(pi/2*cross_prod)
}
