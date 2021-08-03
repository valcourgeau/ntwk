# This script attempts to provide construction functions for 4 main types of networks;
# They cover the following properties
#   - max(degree) = 0, 2, 4/8, d-1

#' Checks that \code{|theta_1| < theta_2}.
#' @param theta_1 Off-diagonal value or network effect.
#' @param theta_2 Diagonal value or momentum effect.
CheckThetas <- function(theta_1, theta_2){
  assertthat::assert_that(theta_2 > abs(theta_1), msg = paste('In', match.call()[1], 'must verify theta_2 > abs(theta_1)'))
}

#' Creates adjacency matrix with max degree zewro.
#' @param d Number of nodes.
#' @param theta_1 Off-diagonal value or network effect.
#' @param theta_2 Diagonal value or momentum effect.
#' @return adjacency matrix of a graph with max degree zero.
#' @example IsolatedNetwork(d=10)
#' @export
IsolatedNetwork <- function(d, theta_1, theta_2){
  # theta_2 is the diagonal element
  return(theta_2*diag(d))
}

#' Creates adjacency matrix with min degree one and max degree two.
#' @param d Number of nodes.
#' @param theta_1 Off-diagonal value or network effect.
#' @param theta_2 Diagonal value or momentum effect.
#' @param directed Boolean if lower triangular is the opposite to the upper triangular matrix.
#' @return Adjacency matrix of a graph with max degree two and minimum one.
#' @example PolymerNetwork(d=10)
#' @export
PolymerNetwork <- function(d, theta_1, theta_2, directed=F){
  CheckThetas(theta_1, theta_2)
  # theta_2 is the diagonal element
  mat_temp <- AugmentedDiag(d = d, offset = 1) + if(directed){-AugmentedDiag(d = d, offset = 1)} else{AugmentedDiag(d=d,offset=-1)}
  mat_temp <- theta_1 * mat_temp / rowSums(abs(mat_temp))
  return(mat_temp + IsolatedNetwork(d = d, theta_2 = theta_2))
}

#' Lattice/grid-like network with max degree of four across all nodes.
#' @param d Number of nodes.
#' @param theta_1 Off-diagonal value or network effect.
#' @param theta_2 Diagonal value or momentum effect.
#' @param directed Boolean if lower triangular is the opposite to the upper triangular matrix.
#' @example LatticeNetwork(10, theta_1 = 1, theta_2 = 2)
#' @importFrom copCAR adjacency.matrix
#' @export
LatticeNetwork <- function(d, theta_1, theta_2, directed=F){
  # TODO add directed for it
  CheckThetas(theta_1, theta_2)

  if(directed) stop('NotImplementedError')
  d_real <- round(sqrt(d))
  net_matrix <- adjacency.matrix(d_real)
  if(d_real^2 < d){
    net_matrix <- rbind(net_matrix, matrix(0, nrow=d-d_real^2, ncol=d_real*d_real))
    net_matrix <- cbind(net_matrix, rbind(matrix(0, ncol=d-d_real^2, nrow=d_real*d_real), diag(1, nrow = d - d_real^2)))
  }

  net_matrix <- theta_1 * RowNormalised(net_matrix)
  diag(net_matrix) <- theta_2

  return(net_matrix)
}

#' Fully-connected graph.
#' @param d Number of nodes.
#' @param theta_1 Off-diagonal value or network effect.
#' @param theta_2 Diagonal value or momentum effect.
#' @param directed Boolean if lower triangular is the opposite to the upper triangular matrix.
#' @example FullyConnectedNetwork(10, 1, 2)
#' @export
FullyConnectedNetwork <- function(d, theta_1, theta_2, directed=F){
  net_temp <- matrix(1, d, d)
  net_temp <- theta_1 * RowNormalised(net_temp)

  if(directed){
    net_temp[lower.tri(net_temp)] <- -net_temp[lower.tri(net_temp)]
  }

  diag(net_temp) <- theta_2
  return(net_temp)
}
