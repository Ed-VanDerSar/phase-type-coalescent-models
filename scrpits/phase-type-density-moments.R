library("expm")

#' Probability density of the
#' time to absorption associated with the given rate matrix. We are assuming the
#' underlying markov chain starts in the first state and is absorbed in
#' in the last state.
#'
#' @param rate_matrix the rate matrix.
#' @return the density function.
tmrca_density <- function(rate_matrix) {
  ## Restrict the rate matrix and invert it
  rest_rate <- rate_matrix[1:(ncol(rate_matrix) - 1),
                           1:(ncol(rate_matrix) - 1)]
  id_matrix <- diag(1, (ncol(rate_matrix) - 1))
  e <- rep(1, ncol(rate_matrix) - 1)
  exit_rate <- - rest_rate %*% e
  function(x) {
    id_matrix[1, ] %*% (expm(rest_rate * x)) %*% exit_rate
  }
}

#' Function that receives an integer n
#' and computes the n-th moment of the
#' time to abortion associated with the given rate matrix. We are assuming the
#' underlying markov chain starts in the first state and is absorbed in
#' in the last state.
#'
#' @param rate_matrix the rate matrix.
#' @return the moment function.
tmrca_moments <- function(rate_matrix) {
  ## Restrict the rate matrix and invert it
  inv_rate <- solve(-rate_matrix[1:(ncol(rate_matrix) - 1),
                                 1:(ncol(rate_matrix) - 1)])
  id_matrix <- diag(1, (ncol(rate_matrix) - 1))
  e <- rep(1, ncol(rate_matrix) - 1)
  ## Obtain the nth moment of the tree height
  function(power) {
    id_matrix[1, ] %*% ((inv_rate) %^% power) %*% e
  }
}