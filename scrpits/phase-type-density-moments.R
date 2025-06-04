library("expm")

#' Retrieves the probability density of the
#' time to absortion asssociated with the given rate matrix. We are assuming the
#' underliying markov chain starts in the first state and is aborbed in
#' in the last state.
#'
#' @param rate_matrix the rate matrix.
#' @return the density funtion.
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

#' Retrieves a funtion that receives an integer n
#' and computes the n-th moment of the
#' time to absortion asssociated with the given rate matrix. We are assuming the
#' underliying markov chain starts in the first state and is aborbed in
#' in the last state.
#'
#' @param rate_matrix the rate matrix.
#' @return the moment funtion.
tmrca_moments <- function(rate_matrix) {
  ## Restrict the rate matrix and invert it
  inv_rate <- solve(-rate_matrix[1:(ncol(rate_matrix) - 1),
                                 1:(ncol(rate_matrix) - 1)])
  id_matrix <- diag(1, (ncol(rate_matrix) - 1))
  e <- rep(1, ncol(rate_matrix) - 1)
  ## Obtain the nth moment of the tree hight
  function(power) {
    id_matrix[1, ] %*% ((inv_rate) %^% power) %*% e
  }
}