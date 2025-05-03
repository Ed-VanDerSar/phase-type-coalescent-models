#' Get all ordered partitions of a non-negative integer into two non-negative
#' integers. This function returns a matrix of all ordered pairs of non-negative
#' integers that sum to the given integer `k`.
#'
#' @param k A non-negative integer.
#'
#' @return A matrix with `k + 1` rows and 2 columns.Each row is a pair (a, b)
#'         such that a + b = k.
#'
#' @examples
#' get_ordered_partitions_with_zero(3)
#' #      [,1] [,2]
#' # [1,]    0    3
#' # [2,]    1    2
#' # [3,]    2    1
#' # [4,]    3    0
#'
get_ordered_partitions <- function(k) {
  stopifnot(is.numeric(k), k >= 0, length(k) == 1)
  matrix(c(0:k, k:0), ncol = 2)
}

#' Given integers n and b, returna a matrix encoding
#' all possible states for a coalescent with migration model structured in two
#' independent islans starting with n and m lienages respectively.
#' Each row of the outpu matix represents represents a sample
#' where the first coordinate registers the number of blocks in the first
#' island and the second the number of blocks in the second one.
#'
#' The states (1,0) and (0,1) represnts the final state of the sistem where
#' only one lineage remains.
#'
#' @param n the size of the initial gene sample in each species.
#' @param b the size of the initial species sample.
#' @return the state space.
state_space_mapper <- function(n, m) {
  valid_states <- list()
  for (i in 1:(n + m)) {
    states <- get_ordered_partitions(i)
    valid_states <- append(valid_states, split(states, row(states)))
  }
  valid_states <- matrix(
                         unlist(valid_states),
                         nrow = length(valid_states),
                         byrow = TRUE)
  # reorder states
  ordered_valid_states <-  valid_states[rev(seq_len(nrow(valid_states))), ]
  return(ordered_valid_states)
}