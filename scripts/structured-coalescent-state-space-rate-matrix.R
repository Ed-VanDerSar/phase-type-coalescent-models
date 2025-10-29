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

#' Given integers n and m, return a matrix encoding
#' all possible states for a coalescent with migration model structured in two
#' independent islands starting with n and m lineages respectively.
#' Each row of the output matrix represents represents a sample
#' where the first coordinate registers the number of blocks in the first
#' island and the second the number of blocks in the second one.
#'
#' The states (1,0) and (0,1) represents the final state of the system where
#' only one lineage remains.
#'
#' @param n the size of the initial gene sample in first island
#' @param m the size of the initial gene sample in second island
#' @return the state space.
two_islands_state_space <- function(n, m) {
  valid_states <- list()
  for (i in 1:(n + m - 1)) {
    states <- get_ordered_partitions(i)
    valid_states <- append(valid_states, split(states, row(states)))
  }
  ## Reorder to choose the first state as (n,m)
  initial_states <- get_ordered_partitions((n + m))
  # Find the index of the row where the first column is n
  start_index <- which(initial_states[, 1] == n)
  # Reorder the matrix
  initial_states <- rbind(
                          initial_states[-start_index, ],
                          initial_states[start_index, ])
  valid_states <- append(
                         valid_states,
                         split(
                               initial_states,
                               row(initial_states)))
  valid_states <- matrix(
                         unlist(valid_states),
                         nrow = length(valid_states),
                         byrow = TRUE)
  # reorder rows to star with (n,m) and end in the absorbing state
  ordered_valid_states <-  valid_states[rev(seq_len(nrow(valid_states))), ]
  return(ordered_valid_states)
}

#' Provided migration rate and merging rates for each island, return the rate 
#' matrix for a coalescent model with migration structured in two
#' independent islands starting with n and m lineages respectively.
#'
#' The states (1,0) and (0,1) represents the final state of the system where
#' only one lineage remains.
#'
#' @param n the size of the initial gene sample in first island
#' @param m the size of the initial gene sample in second island
#' @param mergingRateFirstIsland 
#' @param mergingRateSecondIsland
#' @param migrationRateFirstToSecond
#' @param migrationRateSecondToFirst
#' @return the corresponding rate matrix.
two_islands_rate_matrix <- function(n, m, 
                                    mergingRateFirstIsland, 
                                    mergingRateSecondIsland,
                                    migrationRateFirstToSecond, 
                                    migrationRateSecondToFirst) {
  e <- two_islands_state_space(n, m)
  dim <- NROW(e)
  rate <- matrix(0, ncol = dim, nrow = dim)
  for (i in 2:dim) {
    for (j in 1:dim) {
      ## establishing differences between two states
      c <- e[i, ] - e[j, ]
      ## Identifying if the two states are compatible
      blocks_transformed <- c[1] + c[2]
      # blocks_transformed==0 means that we have a potential migration event
      # blocks_transformed==-1 means that we have a merging event

      ##Fulfilling the rate matrix
      if (blocks_transformed == 0) {
        if (c[1] == -1) {
          rate[j, i] <-  (e[j, 1] * migrationRateFirstToSecond)
        } else if (c[1] == 1) {
          rate[j, i] <-  (e[j, 2] * migrationRateSecondToFirst)
        }
      } else if (blocks_transformed == -1) {
        if (c[1] == -1) {
          rate[j, i] <-  (choose(e[j, 1], 2) * mergingRateFirstIsland)
        } else if (c[2] == -1) {
          rate[j, i] <-  (choose(e[j, 2], 2) * mergingRateSecondIsland)
        }
      }
    }
  }
  ## Diagonal part of the matrix
  for (i in 1:dim){
    rate[i, i] <- - sum(rate[i, ])
  }
  # Step 1: Add the last column to the one before it
  rate[, ncol(rate) - 1] <- rate[, ncol(rate) - 1] + rate[, ncol(rate)]
  # Step 2: Remove the last column
  rate <- rate[, -ncol(rate)]
  # Step 3: Remove the last row
  rate <- rate[-nrow(rate), ]
  return(rate)
}