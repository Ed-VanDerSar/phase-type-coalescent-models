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
#' @param m the size of the initial species sample.
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

rate_matrix <- function(n, m, a1, a2, b1, b2) {
  e <- two_islands_state_space(n, m)
  dim <- NROW(e)
  rate <- matrix(0, ncol = dim, nrow = dim)
  for (i in 2:dim) {
    for (j in 1:dim) {
      ## establishing differences between two states
      c <- e[i, ] - e[j, ]
      ## Identifying if the two states are compatible
      blocks_transformed <- c[1] + c[2]
      # blocks_transformed==0 means that we have a potential migaration event
      # blocks_transformed==-1 means that we have a merging event

      ##Fullfilling the rate matrix
      if (blocks_transformed == 0) {
        if (c[1] == -1) {
          rate[j, i] <-  (e[j, 1] * b1)
        } else if (c[1] == 1) {
          rate[j, i] <-  (e[j, 2] * b2)
        }
      } else if (blocks_transformed == -1) {
        if (c[1] == -1) {
          rate[j, i] <-  (choose(e[j, 1], 2) * a1)
        } else if (c[2] == -1) {
          rate[j, i] <-  (choose(e[j, 2], 2) * a2)
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

library("expm")

tmrca_moments <- function(n, m, power) {

RmRw <- rate_matrix(n, m, 1, 1, 1, 1)
## Restrict the rate matrix and invert it
inv_rate <- solve(-RmRw[1:(ncol(RmRw) - 1), 1:(ncol(RmRw) - 1)])
## Obtain the kth moment of the branch length
id <- diag(1, (ncol(RmRw) - 1))
e <- rep(1, ncol(RmRw) - 1)
dr <- diag(RmRw$Reward)
moment <- (inv_rate %*% dr)%^%power
moment <- id[1, ] %*% moment %*% e
return(moment)
}