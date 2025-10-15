library(partitions)

#' Given a integer n, returns all the possible
#' [arithmetic partitions](\url{https://en.wikipedia.org/wiki/Integer_partition})
#' of n in a vector representation of the
#' [Yung diagram](\url{https://en.wikipedia.org/wiki/Young_tableau#Diagrams})
#'
#' @param n the size of the initial sample.
#' @return all the possible partitions
all_possible_integer_partitions <- function(n) {
  dim <- P(n)
  ## Definition of the state matrix
  r_matrix <- matrix(ncol = n, nrow = dim)
  ## Set of partitions of [n]
  x <- parts(n)
  ## Rewriting the partitions as (a1,...,an)
  for (i in 1:dim) {
    y <- x[, dim - i + 1]
    for (j in 1:n) {
      r_matrix[i, j] <- length(which(y == j))
    }
  }
  ## Reordering
  r_matrix[order(r_matrix[, 1], decreasing = TRUE), ]
}

#' Given integers n and b, return a matrix encoding
#' all possible truncated_ of the nested Kingman model
#' starting with s species with n genes each. Each row of
#' the matrix represents a sample where the i-th coordinate
#' is the number species with i genes.
#'
#' @param n the size of the initial gene sample in each species.
#' @param b the size of the initial species sample.
#' @return the state space.
nested_state_space <- function(n, b) {
  valid_states <- list()
  # Process states all states with total gene mass less than n*b
  for (k in 1:((n * b) - 1)) {
    truncated_states <- all_possible_integer_partitions(k)
    zero_cols <- matrix(0, nrow = NROW(truncated_states), ncol =   (b * n) - k)
    states <- cbind(truncated_states, zero_cols)
    valid_states <- append(valid_states, split(states, row(states)))
  }
  ## Reorder to choose the first state as the one with b species with n genes
  initial_states <- all_possible_integer_partitions((n * b))
  # Find the index of the row with the maximum value in column i
  max_index <- which.max(initial_states[, n])
  # Reorder the matrix: put all rows except max_index first,
  # then add max_index row at the end
  reordered_initial_states <- rbind(
                                    initial_states[-max_index, , drop = FALSE],
                                    initial_states[max_index, , drop = FALSE])
  valid_states <- append(
                         valid_states,
                         split(
                               reordered_initial_states,
                               row(reordered_initial_states)))
  valid_states <- matrix(
                         unlist(valid_states),
                         nrow = length(valid_states),
                         byrow = TRUE)
  # reorder states
  valid_states[rev(seq_len(nrow(valid_states))), ]
}

#' Populates the rate matrix for the state space associated to
#' the nested coalescent model. By no, we only use the rates for a
#' Kingman-in-Kingman model.
#'
#' @param n the size of the initial gene sample in each species.
#' @param b the size of the initial species sample.
#' @return the rate matrix.
nested_rate_matrix <- function(n, b) {
  e <- nested_state_space(n, b)
  dim <- NROW(e)
  total_gene_sample <- NCOL(e)
  rate <- matrix(0, ncol = dim, nrow = dim)
  for (i in 1:dim) {
    for (j in 1:dim) {
      ## establishing differences between two states
      c <- e[i, ] - e[j, ]
      ## Identifying if the two states are compatible
      gene_mass <- c %*% rep(1:total_gene_sample)
      ## gene_mass==0 means that the size of the new blocks
      ## equals the size of the disappearing blocks, so
      ## this means SPECIES MERGING EVENT
      ## gene_mass==-1 means GENE MERGING EVENT
      ## Disappearing blocks
      w2 <- ifelse(c < 0, 1, 0)
      neg_merged_species <- -c * w2
      ##Fullfilling the rate matrix
      if (gene_mass == 0 && sum(c) == -1) {
        provrate <- 1
        for (k in 1:total_gene_sample){
          provrate <- provrate * choose(e[j, k], neg_merged_species[k])
        }
        rate[j, i] <- provrate
      } else if (gene_mass == -1 && sum(c) == 0) {
        provrate <- 0
        for (k in 1:total_gene_sample){
          provrate <- provrate + choose(k * w2[k], 2)
        }
        rate[j, i] <- provrate
      }
    }
  }
  ## Diagonal part of the matrix
  for (i in 1:dim){
    rate[i, i] <- - sum(rate[i, ])
  }
  rate
}