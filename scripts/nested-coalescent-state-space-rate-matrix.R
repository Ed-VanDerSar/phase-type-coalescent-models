library(partitions)

#' Generate all possible integer partitions in Young diagram vector representation
#'
#' @param n the integer to partition (must be positive)
#' @return matrix where each row represents a partition as (a1,...,an),
#'         where ai is the number of parts of size i, ordered by largest part first
#' @examples
#' all_possible_integer_partitions(4)
all_possible_integer_partitions <- function(n) {
  if (n <= 0) {
    stop("n must be a positive integer")
  }
  
  # Get all partitions of n
  partition_list <- parts(n)
  dim_partitions <- ncol(partition_list)
  
  # Initialize result matrix
  result_matrix <- matrix(0L, nrow = dim_partitions, ncol = n)
  
  # Convert partitions to count representation
  for (i in seq_len(dim_partitions)) {
    partition <- partition_list[, i]
    # Count occurrences of each part size
    for (part_size in partition) {
      if (part_size <= n) {
        result_matrix[i, part_size] <- result_matrix[i, part_size] + 1L
      }
    }
  }
  
  # Order by largest part first (column 1), then by next largest, etc.
  if (dim_partitions > 1) {
    # Create ordering key: convert to string representation for stable sorting
    order_keys <- apply(result_matrix, 1, function(row) {
      paste(sprintf("%03d", row), collapse = "")
    })
    result_matrix <- result_matrix[order(order_keys, decreasing = TRUE), , drop = FALSE]
  }
  
  result_matrix
}

#' Generate state space for nested Kingman model
#'
#' @param n the size of initial gene sample in each species
#' @param b the size of initial species sample  
#' @return matrix where each row represents a state, with i-th coordinate being
#'         the number of species with i genes, ordered by total gene mass
#' @examples
#' nested_state_space_mapper(2, 2)
nested_state_space_mapper <- function(n, b) {
  if (n <= 0 || b <= 0) {
    stop("n and b must be positive integers")
  } else if (n == 1 && b == 1) {
     return(matrix(1))
  }
  
  total_genes <- n * b
  valid_states <- list()
  
  # Process all possible total gene masses from 1 to total_genes-1
  for (k in 1:(total_genes - 1)) {
    partitions <- all_possible_integer_partitions(k)
    # Pad with zeros to make all states have length total_genes
    padded_partitions <- cbind(
      partitions,
      matrix(0L, nrow = nrow(partitions), ncol = total_genes - ncol(partitions))
    )
    valid_states <- c(valid_states, lapply(seq_len(nrow(padded_partitions)), 
                                           function(i) padded_partitions[i, ]))
  }
  
  # Add the full mass partitions (n*b)
  full_partitions <- all_possible_integer_partitions(total_genes)
  # Ensure full_partitions has correct dimensions
  if (ncol(full_partitions) < total_genes) {
    full_partitions <- cbind(
      full_partitions,
      matrix(0L, nrow = nrow(full_partitions), ncol = total_genes - ncol(full_partitions))
    )
  }
  
  # Find and reorder initial state (b species with n genes each)
  initial_state <- integer(total_genes)
  initial_state[n] <- b
  
  # Identify which row corresponds to the initial state
  is_initial_state <- apply(full_partitions, 1, function(row) {
    all(row == initial_state)
  })
  
  if (!any(is_initial_state)) {
    stop("Initial state not found in partitions")
  }
  
  # Reorder: put initial state last
  initial_state_index <- which(is_initial_state)
  reordered_full <- rbind(
    full_partitions[-initial_state_index, , drop = FALSE],
    full_partitions[initial_state_index, , drop = FALSE]
  )
  
  # Add all full partitions to valid states
  valid_states <- c(valid_states, lapply(seq_len(nrow(reordered_full)), 
                                         function(i) reordered_full[i, ]))
  
  # Convert to matrix and reverse order (so initial state is first in time)
  state_matrix <- do.call(rbind, valid_states)
  state_matrix[rev(seq_len(nrow(state_matrix))), , drop = FALSE]
}

#' Construct rate matrix for nested Kingman coalescent model
#'
#' @param n the size of initial gene sample in each species
#' @param b the size of initial species sample
#' @return the rate matrix Q where Q[i,j] is transition rate from state i to j
#' @examples
#' nested_rate_matrix(2, 1)
nested_rate_matrix <- function(n, b) {
  e <- nested_state_space_mapper(n, b)
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

# Example usage and verification
if (FALSE) {
  # Test the functions
  partitions_4 <- all_possible_integer_partitions(4)
  print("Partitions of 4:")
  print(partitions_4)
  
  state_space_2_2 <- nested_state_space_mapper(2, 2)
  print("State space for n=2, b=2:")
  print(state_space_2_2)
  
  rate_matrix_2_1 <- nested_rate_matrix(2, 1)
  print("Rate matrix for n=2, b=1:")
  print(rate_matrix_2_1)
}