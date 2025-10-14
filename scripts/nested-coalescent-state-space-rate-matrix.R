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
  state_space <- nested_state_space_mapper(n, b)
  num_states <- nrow(state_space)
  total_gene_sample <- ncol(state_space)
  
  rate_matrix <- matrix(0, nrow = num_states, ncol = num_states)
  
  for (i in seq_len(num_states)) {
    for (j in seq_len(num_states)) {
      if (i == j) next
      
      state_diff <- state_space[i, ] - state_space[j, ]
      gene_mass_change <- sum(state_diff * seq_len(total_gene_sample))
      total_count_change <- sum(state_diff)
      
      # Species merging event (Kingman at species level)
      if (gene_mass_change == 0 && total_count_change == -1) {
        # Identify which species are being merged
        decreasing_sizes <- which(state_diff < 0)
        increasing_sizes <- which(state_diff > 0)
        
        if (length(decreasing_sizes) == 2 && length(increasing_sizes) == 1) {
          # Check if this is a valid species merge: k + l = m
          if (decreasing_sizes[1] + decreasing_sizes[2] == increasing_sizes[1]) {
            # Rate is product of binomial coefficients
            k_count <- state_space[j, decreasing_sizes[1]]
            l_count <- state_space[j, decreasing_sizes[2]]
            rate_matrix[i, j] <- k_count * l_count
          }
        }
      }
      # Gene merging event (Kingman at gene level)  
      else if (gene_mass_change == -1 && total_count_change == 0) {
        # Identify which species size decreased and which increased
        decreased_size <- which(state_diff < 0)
        increased_size <- which(state_diff > 0)
        
        if (length(decreased_size) == 1 && length(increased_size) == 1 &&
            increased_size + 1 == decreased_size) {
          # Gene coalescence within species of size k: rate = choose(k, 2)
          k <- decreased_size
          species_count <- state_space[j, decreased_size]  # Use target state count
          rate_matrix[i, j] <- species_count * choose(k, 2)
        }
      }
    }
  }
  
  # Set diagonal elements
  for (i in seq_len(num_states)) {
    rate_matrix[i, i] <- -sum(rate_matrix[i, ])
  }
  
  rate_matrix
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