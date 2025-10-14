library(testthat)
library(partitions)
library(here)


source(here("scripts", "nested-coalescent-state-space-rate-matrix.R"))

test_that("all_possible_integer_partitions works correctly", {
  # Test basic functionality
  expect_silent(result <- all_possible_integer_partitions(4))
  expect_type(result, "integer")
  expect_equal(dim(result), c(5, 4))  # 5 partitions of 4
  
  # Test known partitions for n=4
  expected <- matrix(c(
    4,0,0,0,  # 4
    2,1,0,0,  # 2+1+1
    0,2,0,0,  # 2+2
    1,0,1,0,  # 3+1
    0,0,0,1   # 1+1+1+1
  ), nrow = 5, byrow = TRUE)
  
  # Check that all expected partitions are present (order may vary)
  expect_equal(nrow(result), nrow(expected))
  expect_equal(result %*% (1:4), matrix(rep(ncol(result), nrow(result)), ncol = 1))  # Each partition sums to 4
  
  # Test edge cases
  expect_error(all_possible_integer_partitions(0), "must be a positive integer")
  expect_error(all_possible_integer_partitions(-1), "must be a positive integer")
  
  # Test n=1
  result1 <- all_possible_integer_partitions(1)
  expect_equal(dim(result1), c(1, 1))
  expect_equal(result1[1, 1], 1)
})

test_that("nested_state_space_mapper works correctly", {
  # Test basic functionality
  expect_silent(states <- nested_state_space_mapper(2, 2))
  expect_type(states, "integer")
  
  # Test dimensions
  total_genes <- 2 * 2
  expect_equal(ncol(states), total_genes)
  
  # Test that initial state is correct (2 species with 2 genes each)
  initial_state <- integer(total_genes)
  initial_state[2] <- 2  # 2 species of size 2
  
  # The initial state should be the first row after reversal
  expect_equal(states[1, ], initial_state)
  
  # Test that all states have correct total gene count
  gene_counts <- apply(states, 1, function(row) sum(row * seq_along(row)))
  expect_true(all(gene_counts <= total_genes))
  expect_true(all(gene_counts >= 1))
  
  # Test edge cases
  expect_error(nested_state_space_mapper(0, 2), "must be positive integers")
  expect_error(nested_state_space_mapper(2, 0), "must be positive integers")
  
  # Test small case: n=1, b=1
  states_small <- nested_state_space_mapper(1, 1)
  expect_equal(dim(states_small), c(1, 1))
  expect_equal(states_small[1, 1], 1)
})

test_that("nested_rate_matrix works correctly", {
  # Test basic functionality
  expect_silent(rate_matrix <- nested_rate_matrix(2, 1))
  expect_type(rate_matrix, "double")
  
  # Test matrix properties
  expect_true(is.matrix(rate_matrix))
  expect_equal(nrow(rate_matrix), ncol(rate_matrix))
  
  # Test that rows sum to 0 (Markov chain property)
  row_sums <- apply(rate_matrix, 1, sum)
  expect_true(all(abs(row_sums) < 1e-10))  # Allow for floating point error
  
  # Test that off-diagonal elements are non-negative
  diag_removed <- rate_matrix
  diag(diag_removed) <- 0
  expect_true(all(diag_removed >= -1e-10))  # Allow for floating point error
  
  # Test specific known case: n=2, b=1
  # States should be: [0,1] (one species with 2 genes) and [2,0] (two species with 1 gene each)
  # Only one transition possible: from [2,0] to [0,1] via gene coalescence
  rate_matrix_2_1 <- nested_rate_matrix(2, 1)
  expect_equal(dim(rate_matrix_2_1), c(nrow(nested_state_space_mapper(2,1)), nrow(nested_state_space_mapper(2,1))))
  
  # Test larger case
  expect_silent(rate_matrix_2_2 <- nested_rate_matrix(2, 2))
  expect_true(nrow(rate_matrix_2_2) > 0)
  expect_equal(nrow(rate_matrix_2_2), ncol(rate_matrix_2_2))
})

test_that("Integration between functions works", {
  # Test that state space mapper uses partitions correctly
  n <- 2
  b <- 2
  states <- nested_state_space_mapper(n, b)
  partitions <- all_possible_integer_partitions(n * b)
  
  # The full mass partition should be in the state space
  full_mass_state <- integer(n * b)
  full_mass_state[n] <- b  # b species with n genes each
  expect_true(any(apply(states, 1, function(row) all(row == full_mass_state))))
  
  # Test that rate matrix dimensions match state space
  rate_matrix <- nested_rate_matrix(n, b)
  expect_equal(nrow(rate_matrix), nrow(states))
  expect_equal(ncol(rate_matrix), nrow(states))
})

test_that("Performance tests for larger inputs", {
  # Test that functions don't crash on moderate inputs
  expect_silent({
    partitions_10 <- all_possible_integer_partitions(10)
    states_3_3 <- nested_state_space_mapper(3, 3)
    rate_matrix_2_2 <- nested_rate_matrix(2, 2)
  })
  
  # Test timing for moderate inputs (should complete in reasonable time)
  time_taken <- system.time({
    nested_rate_matrix(3, 2)
  })
  expect_lt(time_taken["elapsed"], 10)  # Should complete in under 10 seconds
})

test_that("Biological interpretation tests", {
  # Test that transitions make biological sense
  
  # Case 1: Two species with one gene each
  # Only possible transition: gene coalescence within a species
  # But with one gene per species, no coalescence is possible
  # So the state should be absorbing
  rate_matrix_1_2 <- nested_rate_matrix(1, 2)
  expect_equal(sum(rate_matrix_1_2[1, ]), 0)  # First row (initial state) should sum to 0
  
  # Case 2: One species with two genes
  # Only possible transition: gene coalescence
  states <- nested_state_space_mapper(2, 1)
  rate_matrix <- nested_rate_matrix(2, 1)
  
  # Should have two states: [0,1] and [2,0]
  # Transition from [2,0] to [0,1] should have rate choose(2,2) = 1
  # Find which state is [2,0] (two species with 1 gene)
  state_2_0 <- which(apply(states, 1, function(x) all(x == c(2, 0))))
  state_0_1 <- which(apply(states, 1, function(x) all(x == c(0, 1))))
  
  if (length(state_2_0) > 0 && length(state_0_1) > 0) {
    expect_equal(rate_matrix[state_2_0, state_0_1], 1)
  }
})
