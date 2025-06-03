library(partitions)

#' Given a integer n, returns all the possible
#' [artmetic partitions](\url{https://en.wikipedia.org/wiki/Integer_partition})
#' of n in a vector representation of the
#' [Yung diagram](\url{https://en.wikipedia.org/wiki/Young_tableau#Diagrams})
#'
#' @param n the size of the initial sample.
#' @return all the possible partitions
state_space_mapper <- function(n) {
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
  r_matrix <- r_matrix[order(r_matrix[, 1], decreasing = TRUE), ]
  return(r_matrix)
}

#' Given integers n and b, returna a matrix encoding
#' all possible truncated_ of the nested Kingman model
#' starting with s species with n genes each. Each row of
#' the matrix represents a sample where the i-th coordinate
#' is the number species with i genes.
#'
#' @param n the size of the initial gene sample in each species.
#' @param b the size of the initial species sample.
#' @return the state space.
nested_state_space_mapper <- function(n, b) {
  valid_states <- list()
  # Process and add all remaining states by its gene mass
  for (k in 1:((n * b) - 1)) {
    truncated_states <- state_space_mapper(k)
    zero_cols <- matrix(0, nrow = NROW(truncated_states), ncol =   (b * n) - k)
    states <- cbind(truncated_states, zero_cols)
    valid_states <- append(valid_states, split(states, row(states)))
  }
  ## Reorder to choose the first state as the one with b species with n genes
  initial_states <- state_space_mapper((n * b))
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
  ordered_valid_states <-  valid_states[rev(seq_len(nrow(valid_states))), ]
  return(ordered_valid_states)
}

#' Populates the rate matrix for the state space associated to
#' the nested coalescente model. By no, we only use the rates for a 
#' Kingman-in-Kingman model.
#'
#' @param n the size of the initial gene sample in each species.
#' @param b the size of the initial species sample.
#' @return the rate matrix.
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
  return(rate)
}

library("expm")

#' Evaluates the poroibability density...
#' 
#' @param x the point where we evaluate the density
#' @param rate_matrix the rate matrix.
#' @return the value of the density on x.
tmrca_density <- function(x, rate_matrix) {
  ## Restrict the rate matrix and invert it
  rest_rate <- rate_matrix[1:(ncol(rate_matrix) - 1),
                           1:(ncol(rate_matrix) - 1)]
  rest_rate <-  rest_rate * x
  value <- expm(rest_rate)
  id <- diag(1, (ncol(rate_matrix) - 1))
  e <- rep(1, ncol(rate_matrix) - 1)
  exit_rate <- - rest_rate %*% e
  value <- id[1, ] %*% value %*% exit_rate
  return(value)
}


tmrca_moments <- function(rate_matrix, power) {
  ## Restrict the rate matrix and invert it
  inv_rate <- solve(-rate_matrix[1:(ncol(rate_matrix) - 1),
                                 1:(ncol(rate_matrix) - 1)])
  ## Obtain the kth moment of the tree hight
  id <- diag(1, (ncol(rate_matrix) - 1))
  e <- rep(1, ncol(rate_matrix) - 1)
  moment <- (inv_rate) %^% power
  moment <- id[1, ] %*% moment %*% e
  return(moment)
}

library(parallel)

plot_tmrca_density <- function(rate_matrix, limit_interval) {
  f <- function(x) {
    value <- tmrca_density(x, rate_matrix)
    return(value)
  }
  x_vals <- seq(0, limit_interval, length.out = 250)
  y_vals <- unlist(mclapply(x_vals, f, mc.cores = detectCores() - 1))

  plot(x_vals, y_vals, type = "l", col = "darkred", lwd = 2,
       xlab = "x", ylab = "f(x)", main = "density")
}

plot_tmrca_1st_2nd_moments <- function() {
  compute_moments <- function(n) {
    rate_matrix <- nested_rate_matrix(n, n)
    c(tmrca_moments(rate_matrix, 1), tmrca_moments(rate_matrix, 2))
  }
  n_values <- seq(1, 10, by = 1)
  results <- mclapply(n_values, compute_moments, mc.cores = detectCores() - 1)
  moments1 <- sapply(results, `[[`, 1)
  moments2 <- sapply(results, `[[`, 2)

  # Plot the first moment
  plot(n_values, moments1, type = "l", col = "blue", lwd = 2,
       ylim = range(c(moments1, moments2)),
       xlab = "Sample size", ylab = "TMRCA Moment",
       main = "TMRCA Moments vs sample size")

  # Add the second moment
  lines(n_values, moments2, col = "red", lwd = 2)

  # Add a legend
  legend("topleft", legend = c("Moment 1", "Moment 2"),
         col = c("blue", "red"), lwd = 2)
}