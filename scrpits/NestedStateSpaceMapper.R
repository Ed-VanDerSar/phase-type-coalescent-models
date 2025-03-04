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
#' all possible states of the nested Kingman model
#' starting with s species with n genes each. Each row of
#' the matrix represents a sample where the i-th coordinate
#' is the number speciec with i genes.
#'
#' @param n the size of the initial gene sample.
#' @param b the size of the initial species sample.
#' @return the state space.
nested_state_space_mapper <- function(n, b) {
  valid_vectors <- list()
  for (k in 1:(n * b)) {
    a <- state_space_mapper(k)
    zero_cols <- matrix(0, nrow = NROW(a), ncol =   (b * n) - k)
    a_new <- cbind(a, zero_cols)
    valid_vectors <- append(valid_vectors, split(a_new, row(a_new)))
  }
  valid_vectors <- matrix(
                          unlist(valid_vectors),
                          nrow = length(valid_vectors),
                          byrow = TRUE)
  return(valid_vectors)
}

rate_matrix <- function(e) {
  dim <- NROW(e)
  total_gene_sample <- NCOL(e)
  rate <- matrix(0, ncol = dim, nrow = dim)
  for (i in 2:dim) {
    for (j in 1:(i - 1)) {
      ## establishing differences between two states
      c <- e[i, ] - e[j, ]
      ## Identifyung if the two states are compatible
      sum1 <- c %*% rep(1:total_gene_sample)
      gene_mass <- sum1[1, 1]
      ## check1==0 means that the size of the new blocks 
      ## equals the size of the disappearing blocks
      ## Identifying how many new blocks are created
      w1 <- ifelse(c > 0, 1, 0)
      sum2 <- c %*% w1
      created_species <- sum2[1, 1]
      ## check2==1 means that only one new block is created
      ## Disappearing blocks
      w2 <- ifelse(c < 0, 1, 0)
      neg_merged_species <- -c * w2
      ##Fullfilling the rate matrix
      if (gene_mass == 0 && created_species == 1) {
        provrate <- 1
        for (k in 1:total_gene_sample){
          provrate <- provrate * choose(e[j, k], neg_merged_species[k])
        }
        rate[j, i] <- provrate
      } else if (gene_mass == -1 && sum(c) == -1) {
        provrate <- 0
        for (k in 1:total_gene_sample){
          provrate <- provrate + choose(k, 2)
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