library(partitions)


#' Given a integer n, returns all the possible
#' [artmetic partitions](\url{https://en.wikipedia.org/wiki/Integer_partition})
#' of n in a vector representation of the
#' [Yung diagram](\url{https://en.wikipedia.org/wiki/Young_tableau#Diagrams})
#'
#' @param n the size of the initial sample.
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
#' all possible states of the nestes Kingman model
#' starting with s species with n genes each. Each row of
#' the matrix represents a sample where the i-th coordinate
#' is the number speciec with i genes.
#'
#' @param n the size of the initial gene sample.
#' @param b the size of the initial species sample.
nested_state_space_mapper <- function(n, b) {
  valid_vectors <- list()
  for (k in 1:n * b) {
    a <- state_space_mapper(k)
    zero_cols <- matrix(0, nrow = NROW(a), ncol =   b * n - k)
    a_new <- cbind(a, zero_cols)
    valid_vectors <- append(valid_vectors, split(a_new, row(a_new)))
  }
  return(valid_vectors)
}

cardinality_state_space(n, b) {
  e = 0
   for (k in 1:n * b) {
    e <- P()
  }
  return(valid_vectors)
}