library(partitions)

#' Get all ordered partitions of a positive integer into two positive integers
#'
#' This function returns a matrix containing all ordered pairs of positive integers
#' that sum to the given integer `k`. Each row in the matrix is one such partition.
#'
#' @param k A positive integer greater than 1.
#'
#' @return A matrix with `k - 1` rows and 2 columns. Each row is a pair of positive integers (a, b)
#'         such that a + b = k.
#'
#' @examples
#' get_ordered_partitions_matrix(5)
#' #      [,1] [,2]
#' # [1,]    1    4
#' # [2,]    2    3
#' # [3,]    3    2
#' # [4,]    4    1
#'
get_ordered_partitions_matrix <- function(k) {
  stopifnot(k > 1)
  matrix(c(rep(1:(k - 1), each = 1), (k - 1):1), ncol = 2)
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
nested_state_space_mapper <- function(n, m) {

}