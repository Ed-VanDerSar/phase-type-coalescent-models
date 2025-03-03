library(partitions)

#' Given integers n and b, returna the cardinality of the stae space
#' of  of the nested Kingman model starting with s species with n genes each.
#'
#' @param n the size of the initial gene sample.
#' @param b the size of the initial species sample.
cardinality_state_space <- function(n, b) {
  e <- 0
  for (k in 1:(n * b)) {
    e <- e + P(k)
  }
  return(e)
}

#' Given integers n and b, checks that the rows of the
#' state space matrix are in fact admissible states for the
#' nested coalescent model.
#'
#' @param n the size of the initial gene sample.
#' @param b the size of the initial species sample.
#' @param e the state space matrix.
consistency_state_space <- function(n, b, e) {
  for (i in 1:nrow(e)) {
    state <- e[i, ]
    total_gene_mass <- sum((1:length(state)) * state)
    if (total_gene_mass > (n * b)) {
      print(paste("State indexed with", i, "is not admissible"))
      break
    }
  }
}
