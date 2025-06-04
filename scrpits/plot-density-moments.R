library(parallel)
source("scrpits/phase-type-density-moments.R")
source("scrpits/NestedCoalescent/nested-coalescent-state-space-rate-matrix.R")

plot_tmrca_density <- function(rate_matrix, limit_interval) {
  f <- tmrca_density(rate_matrix)
  x_vals <- seq(0, limit_interval, length.out = 250)
  y_vals <- unlist(mclapply(x_vals, f, mc.cores = detectCores() - 1))

  plot(x_vals, y_vals, type = "l", col = "darkred", lwd = 2,
       xlab = "x", ylab = "f(x)", main = "density")
}

plot_tmrca_1st_2nd_moments <- function() {
  compute_moments <- function(n) {
    g <- tmrca_moments(nested_rate_matrix(n, n))
    c(g(1), g(2))
  }
  n_values <- seq(1, 5, by = 1)
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