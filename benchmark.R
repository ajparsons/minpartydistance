#!/usr/bin/env Rscript

# Simple benchmark script to test performance
suppressMessages(library(minPartyDistance))

# Test with different numbers of parties
test_cases <- list(
  small = c(-1, 0, 1),
  medium = c(-2, -1, 0, 1, 2),
  large = seq(-5, 5, by = 0.5)
)

cat("Performance Benchmark Results:\n")
cat("==============================\n")

for (name in names(test_cases)) {
  parties <- test_cases[[name]]
  
  cat(sprintf("\nTest case: %s (%d parties)\n", name, length(parties)))
  cat("Parties:", paste(parties, collapse = ", "), "\n")
  
  # Time the function
  start_time <- Sys.time()
  result <- minPartyDistance::mean_min_abs_distance_midpoints(parties)
  end_time <- Sys.time()
  
  execution_time <- as.numeric(end_time - start_time, units = "secs")
  
  cat("Result:", round(result, 6), "\n")
  cat("Execution time:", round(execution_time * 1000, 3), "ms\n")
}

# Run multiple iterations for a more accurate benchmark
cat("\n\nRepeated benchmark (100 iterations):\n")
cat("=====================================\n")

parties <- c(-2, -1, 0, 1, 2)
n_iterations <- 100

start_time <- Sys.time()
for (i in 1:n_iterations) {
  result <- minPartyDistance::mean_min_abs_distance_midpoints(parties)
}
end_time <- Sys.time()

total_time <- as.numeric(end_time - start_time, units = "secs")
avg_time <- total_time / n_iterations

cat("Total time for", n_iterations, "iterations:", round(total_time * 1000, 3), "ms\n")
cat("Average time per call:", round(avg_time * 1000000, 3), "μs\n")
cat("Final result:", round(result, 6), "\n")
