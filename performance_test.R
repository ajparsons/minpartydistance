#!/usr/bin/env Rscript

# Comprehensive performance test
suppressMessages(library(minPartyDistance))

cat("Advanced Performance Analysis\n")
cat("=============================\n\n")

# Test cases with different party configurations
test_cases <- list(
  single = c(0),
  two_parties = c(-1, 1),
  small = c(-1, 0, 1),
  medium = c(-2, -1, 0, 1, 2),
  large = seq(-5, 5, by = 0.5),
  very_large = seq(-10, 10, by = 0.2),
  asymmetric = c(-3, -1, 0, 1, 5),
  clustered = c(-2, -1.9, -1.8, 1.8, 1.9, 2)
)

# Function to benchmark with multiple iterations
benchmark_function <- function(parties, iterations = 1000) {
  times <- numeric(iterations)
  
  for (i in 1:iterations) {
    start_time <- Sys.time()
    result <- minPartyDistance::mean_min_abs_distance_midpoints(parties)
    end_time <- Sys.time()
    times[i] <- as.numeric(end_time - start_time, units = "secs") * 1e6  # microseconds
  }
  
  list(
    result = result,
    mean_time = mean(times),
    median_time = median(times),
    min_time = min(times),
    max_time = max(times),
    sd_time = sd(times)
  )
}

# Run benchmarks
for (name in names(test_cases)) {
  parties <- test_cases[[name]]
  cat(sprintf("Test case: %s (%d parties)\n", name, length(parties)))
  cat("Parties:", paste(round(parties, 2), collapse = ", "), "\n")
  
  benchmark_result <- benchmark_function(parties, 1000)
  
  cat("Result:", round(benchmark_result$result, 6), "\n")
  cat("Performance (1000 iterations):\n")
  cat("  Mean time:", round(benchmark_result$mean_time, 2), "μs\n")
  cat("  Median time:", round(benchmark_result$median_time, 2), "μs\n")
  cat("  Min time:", round(benchmark_result$min_time, 2), "μs\n")
  cat("  Max time:", round(benchmark_result$max_time, 2), "μs\n")
  cat("  Std dev:", round(benchmark_result$sd_time, 2), "μs\n")
  cat("\n")
}

# Scaling test
cat("Scaling Analysis\n")
cat("================\n")
party_counts <- c(2, 5, 10, 20, 50, 100)

for (n in party_counts) {
  parties <- seq(-n/2, n/2, length.out = n)
  benchmark_result <- benchmark_function(parties, 100)
  cat(sprintf("Parties: %d, Mean time: %.2f μs, Result: %.6f\n", 
              n, benchmark_result$mean_time, benchmark_result$result))
}
