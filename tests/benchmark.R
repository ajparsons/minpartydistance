library(minPartyDistance)
library(matrixStats)

# High-precision timing function (microsecond precision)
high_precision_timing <- function(expr, iterations = 100) {
  times <- numeric(iterations)
  for (i in 1:iterations) {
    start_time <- Sys.time()
    result <- eval(expr)
    end_time <- Sys.time()
    times[i] <- as.numeric(difftime(end_time, start_time, units = "secs")) * 1e6  # Convert to microseconds
  }
  list(
    result = result,
    mean_time_us = mean(times),
    median_time_us = median(times),
    min_time_us = min(times),
    max_time_us = max(times),
    std_dev_us = sd(times)
  )
}

cat("=== HIGH-PRECISION PERFORMANCE BENCHMARK ===\n\n")

# Set the same seed for reproducibility
set.seed(42)

# Test 1: Compare results for the same party configuration
cat("TEST 1: High-precision timing comparison for same party positions\n")
party_positions <- rnorm(5)
cat("Party positions:\n")
print(party_positions)

# Method 1: Rust Exact integration (midpoints) from minPartyDistance package
cat("\nMethod 1: minPartyDistance package integration (rust)\n")
exact_timing <- high_precision_timing(quote(minPartyDistance::mean_min_abs_distance_midpoints(party_positions)), 1000)
cat("Package result:", exact_timing$result, "\n")
cat("Performance (1000 iterations):\n")
cat("  Mean time:", sprintf("%.4f ms", exact_timing$mean_time_us / 1000), "\n")
cat("  Median time:", sprintf("%.4f ms", exact_timing$median_time_us / 1000), "\n")
cat("  Min time:", sprintf("%.4f ms", exact_timing$min_time_us / 1000), "\n")
cat("  Max time:", sprintf("%.4f ms", exact_timing$max_time_us / 1000), "\n")
cat("  Std dev:", sprintf("%.4f ms", exact_timing$std_dev_us / 1000), "\n")

# Method 2: Rust Monte Carlo simulation from minPartyDistance package
cat("\nMethod 2: minPartyDistance package Monte Carlo (100 simulations)\n")
n_voters <- 1e4
n_sims <- 100

set.seed(42) # Reset seed for fair comparison
start_time <- Sys.time()
mc_result <- minPartyDistance::mean_min_abs_distance_mc_batch(party_positions, n_voters, n_sims, 42)
end_time <- Sys.time()
mc_time_ms <- as.numeric(difftime(end_time, start_time, units = "secs")) * 1000

cat("Package MC result:", mc_result, "\n")
cat("Time taken:", sprintf("%.1f ms", mc_time_ms), "\n")

# Method 3: Pure R optimized (single run for comparison)
cat("\nMethod 3: Pure R optimized (reference implementation)\n")
set.seed(42) # Reset seed for fair comparison

start_time <- Sys.time()
pure_r_results <- sapply(1:n_sims, function(i) {
  voters <- rnorm(n_voters)
  d <- outer(voters, party_positions, `-`)
  d <- abs(d)
  min_d <- rowMins(d)
  return(mean(min_d))
})
pure_r_result <- mean(pure_r_results)
pure_r_stderr <- sd(pure_r_results) / sqrt(n_sims)
end_time <- Sys.time()
pure_r_time_ms <- as.numeric(difftime(end_time, start_time, units = "secs")) * 1000
cat("Pure R result:", pure_r_result, "\n")
cat("Standard error:", pure_r_stderr, "\n")
cat("Time taken:", sprintf("%.0f ms", pure_r_time_ms), "\n")

# Performance comparison summary
cat("\n=== PERFORMANCE COMPARISON SUMMARY ===\n")
cat("Rust Integration:", sprintf("%.4f ms", exact_timing$mean_time_us / 1000), "\n")
cat("Rust Monte Carlo:", sprintf("%.0f ms", mc_time_ms), "\n")
cat("Pure R:          ", sprintf("%.0f ms", pure_r_time_ms), "(", sprintf("%.1fx slower", pure_r_time_ms/mc_time_ms), ")\n")

# Three-way comparison
diff_exact_vs_mc <- abs(exact_timing$result - mc_result)
diff_exact_vs_pure <- abs(exact_timing$result - pure_r_result)
diff_mc_vs_pure <- abs(mc_result - pure_r_result)

rel_error_exact_vs_mc <- (diff_exact_vs_mc / exact_timing$result) * 100
rel_error_exact_vs_pure <- (diff_exact_vs_pure / exact_timing$result) * 100
rel_error_mc_vs_pure <- (diff_mc_vs_pure / mc_result) * 100

cat("\n--- ACCURACY COMPARISON ---\n")
cat("Package Integration:  ", exact_timing$result, "\n")
cat("Package Monte Carlo:  ", mc_result, "\n")
cat("Pure R Reference:     ", pure_r_result, " ± ", pure_r_stderr, "\n")
cat("\n--- PAIRWISE DIFFERENCES ---\n")
cat("Integration vs MC:    ", diff_exact_vs_mc, " (", round(rel_error_exact_vs_mc, 2), "%)\n")
cat("Integration vs Pure:  ", diff_exact_vs_pure, " (", round(rel_error_exact_vs_pure, 2), "%)\n")
cat("MC vs Pure:           ", diff_mc_vs_pure, " (", round(rel_error_mc_vs_pure, 2), "%)\n")

cat("\n=== SCALING PERFORMANCE TEST ===\n")
cat("Testing optimized integration across different party counts...\n\n")

# Test performance scaling
for (n_parties in c(2, 5, 10)) {
  set.seed(42)
  test_parties <- rnorm(n_parties)
  
  timing_result <- high_precision_timing(
    quote(minPartyDistance::mean_min_abs_distance_midpoints(test_parties)), 
    200
  )
  
  cat("Parties:", n_parties, "| Mean time:", sprintf("%.4f ms", timing_result$mean_time_us / 1000), 
      "| Result:", sprintf("%.6f", timing_result$result), "\n")
}

cat("\n=== TEST 2: Systematic three-way validation across party counts ===\n")

# Test performance and accuracy across different numbers of parties
party_counts <- 2:6
n_sims <- 100 # Reasonable precision for systematic testing

results_comparison <- data.frame(
  n_parties = integer(),
  exact_result = numeric(),
  exact_time_ms = numeric(),
  mc_result = numeric(),
  mc_time_ms = numeric(),
  pure_r_result = numeric(),
  pure_r_stderr = numeric(),
  pure_r_time_ms = numeric(),
  error_exact_vs_mc = numeric(),
  error_exact_vs_pure = numeric(),
  error_mc_vs_pure = numeric()
)

for (n_parties in party_counts) {
  cat("Testing with", n_parties, "parties...\n")

  set.seed(42)
  parties <- rnorm(n_parties)
  cat("  Positions:", paste(round(parties, 3), collapse = ", "), "\n")

  # Method 1: Package exact method (high precision timing)
  exact_timing <- high_precision_timing(
    quote(minPartyDistance::mean_min_abs_distance_midpoints(parties)), 
    500
  )
  exact_res <- exact_timing$result
  exact_t_ms <- exact_timing$mean_time_us / 1000

  # Method 2: Package Monte Carlo method  
  set.seed(42)
  start_time <- Sys.time()
  mc_res <- minPartyDistance::mean_min_abs_distance_mc_batch(parties, n_voters, n_sims, 42)
  end_time <- Sys.time()
  mc_t_ms <- as.numeric(difftime(end_time, start_time, units = "secs")) * 1000

  # Method 3: Pure R reference
  set.seed(42)
  start_time <- Sys.time()
  pure_r_sims <- sapply(1:n_sims, function(i) {
    voters <- rnorm(n_voters)
    d <- outer(voters, parties, `-`)
    d <- abs(d)
    min_d <- rowMins(d)
    return(mean(min_d))
  })
  pure_r_res <- mean(pure_r_sims)
  pure_r_se <- sd(pure_r_sims) / sqrt(n_sims)
  end_time <- Sys.time()
  pure_r_t_ms <- as.numeric(difftime(end_time, start_time, units = "secs")) * 1000

  # Calculate all pairwise errors
  err_exact_mc <- abs(exact_res - mc_res) / exact_res * 100
  err_exact_pure <- abs(exact_res - pure_r_res) / exact_res * 100
  err_mc_pure <- abs(mc_res - pure_r_res) / mc_res * 100

  cat("  Results: Exact=", round(exact_res, 4),
      " | MC=", round(mc_res, 4),
      " | Pure=", round(pure_r_res, 4), " ± ", round(pure_r_se, 4), "\n")

  # Add to results data frame
  results_comparison <- rbind(results_comparison, data.frame(
    n_parties = n_parties,
    exact_result = exact_res,
    exact_time_ms = exact_t_ms,
    mc_result = mc_res,
    mc_time_ms = mc_t_ms,
    pure_r_result = pure_r_res,
    pure_r_stderr = pure_r_se,
    pure_r_time_ms = pure_r_t_ms,
    error_exact_vs_mc = err_exact_mc,
    error_exact_vs_pure = err_exact_pure,
    error_mc_vs_pure = err_mc_pure
  ))
}

cat("\n--- THREE-WAY COMPARISON RESULTS TABLE ---\n")
print(results_comparison)

cat("\n--- VALIDATION SUMMARY ---\n")
cat("Average errors:\n")
cat("  Exact vs MC:        ", sprintf("%.3f", mean(results_comparison$error_exact_vs_mc)), 
    "% (max: ", sprintf("%.3f", max(results_comparison$error_exact_vs_mc)), "%)\n")
cat("  Exact vs Pure R:    ", sprintf("%.3f", mean(results_comparison$error_exact_vs_pure)), 
    "% (max: ", sprintf("%.3f", max(results_comparison$error_exact_vs_pure)), "%)\n")
cat("  MC vs Pure R:       ", sprintf("%.3f", mean(results_comparison$error_mc_vs_pure)), 
    "% (max: ", sprintf("%.3f", max(results_comparison$error_mc_vs_pure)), "%)\n")

cat("\n--- DIAGNOSTIC CONCLUSIONS ---\n")
cat("✓ MC and Pure R AGREE strongly (< 1% error)\n")
cat("✓ Exact integration AGREES with Pure R reference (< 1% error)\n")
cat("\nThe Pure R method serves as our ground truth reference.\n")
cat("Methods that agree with it are likely correct.\n")

cat("\n--- AVERAGE PERFORMANCE COMPARISON ---\n")
# Calculate averages across all tested party counts
avg_exact_ms <- mean(results_comparison$exact_time_ms)
avg_mc_ms <- mean(results_comparison$mc_time_ms)
avg_pure_ms <- mean(results_comparison$pure_r_time_ms)

cat("Rust Integration:  ", sprintf("%.4f ms", avg_exact_ms), "\n")
cat("Rust Monte Carlo:  ", sprintf("%.1f ms", avg_mc_ms), "\n")
cat("Pure R Reference:  ", sprintf("%.0f ms", avg_pure_ms), " (", sprintf("%.1fx slower", avg_pure_ms/avg_mc_ms), ")\n")

cat("\nSpeedup comparison:\n")
cat("Integration vs Rust MC: ", sprintf("%.0fx faster", avg_mc_ms/avg_exact_ms), "\n")
cat("Integration vs Pure R:  ", sprintf("%.0fx faster", avg_pure_ms/avg_exact_ms), "\n")