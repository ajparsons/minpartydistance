#![allow(non_snake_case)]

use extendr_api::prelude::*;
use wide::f64x4;

// Pre-computed constants for standard normal PDF
const INV_SQRT_2PI: f64 = 0.3989422804014327; // 1.0 / sqrt(2.0 * PI)

// Fast standard normal probability density function calculation
#[inline]
fn fast_standard_normal_pdf(x: f64) -> f64 {
    INV_SQRT_2PI * (-0.5 * x * x).exp()
}

// Even faster PDF for cases where we know x is in a reasonable range
#[inline]
fn fast_standard_normal_pdf_unchecked(x: f64) -> f64 {
    // For x in [-8, 8], this is numerically stable
    // Beyond that range, PDF is essentially 0
    if x.abs() > 8.0 {
        0.0
    } else {
        fast_standard_normal_pdf(x)
    }
}

/// Analytical solution for the single party case
#[inline]
fn analytical_single_party_distance(p: f64) -> f64 {
    let sqrt_2_over_pi = 0.7978845608028654; // sqrt(2/π) ≈ 0.7979

    // First term: sqrt(2/π) * exp(-p²/2)
    // This represents the base expected distance modified by the party's position
    let exp_term = sqrt_2_over_pi * (-0.5 * p * p).exp();

    // Second term: p * erf(p/sqrt(2))
    // This accounts for the asymmetric shift when the party is not at the mean (p ≠ 0)
    // The factor 1/sqrt(2) ≈ 0.7071 standardizes p for the error function
    let erf_term = p * erf_approx(p * 0.7071067811865476); // p/sqrt(2)

    exp_term + erf_term
}

/// High-quality error function approximation using Abramowitz and Stegun method (formula 7.1.26)

#[inline]
fn erf_approx(x: f64) -> f64 {
    // Coefficients from Abramowitz and Stegun
    let a1 = 0.254829592;
    let a2 = -0.284496736;
    let a3 = 1.421413741;
    let a4 = -1.453152027;
    let a5 = 1.061405429;
    let p = 0.3275911;

    // Handle the sign: erf(-x) = -erf(x)
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();

    // Calculate t = 1/(1 + px)
    let t = 1.0 / (1.0 + p * x);

    // Apply the polynomial approximation
    // y = 1 - (((((a₅t + a₄)t + a₃)t + a₂)t + a₁)t * exp(-x²)
    let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x * x).exp();

    // Restore the sign
    sign * y
}

// SIMD-optimized integration for one segment
// Processes 4 integration points simultaneously
#[inline]
fn vectorized_integration_segment(left: f64, right: f64, party_pos: f64, num_points: usize) -> f64 {
    let step = (right - left) / num_points as f64;
    let mut segment_sum = 0.0;

    // SIMD constants - broadcast single values to 4-element vectors
    let inv_sqrt_2pi_vec = f64x4::splat(INV_SQRT_2PI);
    let party_pos_vec = f64x4::splat(party_pos);
    let step_vec = f64x4::splat(step);
    let left_vec = f64x4::splat(left);
    let neg_half = f64x4::splat(-0.5);

    // Process 4 points at once (main SIMD loop)
    let chunks = num_points / 4;
    for chunk in 0..chunks {
        // Calculate 4 integration point indices simultaneously
        let base_indices = f64x4::new([
            (chunk * 4) as f64 + 0.5,
            (chunk * 4) as f64 + 1.5,
            (chunk * 4) as f64 + 2.5,
            (chunk * 4) as f64 + 3.5,
        ]);

        // Calculate 4 x positions: left + (index * step)
        let x_vals = left_vec + base_indices * step_vec;

        // SIMD PDF calculation: calculate 4 exponential terms simultaneously
        let x_squared = x_vals * x_vals;
        let exp_vals = (neg_half * x_squared).exp();
        let pdf_vals = inv_sqrt_2pi_vec * exp_vals;

        // SIMD distance calculation: 4 absolute values simultaneously
        let distance_vals = (x_vals - party_pos_vec).abs();

        // Multiply distances * PDFs * step for all 4 points
        let contributions = distance_vals * pdf_vals * step_vec;

        // Sum the 4 contributions and add to total
        segment_sum += contributions.reduce_add();
    }

    // Handle leftover points (when num_points isn't divisible by 4)
    for j in (chunks * 4)..num_points {
        let x = left + (j as f64 + 0.5) * step;
        let pdf = fast_standard_normal_pdf_unchecked(x);
        let distance = (x - party_pos).abs();
        segment_sum += distance * pdf * step;
    }

    segment_sum
}

/// Given an array of voter positions and party positions,
/// computes the mean of the minimum absolute distances from each voter to the nearest party.
/// @export
#[extendr]
fn mean_min_abs_distance(voters: Vec<f64>, parties: Vec<f64>) -> f64 {
    let mut total_min = 0.0;

    for &v in &voters {
        let min_dist = parties
            .iter()
            .map(|&p| (v - p).abs())
            .fold(f64::INFINITY, |a, b| a.min(b));
        total_min += min_dist;
    }

    total_min / (voters.len() as f64)
}

/// Monte Carlo simulation entirely in Rust for fair benchmarking
/// Runs multiple simulations and returns the mean result
/// @export
#[extendr]
fn mean_min_abs_distance_mc_batch(
    parties: Vec<f64>,
    n_voters: i32,
    n_simulations: i32,
    seed: i32,
) -> f64 {
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;

    let mut rng = ChaCha8Rng::seed_from_u64(seed as u64);
    let mut total_result = 0.0;

    for _ in 0..n_simulations {
        let mut sim_total = 0.0;

        for _ in 0..n_voters {
            // Generate one voter position from standard normal
            let voter = rng.sample::<f64, _>(rand_distr::StandardNormal);

            // Find minimum distance to any party
            let min_dist = parties
                .iter()
                .map(|&p| (voter - p).abs())
                .fold(f64::INFINITY, |a, b| a.min(b));

            sim_total += min_dist;
        }

        total_result += sim_total / (n_voters as f64);
    }

    total_result / (n_simulations as f64)
}

/// Midpoint-based approach to computing expected distance to nearest party
///
/// This divides the voter space into segments between party midpoints,
/// then integrates the expected distance to the nearest party in each segment.
/// @export
#[extendr]
fn mean_min_abs_distance_midpoints(parties: Vec<f64>) -> f64 {
    // Return NaN if no parties provided
    if parties.is_empty() {
        return f64::NAN;
    }

    // Special case for single party - use analytical solution
    if parties.len() == 1 {
        return analytical_single_party_distance(parties[0]);
    }

    let range = 6.0; // beyond ±6σ, contribution is negligible

    // Sort the party positions
    let mut sorted_parties = parties;
    sorted_parties.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Pre-allocate boundaries vector
    let mut boundaries = Vec::with_capacity(sorted_parties.len() + 1);

    // Add a lower bound far to the left
    boundaries.push(-range);
    // Calculate midpoints between each pair of parties
    for i in 0..(sorted_parties.len() - 1) {
        // Midpoint between party i and i+1
        let mid = (sorted_parties[i] + sorted_parties[i + 1]) / 2.0;
        boundaries.push(mid);
    }
    // Add upper bound far to the right
    boundaries.push(range);

    // Sequential processing with SIMD optimization
    let mut total = 0.0;
    for i in 0..sorted_parties.len() {
        let left = boundaries[i];
        let right = boundaries[i + 1];
        let party_pos = sorted_parties[i];
        let segment_width = right - left;

        // Skip segments with negligible contribution
        if segment_width < 1e-6 {
            continue;
        }

        // Adaptive number of points based on where segement is in distribution
        let center_distance = ((left + right) / 2.0).abs();
        let base_points = match center_distance {
            x if x < 1.5 => 1500, // High precision for center region
            x if x < 3.0 => 800,  // Medium precision for mid region
            _ => 400,             // Lower precision for tail regions
        };

        // Scale by segment width, but cap the maximum
        let num_points = ((base_points as f64 * (segment_width / 2.0).min(1.0)).ceil() as usize)
            .max(50) // Minimum value
            .min(2000); // Maximum value

        // Use SIMD-optimized integration
        let segment_sum = vectorized_integration_segment(left, right, party_pos, num_points);
        total += segment_sum;
    }

    total
}

extendr_module! {
    mod minPartyDistance;
    fn mean_min_abs_distance_midpoints;
    fn mean_min_abs_distance;
    fn mean_min_abs_distance_mc_batch;
}
