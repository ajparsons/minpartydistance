#![allow(non_snake_case)]

use extendr_api::prelude::*;
use statrs::distribution::{Continuous, Normal};

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

    // Special case for single party
    if parties.len() == 1 {
        // For a single party at position p, the expected distance is just
        // the expected absolute deviation from p under a standard normal
        // E[|X - p|] where X ~ N(0,1)
        let p = parties[0];
        let norm = Normal::new(0.0, 1.0).unwrap();

        // Numerical integration over the entire real line
        let mut total = 0.0;
        let range = 8.0; // Integrate from -8 to 8 (covers ~99.999% of normal distribution)
        let num_points = 10000;
        let step = (2.0 * range) / num_points as f64;

        for i in 0..num_points {
            let x = -range + i as f64 * step;
            let pdf = norm.pdf(x);
            let distance = (x - p).abs();
            total += distance * pdf * step;
        }

        return total;
    }

    // Create the standard normal distribution
    let norm = Normal::new(0.0, 1.0).unwrap();

    // Sort the party positions
    let mut sorted_parties = parties.clone();
    sorted_parties.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut total = 0.0;

    // Special case for 2 parties - more careful boundary handling
    if sorted_parties.len() == 2 {
        let p1 = sorted_parties[0];
        let p2 = sorted_parties[1];
        let midpoint = (p1 + p2) / 2.0;

        // Region 1: (-∞, midpoint] - closest to p1
        let range = 8.0;
        let num_points = 5000;
        let left_bound = -range;
        let step1 = (midpoint - left_bound) / num_points as f64;

        for i in 0..num_points {
            let x = left_bound + i as f64 * step1;
            let pdf = norm.pdf(x);
            let distance = (x - p1).abs();
            total += distance * pdf * step1;
        }

        // Region 2: (midpoint, ∞) - closest to p2
        let right_bound = range;
        let step2 = (right_bound - midpoint) / num_points as f64;

        for i in 0..num_points {
            let x = midpoint + i as f64 * step2;
            let pdf = norm.pdf(x);
            let distance = (x - p2).abs();
            total += distance * pdf * step2;
        }

        return total;
    }

    // General case for 3+ parties
    let mut boundaries = Vec::new();

    // Add a lower bound far to the left
    boundaries.push(-8.0); // Instead of -∞, use -8 (covers 99.999% of normal)
    for i in 0..(sorted_parties.len() - 1) {
        // Midpoint between party i and i+1
        let mid = (sorted_parties[i] + sorted_parties[i + 1]) / 2.0;
        boundaries.push(mid);
    }
    // Add upper bound far to the right
    boundaries.push(8.0); // Instead of ∞, use 8

    // Now integrate each segment where one party is closest
    for i in 0..sorted_parties.len() {
        let left = boundaries[i];
        let right = boundaries[i + 1];
        let party_pos = sorted_parties[i];

        // Higher precision numerical integration
        let num_points = 5000; // Increased from 1000
        let step = (right - left) / num_points as f64;
        let mut segment_sum = 0.0;

        for j in 0..num_points {
            let x = left + j as f64 * step;
            let pdf = norm.pdf(x);
            let distance = (x - party_pos).abs();
            segment_sum += distance * pdf * step;
        }

        total += segment_sum;
    }

    total
}

extendr_module! {
    mod MinPartyDistance;
    fn mean_min_abs_distance_midpoints;
    fn mean_min_abs_distance;
}
