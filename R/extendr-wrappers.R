# Generated by extendr: Do not edit by hand

# nolint start

#
# This file was created with the following call:
#   .Call("wrap__make_MinPartyDistance_wrappers", use_symbols = TRUE, package_name = "MinPartyDistance")

#' @usage NULL
#' @useDynLib MinPartyDistance, .registration = TRUE
NULL

#' Midpoint-based approach to computing expected distance to nearest party
#'
#' This divides the voter space into segments between party midpoints,
#' then integrates the expected distance to the nearest party in each segment.
#' @export
mean_min_abs_distance_midpoints <- function(parties) .Call(wrap__mean_min_abs_distance_midpoints, parties)

#' Given an array of voter positions and party positions,
#' computes the mean of the minimum absolute distances from each voter to the nearest party.
#' @export
mean_min_abs_distance <- function(voters, parties) .Call(wrap__mean_min_abs_distance, voters, parties)


# nolint end
