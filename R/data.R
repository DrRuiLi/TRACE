#' TRACE low-C cutoff table
#'
#' Lookup table used by TRACE CN-pattern filtering to constrain possible
#' carbon counts by precursor mass range.
#'
#' @format A data frame with columns `c.count`, `mass_min`, and `mass_max`.
#' @source Internal TRACE preprocessing workflow.
"TRACE_LowC_cutoff"
