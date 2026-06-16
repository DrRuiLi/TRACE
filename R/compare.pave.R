#' Test whether a PAVE/TRACE annotation type is non-blank
#' @noRd
.compare_pave_nonblank <- function(x) {
  !is.na(x) & x != "Blank/noise"
}

#' Classify PAVE/TRACE features into CN comparison seed groups
#'
#' @param pave.trace.merged Data frame from a PAVE/TRACE feature merge.
#'
#' @returns A data frame with one row per seed `feature_id` and comparison
#'   group `cn_cmp`.
#' @noRd
.compare_pave_cn_seed_groups <- function(pave.trace.merged) {
  pave.trace.merged %>%
    dplyr::mutate(
      cn_cmp = dplyr::case_when(
        .compare_pave_nonblank(type_pave) & !.compare_pave_nonblank(type_trace) ~
          "unmatched_blank",
        .compare_pave_nonblank(type_pave) & .compare_pave_nonblank(type_trace) ~
          "matched_cn",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(cn_cmp)) %>%
    dplyr::distinct(feature_id, cn_cmp)
}

#' Aggregate absolute m/z error across CN bundle edges
#' @noRd
.compare_pave_mz_error <- function(x, mz_error = c("max", "mean")) {
  mz_error <- match.arg(mz_error)
  x <- abs(x)
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(NA_real_)
  }
  if (mz_error == "max") {
    max(x)
  } else {
    mean(x)
  }
}

#' Summarize CN bundles for selected seeds from TRACE CNfinder results
#'
#' Extracts CN pattern bundles whose seed (`from`) is in `seed_ids`. Bundle
#' edges are the CN-label rows with a finite `TRACE_cor`. The m/z error across
#' bundle edges is summarized by `mz_error`.
#'
#' @param object MSdev object with TRACE results in `object@advancedAna$TRACE_temp`.
#' @param seed_ids Integer vector of CN seed feature IDs.
#' @param i.pol Polarity index. `0` for negative and `1` for positive.
#' @param mz_error Character. Use `"max"` or `"mean"` absolute m/z error (ppm)
#'   across bundle edges.
#'
#' @returns A data frame with one row per CN seed (`bundle_seed`) and columns
#'   `cn_bundle_cor` and `cn_bundle_mz_ppm`.
#' @noRd
.compare_pave_cn_bundles <- function(
    object,
    seed_ids,
    i.pol = 0,
    mz_error = c("max", "mean")) {
  mz_error <- match.arg(mz_error)
  pol <- .trace_get_pol(i.pol)
  cn.finder <- object@advancedAna$TRACE_temp[[pol]]$CNfinder
  if (is.null(cn.finder) || !nrow(cn.finder)) {
    stop(
      "CNfinder not found in object@advancedAna$TRACE_temp$", pol,
      ". Run TRACE_workflow() first."
    )
  }

  seed_ids <- unique(as.integer(seed_ids))
  seed_ids <- seed_ids[!is.na(seed_ids)]
  if (!length(seed_ids)) {
    return(data.frame(
      bundle_seed = integer(),
      cn_bundle_cor = numeric(),
      cn_bundle_mz_ppm = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  cn.finder %>%
    as.data.frame() %>%
    dplyr::filter(
      as.integer(from) %in% seed_ids,
      is.finite(TRACE_cor)
    ) %>%
    dplyr::group_by(from) %>%
    dplyr::summarise(
      cn_bundle_cor = max(TRACE_cor),
      cn_bundle_mz_ppm = .compare_pave_mz_error(mz.ppm, mz_error = mz_error),
      .groups = "drop"
    ) %>%
    dplyr::transmute(
      bundle_seed = as.integer(as.character(from)),
      cn_bundle_cor = cn_bundle_cor,
      cn_bundle_mz_ppm = cn_bundle_mz_ppm
    )
}

#' Build PAVE vs TRACE CN bundle comparison data
#'
#' Classify CN seeds as `unmatched_blank` (PAVE non-blank, TRACE blank) or
#' `matched_cn` (non-blank in both), then extract each seed's CN bundle from the
#' MSdev object. Bundle correlation and m/z error are calculated from the CN
#' pattern edges in each bundle.
#'
#' @param pave.trace.merged Data frame from a PAVE/TRACE feature merge.
#' @param object MSdev object with TRACE results in `object@advancedAna$TRACE_temp`.
#' @param i.pol Polarity index. `0` for negative and `1` for positive.
#' @param mz_error Character. Use `"max"` or `"mean"` absolute m/z error (ppm)
#'   across bundle edges.
#'
#' @returns A data frame with one row per CN seed (`feature_id`), comparison
#'   group (`cn_cmp`), and bundle metrics `cn_bundle_cor` / `cn_bundle_mz_ppm`.
#' @export
compare_pave_cn_bundle_data <- function(
    pave.trace.merged,
    object,
    i.pol = 0,
    mz_error = c("max", "mean")) {
  mz_error <- match.arg(mz_error)
  seed.df <- .compare_pave_cn_seed_groups(pave.trace.merged)
  cn.bundle <- .compare_pave_cn_bundles(
    object = object,
    seed_ids = seed.df$feature_id,
    i.pol = i.pol,
    mz_error = mz_error
  )

  seed.df %>%
    dplyr::inner_join(cn.bundle, by = c("feature_id" = "bundle_seed")) %>%
    dplyr::mutate(
      cn_cmp = factor(
        cn_cmp,
        levels = c("unmatched_blank", "matched_cn"),
        labels = c("PAVE non-noise / TRACE noise", "Both non-noise")
      )
    )
}

#' Plot CN bundle correlation and m/z error for PAVE vs TRACE groups
#'
#' @param pave.trace.merged Data frame from a PAVE/TRACE feature merge.
#' @param object MSdev object with TRACE results in `object@advancedAna$TRACE_temp`.
#' @param i.pol Polarity index. `0` for negative and `1` for positive.
#' @param mz_error Character. Use `"max"` or `"mean"` absolute m/z error (ppm)
#'   across bundle edges.
#' @param return.data Logical. If `TRUE`, return a list with `plot`, `p_cor`,
#'   `p_mz`, and `data`.
#'
#' @returns A patchwork object combining correlation and m/z error boxplots, or
#'   a list when `return.data = TRUE`.
#' @export
plot_compare_pave_cn_bundle <- function(
    pave.trace.merged,
    object,
    i.pol = 0,
    mz_error = c("max", "mean"),
    return.data = FALSE) {
  mz_error <- match.arg(mz_error)
  bundle.plot.df <- compare_pave_cn_bundle_data(
    pave.trace.merged = pave.trace.merged,
    object = object,
    i.pol = i.pol,
    mz_error = mz_error
  )

  p.cor <- bundle.plot.df %>%
    dplyr::filter(is.finite(cn_bundle_cor)) %>%
    ggplot2::ggplot(ggplot2::aes(x = cn_cmp, y = cn_bundle_cor)) +
    ggplot2::geom_boxplot(outlier.size = 0.4, linewidth = 0.3) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.25, size = 0.4) +
    ggplot2::labs(x = NULL, y = "Cor") +
    ggplot2::theme_bw(base_size = 6)

  p.mz <- ggplot2::ggplot(
    bundle.plot.df,
    ggplot2::aes(x = cn_cmp, y = cn_bundle_mz_ppm)
  ) +
    ggplot2::geom_boxplot(linewidth = 0.3) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.25, size = 0.4) +
    ggplot2::labs(x = NULL, y = "mz error (ppm)") +
    ggplot2::theme_bw(base_size = 6)

  p.all <- p.cor + p.mz

  if (return.data) {
    return(list(plot = p.all, p_cor = p.cor, p_mz = p.mz, data = bundle.plot.df))
  }

  p.all
}
