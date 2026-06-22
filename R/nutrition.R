#' Match m/z / RT and extract all CN labeling forms for a metabolite
#'
#' @param pave_formula Metabolite CN formula (e.g. `"C10N5"`).
#' @param seed_mz Seed feature m/z.
#' @param seed_rt Seed feature RT (seconds).
#' @param xcms_fdf xcms feature definition data.frame with `mzmed`, `rtmed`.
#' @param xcms_mat Intensity matrix (features × samples).
#' @param xcms_se SummarizedExperiment with `sample.name` and `group`.
#' @param ppm m/z tolerance in ppm.
#' @param rt_tol RT tolerance in seconds.
#'
#' @returns A list with `cn_table`, `exp_mat`, `c_max`, and `n_max`, or
#'   `NULL` when fewer than two labeling forms are matched.
#' @export
extract_cn_labeling_forms <- function(
    pave_formula,
    seed_mz,
    seed_rt,
    xcms_fdf,
    xcms_mat,
    xcms_se,
    ppm = 3,
    rt_tol = 5) {
  cn_count <- MSCC::chemform_parse(pave_formula)
  c_max <- cn_count[1, "C"]
  n_max <- cn_count[1, "N"]

  cn_diff <- get_CN_mass_diff_table(c_max, n_max)
  cn_diff[, mz := mass_diff + seed_mz]

  cn_match <- MSdev::match_mz_foverlaps(
    cn_diff$mz,
    xcms_fdf$mzmed,
    ppm = ppm
  )
  cn_match[, fid := ion2][, rt := xcms_fdf$rtmed[fid]]

  cn_table <- cbind(cn_match, cn_diff[cn_match$ion1])
  cn_table[, label_pattern := paste0("C", C_count, "N", N_count)][!is.na(fid)]

  cn_table <- cn_table[abs(rt - seed_rt) < rt_tol]
  cn_table <- cn_table %>%
    dplyr::group_by(label_pattern) %>%
    dplyr::slice_min(mz.ppm, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(N_count, C_count)

  if (nrow(cn_table) <= 1) {
    return(NULL)
  }

  exp_mat <- xcms_mat[cn_table$fid, xcms_se$sample.name, drop = FALSE]
  rownames(exp_mat) <- cn_table$label_pattern

  list(
    cn_table = as.data.frame(cn_table),
    exp_mat = exp_mat,
    c_max = c_max,
    n_max = n_max
  )
}

#' Calculate absorb-form fractions per sample
#'
#' For each labeling form `CxNy`, the absorb form is `C(Cmax - x)N(Nmax - y)`.
#' Intensities are normalized to a reference group and converted to fractions
#' within each sample.
#'
#' @param cn_result Output of [extract_cn_labeling_forms()].
#' @param xcms_se SummarizedExperiment with `sample.name`, `sample.type`, `group`.
#' @param norm_group Sample group used for intensity normalization.
#' @param norm_replicates Divisor for the normalization sum (replicate count).
#'
#' @returns A data.frame with absorb fractions and patterns per sample.
#' @export
compute_absorb_forms <- function(
    cn_result,
    xcms_se,
    norm_group = "NSA",
    norm_replicates = 10) {
  cn_table <- cn_result$cn_table
  exp_mat <- cn_result$exp_mat
  c_max <- cn_result$c_max
  n_max <- cn_result$n_max

  norm_samples <- xcms_se$sample.name[xcms_se$group == norm_group]
  norm_to <- sum(exp_mat[, norm_samples, drop = FALSE], na.rm = TRUE) / norm_replicates
  if (!is.finite(norm_to) || norm_to <= 0) {
    norm_to <- mean(exp_mat, na.rm = TRUE)
  }
  exp_mat <- exp_mat / norm_to

  cn_table$c_absorb <- c_max - cn_table$C_count
  cn_table$n_absorb <- n_max - cn_table$N_count
  cn_table$absorb_pattern <- paste0("C", cn_table$c_absorb, "N", cn_table$n_absorb)

  exp_mat %>%
    as.data.frame() %>%
    tibble::rownames_to_column("label_pattern") %>%
    tidyr::pivot_longer(-label_pattern, names_to = "sample", values_to = "intensity") %>%
    dplyr::left_join(cn_table, by = "label_pattern") %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(
      group = unname(stats::setNames(xcms_se$sample.type, xcms_se$sample.name)[sample]),
      int_total = sum(intensity, na.rm = TRUE),
      absorb_fraction = ifelse(int_total > 0, intensity / int_total, NA_real_)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      c_max = c_max,
      n_max = n_max
    )
}

#' Drop fully labeled absorb form C0N0
#'
#' @param df Absorb-form data.frame from [compute_absorb_forms()].
#' @param renormalize If `TRUE`, rescale remaining fractions to sum to 1
#'   within each metabolite sample (composition plots). If `FALSE`, keep
#'   original absorb ratios (global sum plots).
#'
#' @returns Filtered data.frame.
#' @export
drop_absorb_c0n0 <- function(df, renormalize = TRUE) {
  df <- df %>%
    dplyr::filter(absorb_pattern != "C0N0")

  if (!renormalize || !nrow(df)) {
    return(df)
  }

  renorm_keys <- intersect(c("feature_id", "polarity", "sample"), names(df))
  if (!length(renorm_keys)) {
    renorm_keys <- intersect(c("feature_id", "polarity", "group"), names(df))
  }

  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(renorm_keys))) %>%
    dplyr::mutate(
      absorb_fraction = absorb_fraction / sum(absorb_fraction, na.rm = TRUE)
    ) %>%
    dplyr::ungroup()
}

#' Run absorb-form analysis for one polarity
#'
#' @param obj MSdev object with TRACE results and xcms data.
#' @param i.pol Polarity index (`0` = negative, `1` = positive).
#' @param ppm m/z tolerance for CN form matching.
#' @param rt_tol RT tolerance in seconds.
#' @param min_peak_maxo Minimum peak intensity filter.
#' @param nutrient_groups Nutrient sample groups to retain.
#'
#' @returns A data.frame of absorb-form fractions for all metabolite seeds.
#' @export
run_absorb_form_analysis <- function(
    obj,
    i.pol = 1,
    ppm = 3,
    rt_tol = 5,
    min_peak_maxo = 1e4,
    nutrient_groups = c("NSA", "NSB", "NSC", "NSD", "NSE", "NSF", "NSG", "NSH")) {
  pol <- if (i.pol == 0) "Negative" else "Positive"
  xcms_obj <- obj@xcmsData[[paste0(pol, "MS1")]]
  pave_seed <- obj@statData$TRACE[[pol]]
  if (!is.data.frame(pave_seed)) {
    pave_seed <- as.data.frame(pave_seed)
  }
  pave_seed <- pave_seed %>%
    dplyr::filter(type == "metabolite")

  xcms_se_raw <- MSdev::get_xcms_feature_se(xcms_obj, missing = 1)
  xcms_se <- xcms_se_raw[, grepl("NS", xcms_se_raw$sample.type)]
  xcms_se <- MSdev::se_adjuset_by_weight(xcms_se)

  xcms_fdf <- as.data.frame(SummarizedExperiment::rowData(xcms_se))
  xcms_fdf$peakMaxo <- rowMeans(
    SummarizedExperiment::assay(
      xcms_se_raw[, xcms_se_raw$group == "S12C14N", drop = FALSE]
    )
  )
  xcms_mat <- SummarizedExperiment::assay(xcms_se)

  results <- vector("list", nrow(pave_seed))
  pb <- MSdev:::get_progress_bar(nrow(pave_seed))

  for (i in seq_len(nrow(pave_seed))) {
    pb$tick()
    i_fid <- as.numeric(pave_seed$feature_id[i])
    if (!is.finite(i_fid) || i_fid < 1 || i_fid > nrow(xcms_fdf)) next

    cn_result <- extract_cn_labeling_forms(
      pave_formula = pave_seed$pave_formula[i],
      seed_mz = xcms_fdf$mzmed[i_fid],
      seed_rt = pave_seed$rt[i],
      xcms_fdf = xcms_fdf,
      xcms_mat = xcms_mat,
      xcms_se = xcms_se,
      ppm = ppm,
      rt_tol = rt_tol
    )
    if (is.null(cn_result)) next

    absorb_df <- compute_absorb_forms(cn_result, xcms_se)
    absorb_df$feature_id <- pave_seed$feature_id[i]
    absorb_df$compound_id <- if ("compound_id" %in% names(pave_seed)) {
      pave_seed$compound_id[i]
    } else {
      NA_character_
    }
    absorb_df$fid <- paste0(tolower(substr(pol, 1, 3)), i_fid)
    absorb_df$metabolite <- pave_seed$name[i]
    absorb_df$kegg_id <- if ("kegg_id" %in% names(pave_seed)) {
      pave_seed$kegg_id[i]
    } else {
      NA_character_
    }
    absorb_df$pave_formula <- pave_seed$pave_formula[i]
    absorb_df$peak_maxo <- xcms_fdf$peakMaxo[i_fid]
    absorb_df$polarity <- pol
    results[[i]] <- absorb_df
  }

  dplyr::bind_rows(results) %>%
    dplyr::filter(
      group %in% nutrient_groups,
      peak_maxo > min_peak_maxo,
      is.finite(absorb_fraction)
    )
}

#' Build a display label for a nutrient-analysis metabolite
#'
#' @param metabolite Metabolite name.
#' @param pave_formula CN formula.
#' @param feature_id TRACE feature id.
#' @param polarity `"Positive"` or `"Negative"`.
#'
#' @returns Character label.
#' @export
nutrition_met_label <- function(metabolite, pave_formula, feature_id, polarity) {
  nm <- dplyr::coalesce(as.character(metabolite), as.character(pave_formula))
  paste0(nm, " (", polarity, " #", feature_id, ")")
}

#' Summarize absorb forms with atom-count-weighted scores
#'
#' Per metabolite, weighted contribution is `absorb_ratio * (c_absorb +
#' n_absorb)`. For absorb form `C5N5` with ratio `0.2`, the weight is
#' `10 * 0.2 = 2`. Values are summed across metabolites per group and
#' absorb form.
#'
#' @param absorb_df Output from [run_absorb_form_analysis()] after
#'   [drop_absorb_c0n0()].
#' @param nutrient_groups Optional character vector of groups to retain.
#'
#' @returns A data.frame with `sum_absorb_ratio` and `sum_weighted_absorb`.
#' @export
nutrition_summarize_absorb_global <- function(absorb_df, nutrient_groups = NULL) {
  if (!is.null(nutrient_groups)) {
    absorb_df <- absorb_df %>%
      dplyr::filter(group %in% nutrient_groups)
  }

  absorb_df %>%
    dplyr::mutate(atom_count = c_absorb + n_absorb) %>%
    dplyr::group_by(feature_id, polarity, group, absorb_pattern, c_absorb, n_absorb, atom_count) %>%
    dplyr::summarize(
      mean_absorb_ratio = mean(absorb_fraction, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(weighted_absorb = mean_absorb_ratio * atom_count) %>%
    dplyr::group_by(group, absorb_pattern, c_absorb, n_absorb, atom_count) %>%
    dplyr::summarize(
      sum_absorb_ratio = sum(mean_absorb_ratio, na.rm = TRUE),
      sum_weighted_absorb = sum(weighted_absorb, na.rm = TRUE),
      .groups = "drop"
    )
}

#' Select top absorb forms per nutrient group
#'
#' @param absorb_global Output of [nutrition_summarize_absorb_global()].
#' @param groups Group ids to select from (e.g. supplement conditions).
#' @param top_n Number of forms per group.
#' @param rank_by Rank by `"weighted"` (`sum_weighted_absorb`) or `"ratio"`
#'   (`sum_absorb_ratio`).
#'
#' @returns A data.frame of selected group–form rows.
#' @export
nutrition_select_top_absorb_forms <- function(
    absorb_global,
    groups,
    top_n = 10,
    rank_by = c("weighted", "ratio")) {
  rank_by <- match.arg(rank_by)
  score_col <- if (rank_by == "weighted") "sum_weighted_absorb" else "sum_absorb_ratio"

  absorb_global %>%
    dplyr::filter(group %in% groups) %>%
    dplyr::group_by(group) %>%
    dplyr::arrange(dplyr::desc(.data[[score_col]]), absorb_pattern) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()
}

#' Merge top absorb forms selected per group
#'
#' @param top_per_group Output of [nutrition_select_top_absorb_forms()].
#' @param absorb_global Output of [nutrition_summarize_absorb_global()].
#' @param rank_by Score used to order merged forms.
#'
#' @returns A list with `forms` (ordered absorb patterns) and `detail`
#'   (full rows for plotting).
#' @export
nutrition_merge_top_absorb_forms <- function(
    top_per_group,
    absorb_global,
    rank_by = c("weighted", "ratio")) {
  rank_by <- match.arg(rank_by)
  score_col <- if (rank_by == "weighted") "sum_weighted_absorb" else "sum_absorb_ratio"

  merged.forms <- absorb_global %>%
    dplyr::filter(absorb_pattern %in% top_per_group$absorb_pattern) %>%
    dplyr::group_by(absorb_pattern, c_absorb, n_absorb, atom_count) %>%
    dplyr::summarize(
      total = sum(.data[[score_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(total))

  list(
    forms = merged.forms,
    patterns = merged.forms$absorb_pattern
  )
}
