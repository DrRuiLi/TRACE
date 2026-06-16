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

#' Map PAVE annotation to metabolite (M) or adduct (A)
#' @noRd
.compare_pave_ma_pave <- function(type_pave) {
  dplyr::case_when(
    type_pave == "Metabolite" ~ "M",
    type_pave %in% c("Adduct", "Fragment", "Isotope", "Dimer", "Multicharge", "Low_c", "Low_score") ~ "A",
    TRUE ~ NA_character_
  )
}

#' Map TRACE annotation to metabolite (M) or adduct (A)
#' @noRd
.compare_pave_ma_trace <- function(type_trace) {
  dplyr::case_when(
    tolower(type_trace) == "metabolite" ~ "M",
    tolower(type_trace) == "adduct" ~ "A",
    TRUE ~ NA_character_
  )
}

#' Default protonated adduct label for PAVE barplots
#' @noRd
.compare_pave_default_adduct <- function(i.pol = 0) {
  if (i.pol == 0) "[M-H]-" else "[M+H]+"
}

#' Feature intensity from xcms labeled samples
#' @noRd
.compare_pave_feature_intensity <- function(object, i.pol = 0) {
  pol <- .trace_get_pol(i.pol)
  xcms.xcms <- object@xcmsData[[paste0(pol, "MS1")]]
  if (is.null(xcms.xcms)) {
    stop("xcms data not found for polarity ", pol, ".")
  }

  xcms.val <- xcms::featureValues(xcms.xcms, missing = 0, value = "maxo")
  sample.types <- c("S12C14N", "S12C15N", "S13C14N", "S13C15N")
  sample.idx <- Biobase::pData(xcms.xcms)$sample.type %in% sample.types
  if (!any(sample.idx)) {
    stop("No labeled sample groups found in xcms metadata.")
  }

  stats::setNames(
    rowMeans(xcms.val[, sample.idx, drop = FALSE], na.rm = TRUE),
    as.character(seq_len(nrow(xcms.val)))
  )
}

#' Build PAVE vs TRACE adduct selection comparison data
#'
#' For features in the same TRACE CN group (`seed_trace`), classify each
#' feature as `M/M`, `M/A`, `A/M`, or `A/A` based on TRACE and PAVE metabolite
#' versus adduct annotation (`TRACE/PAVE`).
#'
#' @param pave.trace.merged Data frame from a PAVE/TRACE feature merge.
#' @param object MSdev object with xcms results.
#' @param i.pol Polarity index. `0` for negative and `1` for positive.
#'
#' @returns A data frame with CN-group features, `ma_group`, adduct labels, and
#'   `intensity`.
#' @export
compare_pave_adduct_data <- function(pave.trace.merged, object, i.pol = 0) {
  intensity <- .compare_pave_feature_intensity(object, i.pol = i.pol)
  pave_adduct <- .compare_pave_default_adduct(i.pol = i.pol)

  pave.trace.merged %>%
    dplyr::filter(
      !is.na(seed_trace),
      seed_trace != "",
      !is.na(type_trace),
      type_trace != "",
      type_trace != "Blank/noise"
    ) %>%
    dplyr::mutate(
      pave_ma = .compare_pave_ma_pave(type_pave),
      trace_ma = .compare_pave_ma_trace(type_trace),
      ma_group = dplyr::if_else(
        is.na(pave_ma) | is.na(trace_ma),
        NA_character_,
        paste0(trace_ma, "/", pave_ma)
      ),
      intensity = as.numeric(intensity[as.character(feature_id)]),
      pave_adduct = pave_adduct,
      trace_adduct = as.character(adduct)
    ) %>%
    dplyr::filter(!is.na(ma_group)) %>%
    dplyr::mutate(
      ma_group = factor(
        ma_group,
        levels = c("M/M", "M/A", "A/M", "A/A")
      )
    )
}

#' Adduct type counts for PAVE and TRACE barplots
#'
#' @param adduct.df Output of [compare_pave_adduct_data()].
#'
#' @returns A data frame with columns `source`, `adduct`, and `n`.
#' @export
compare_pave_adduct_distribution <- function(adduct.df) {
  pave.dist <- adduct.df %>%
    dplyr::count(pave_adduct, name = "n") %>%
    dplyr::transmute(
      source = "PAVE",
      adduct = pave_adduct,
      n = n
    )

  trace.dist <- adduct.df %>%
    dplyr::mutate(
      trace_adduct = dplyr::coalesce(
        dplyr::na_if(trace_adduct, ""),
        pave_adduct
      )
    ) %>%
    dplyr::count(trace_adduct, name = "n") %>%
    dplyr::transmute(
      source = "TRACE",
      adduct = trace_adduct,
      n = n
    )

  dplyr::bind_rows(pave.dist, trace.dist) %>%
    dplyr::mutate(
      source = factor(source, levels = c("PAVE", "TRACE")),
      adduct = factor(adduct)
    )
}

#' Plot PAVE vs TRACE adduct selection comparison
#'
#' @param pave.trace.merged Data frame from a PAVE/TRACE feature merge.
#' @param object MSdev object with xcms results.
#' @param i.pol Polarity index. `0` for negative and `1` for positive.
#' @param top.trace.adducts Maximum number of TRACE adduct labels to show.
#' @param return.data Logical. If `TRUE`, return a list with `plot`, `p_adduct`,
#'   `p_intensity`, `data`, and `adduct_distribution`.
#'
#' @returns A patchwork object with a stacked adduct barplot (PAVE and TRACE
#'   columns) and an intensity boxplot, or a list when `return.data = TRUE`.
#' @export
plot_compare_pave_adduct <- function(
    pave.trace.merged,
    object,
    i.pol = 0,
    top.trace.adducts = 10L,
    return.data = FALSE) {
  adduct.df <- compare_pave_adduct_data(
    pave.trace.merged = pave.trace.merged,
    object = object,
    i.pol = i.pol
  )
  adduct.dist <- compare_pave_adduct_distribution(adduct.df)

  trace.top <- adduct.dist %>%
    dplyr::filter(source == "TRACE") %>%
    dplyr::slice_max(n, n = top.trace.adducts, with_ties = FALSE) %>%
    dplyr::pull(adduct) %>%
    as.character()

  adduct.dist.plot <- adduct.dist %>%
    dplyr::mutate(
      adduct_plot = dplyr::case_when(
        source == "PAVE" ~ as.character(adduct),
        adduct %in% trace.top ~ as.character(adduct),
        TRUE ~ "Other"
      )
    ) %>%
    dplyr::group_by(source, adduct_plot) %>%
    dplyr::summarise(n = sum(n), .groups = "drop") %>%
    dplyr::mutate(
      source = factor(source, levels = c("PAVE", "TRACE")),
      adduct_plot = factor(adduct_plot)
    )

  p.adduct <- ggplot2::ggplot(
    adduct.dist.plot,
    ggplot2::aes(x = source, y = n, fill = adduct_plot)
  ) +
    ggplot2::geom_col(position = "stack", width = 0.6, colour = "white", linewidth = 0.2) +
    ggplot2::labs(x = NULL, y = "Count", fill = "Adduct") +
    ggplot2::theme_bw(base_size = 6) +
    ggplot2::theme(
      legend.position = "right",
      legend.key.size = grid::unit(0.15, "inch"),
      legend.text = ggplot2::element_text(size = 5)
    )

  p.intensity <- adduct.df %>%
    dplyr::filter(is.finite(intensity), intensity > 0) %>%
    ggplot2::ggplot(ggplot2::aes(x = ma_group, y = intensity, fill = ma_group, color = ma_group)) +
    ggplot2::geom_boxplot(
      outlier.size = 0.4,
      linewidth = 0.3,
      alpha = 0.4,
      show.legend = FALSE
    ) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.35, size = 0.3, show.legend = FALSE) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "TRACE/PAVE", y = "Intensity") +
    ggplot2::theme_bw(base_size = 6)

  p.all <- p.adduct + p.intensity +
    patchwork::plot_layout(widths  = c(1.2, 1))

  if (return.data) {
    return(list(
      plot = p.all,
      p_adduct = p.adduct,
      p_intensity = p.intensity,
      data = adduct.df,
      adduct_distribution = adduct.dist
    ))
  }

  p.all
}

#' Circular barplot from a count table using ggplot2
#' @noRd
.compare_pave_polar_barplot <- function(count.df, title = NULL) {
  if (!nrow(count.df)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, label = "No data", size = 3) +
        ggplot2::labs(title = title) +
        ggplot2::theme_void(base_size = 6)
    )
  }

  count.df <- count.df[order(-count.df$n), , drop = FALSE]
  count.df$label <- factor(count.df$label, levels = count.df$label)

  ggplot2::ggplot(count.df, ggplot2::aes(x = 1, y = n, fill = label)) +
    ggplot2::geom_bar(
      stat = "identity",
      position = "stack",
      width = 0.8,
      color = "white",
      linewidth = 0.2
    ) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::labs(title = title, x = NULL, y = NULL, fill = NULL) +
    ggplot2::theme_void(base_size = 6) +
    ggplot2::theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.key.size = grid::unit(0.12, "inch"),
      legend.text = ggplot2::element_text(size = 4)
    )
}

#' Build PAVE fragment vs TRACE annotation comparison data
#'
#' Selects features annotated as fragment in PAVE but not in TRACE, then
#' summarizes TRACE adduct labels and fragment-network candidacy in `fg.net`.
#'
#' @param pave.trace.merged Data frame from a PAVE/TRACE feature merge.
#' @param object MSdev object with TRACE results in `object@advancedAna$TRACE_temp`.
#' @param i.pol Polarity index. `0` for negative and `1` for positive.
#'
#' @returns A list with `all`, `trace_adduct`, `trace_metabolite`, and merged
#'   `trace_adduct_metabolite` data frames.
#' @export
compare_pave_fragment_data <- function(pave.trace.merged, object, i.pol = 0) {
  pol <- .trace_get_pol(i.pol)
  fg.net <- object@advancedAna$TRACE_temp[[pol]]$fg.net
  if (is.null(fg.net) || !nrow(fg.net)) {
    stop(
      "fg.net not found in object@advancedAna$TRACE_temp$", pol,
      ". Run TRACE_workflow() first."
    )
  }

  fg.ids <- unique(c(
    as.integer(fg.net$from),
    as.integer(fg.net$to)
  ))
  pave_adduct <- .compare_pave_default_adduct(i.pol = i.pol)

  frag.df <- pave.trace.merged %>%
    dplyr::filter(
      type_pave == "Fragment",
      is.na(type_trace) | type_trace != "Fragment"
    ) %>%
    dplyr::mutate(
      trace_adduct = dplyr::coalesce(
        dplyr::na_if(as.character(adduct), ""),
        pave_adduct
      ),
      in_fg_net = feature_id %in% fg.ids,
      fg_net_status = dplyr::if_else(
        in_fg_net,
        "In fg.net",
        "Not in fg.net"
      )
    )

  list(
    all = frag.df,
    trace_adduct = frag.df %>%
      dplyr::filter(tolower(type_trace) == "adduct"),
    trace_metabolite = frag.df %>%
      dplyr::filter(tolower(type_trace) == "metabolite"),
    trace_adduct_metabolite = frag.df %>%
      dplyr::filter(tolower(type_trace) %in% c("adduct", "metabolite")) %>%
      dplyr::mutate(
        trace_group = dplyr::if_else(
          tolower(type_trace) == "adduct",
          "TRACE adduct",
          "TRACE metabolite"
        )
      )
  )
}

#' Summarize TRACE adduct counts for polar barplots
#' @noRd
.compare_pave_trace_adduct_count <- function(df, top.trace.adducts = 8L) {
  count.df <- df %>%
    dplyr::count(trace_adduct, name = "n") %>%
    dplyr::arrange(dplyr::desc(n))

  if (!nrow(count.df)) {
    return(data.frame(label = character(), n = integer(), stringsAsFactors = FALSE))
  }

  if (nrow(count.df) > top.trace.adducts) {
    top.labels <- count.df$trace_adduct[seq_len(top.trace.adducts)]
    count.df <- count.df %>%
      dplyr::mutate(
        label = dplyr::if_else(trace_adduct %in% top.labels, trace_adduct, "Other")
      ) %>%
      dplyr::group_by(label) %>%
      dplyr::summarise(n = sum(n), .groups = "drop")
  } else {
    count.df$label <- count.df$trace_adduct
  }

  count.df
}

#' Plot PAVE fragment vs TRACE annotation comparison
#'
#' @param pave.trace.merged Data frame from a PAVE/TRACE feature merge.
#' @param object MSdev object with TRACE results in `object@advancedAna$TRACE_temp`.
#' @param i.pol Polarity index. `0` for negative and `1` for positive.
#' @param top.trace.adducts Maximum number of TRACE adduct labels to show in
#'   each polar barplot.
#' @param return.data Logical. If `TRUE`, return a list with plots and data.
#'
#' @returns A patchwork object with one polar barplot of TRACE adduct type from
#'   merged `trace_adduct_metabolite` features and one polar barplot of
#'   `fg.net` candidacy for `trace_metabolite`, or a list when `return.data`
#'   is `TRUE`.
#' @export
plot_compare_pave_fragment <- function(
    pave.trace.merged,
    object,
    i.pol = 0,
    top.trace.adducts = 8L,
    return.data = FALSE) {
  frag.data <- compare_pave_fragment_data(
    pave.trace.merged = pave.trace.merged,
    object = object,
    i.pol = i.pol
  )

  adduct.count <- .compare_pave_trace_adduct_count(
    frag.data$trace_adduct_metabolite,
    top.trace.adducts = top.trace.adducts
  )
  fg.net.count <- frag.data$trace_metabolite %>%
    dplyr::count(fg_net_status, name = "n") %>%
    dplyr::mutate(
      label = dplyr::case_when(
        fg_net_status == "In fg.net" ~ "Re-assigned candidate",
        fg_net_status == "Not in fg.net" ~ "Non candidate",
        TRUE ~ as.character(fg_net_status)
      ),
      label = factor(label, levels = c("Re-assigned candidate", "Non candidate"))
    )

  p.adduct <- .compare_pave_polar_barplot(
    adduct.count,
    title = "PAVE fragment >> TRACE adduct / metabolite"
  )

  p.metabolite.fg <- .compare_pave_polar_barplot(
    fg.net.count,
    title = "PAVE fragment >> TRACE metabolite"
  ) +
    ggplot2::labs(fill = "Fragment candidate")

  p.all <- p.adduct + p.metabolite.fg +
    patchwork::plot_layout(width = c(1, 1))

  if (return.data) {
    return(list(
      plot = p.all,
      p_adduct = p.adduct,
      p_metabolite_fg = p.metabolite.fg,
      data = frag.data,
      trace_adduct_metabolite = frag.data$trace_adduct_metabolite,
      adduct_count = adduct.count,
      fg_net_count = fg.net.count
    ))
  }

  p.all
}

