#' List all two-group network split demo candidates
#'
#' A candidate is two assignment groups within the same connectivity component,
#' each with more than `min_group_size - 1` nodes (default: more than 3 nodes).
#'
#' @noRd
.trace_list_network_split_candidates <- function(
    object,
    i.pol = 0L,
    min_group_size = 4L) {
  pol <- .trace_get_pol(i.pol)
  assign <- .trace_get_temp(object, pol, "cn.seed.assign")
  ig <- .trace_get_temp(object, pol, "cn.seed.ig.full")
  if (is.null(assign) || !nrow(assign) || is.null(ig)) {
    return(list())
  }

  eda <- as.data.frame(MSdev::edata(ig))
  vgrp <- setNames(assign$assign.group, assign$name)
  vcc <- setNames(assign$conn.component, assign$name)
  gf <- vgrp[as.character(eda$from)]
  gt <- vgrp[as.character(eda$to)]
  ccf <- vcc[as.character(eda$from)]
  cct <- vcc[as.character(eda$to)]
  cross_idx <- !is.na(gf) & !is.na(gt) & gf != gt & !is.na(ccf) & ccf == cct

  candidates <- list()
  for (cc in unique(assign$conn.component)) {
    sub <- assign[assign$conn.component == cc, , drop = FALSE]
    sz <- as.integer(table(sub$assign.group))
    if (length(sz) == 2L && all(sz >= min_group_size)) {
      g1 <- names(sz)[1]
      g2 <- names(sz)[2]
      n_cross <- sum(cross_idx & ccf == cc & ((gf == g1 & gt == g2) | (gf == g2 & gt == g1)))
      candidates[[length(candidates) + 1L]] <- list(
        conn.component = cc,
        assign.groups = names(sz),
        group_sizes = sz,
        n_cross_edges = n_cross,
        n_nodes = sum(sz),
        reason = "exact_two_group_split"
      )
    }
  }

  for (cc in unique(assign$conn.component)) {
    sub <- assign[assign$conn.component == cc, , drop = FALSE]
    sz <- table(sub$assign.group)
    big <- names(sz)[sz >= min_group_size]
    if (length(big) < 2L) next
    for (i in seq_along(big)) {
      for (j in seq_along(big)) {
        if (j <= i) next
        g1 <- big[i]
        g2 <- big[j]
        groups <- c(g1, g2)
        n_cross <- sum(
          cross_idx & ccf == cc &
            ((gf == g1 & gt == g2) | (gf == g2 & gt == g1))
        )
        candidates[[length(candidates) + 1L]] <- list(
          conn.component = cc,
          assign.groups = groups,
          group_sizes = as.integer(sz[c(g1, g2)]),
          n_cross_edges = n_cross,
          n_nodes = sum(sz[c(g1, g2)]),
          reason = "two_group_split"
        )
      }
    }
  }

  if (!length(candidates)) {
    return(list())
  }

  cand_key <- vapply(candidates, function(x) {
    paste(x$conn.component, paste(sort(x$assign.groups), collapse = "|"), sep = ":")
  }, character(1))
  candidates <- candidates[!duplicated(cand_key)]

  ord <- order(
    vapply(candidates, function(x) x$reason != "exact_two_group_split", logical(1)),
    vapply(candidates, function(x) x$n_nodes, integer(1)),
    -vapply(candidates, function(x) x$n_cross_edges, integer(1))
  )
  candidates[ord]
}

#' Compound metadata per assignment group for split-demo subtitles
#' @noRd
.trace_lookup_group_compounds <- function(object, pol, sub_assign) {
  trace_res <- object@advancedAna$TRACE[[pol]]
  if (is.null(trace_res) || !nrow(trace_res)) {
    return(NULL)
  }
  trace_res <- as.data.frame(trace_res)
  fid <- as.character(trace_res$feature_id)

  groups <- unique(sub_assign$assign.group)
  rows <- lapply(groups, function(g) {
    sub_g <- sub_assign[sub_assign$assign.group == g, , drop = FALSE]
    seed <- as.character(sub_g$assign.seed[1])
    hit <- trace_res[fid == seed, , drop = FALSE]
    if (!nrow(hit) || is.na(hit$compound_id[1])) {
      nodes <- as.character(sub_g$name)
      hit <- trace_res[
        fid %in% nodes &
          trace_res$type == "metabolite" &
          !is.na(trace_res$compound_id),
        ,
        drop = FALSE
      ]
      if (nrow(hit)) {
        hit <- hit[order(hit$TRACE_anno_score, decreasing = TRUE), , drop = FALSE]
        hit <- hit[1, , drop = FALSE]
      }
    }
    if (!nrow(hit)) {
      return(data.frame(
        assign.group = g,
        assign.seed = seed,
        compound_id = NA_character_,
        compound_name = NA_character_,
        formula = NA_character_,
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      assign.group = g,
      assign.seed = seed,
      compound_id = as.character(hit$compound_id[1]),
      compound_name = as.character(hit$name[1]),
      formula = as.character(hit$compound.formula[1]),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

#' Format assignment-group compound lines for plot subtitle
#' @noRd
.trace_format_split_demo_subtitle <- function(group_compounds) {
  if (is.null(group_compounds) || !nrow(group_compounds)) {
    return(NULL)
  }
  lines <- apply(group_compounds, 1, function(r) {
    fmt <- function(x) {
      x <- as.character(x)
      if (length(x) == 0L || is.na(x) || !nzchar(x)) "NA" else x
    }
    sprintf(
      "%s: compound_id=%s, compound_name=%s, formula=%s",
      fmt(r["assign.group"]),
      fmt(r["compound_id"]),
      fmt(r["compound_name"]),
      fmt(r["formula"])
    )
  })
  paste(lines, collapse = "  |  ")
}

#' Build demo object from a split candidate pick
#' @noRd
.trace_build_network_split_demo <- function(object, pick, demo_id = 1L, i.pol = 0L) {
  pol <- .trace_get_pol(i.pol)
  assign <- .trace_get_temp(object, pol, "cn.seed.assign")
  ig <- .trace_get_temp(object, pol, "cn.seed.ig.full")
  sub_assign <- assign[
    assign$conn.component == pick$conn.component &
      assign$assign.group %in% pick$assign.groups,
    ,
    drop = FALSE
  ]
  nodes <- as.character(sub_assign$name)
  sub_ig <- MSdev::igraph_filter_vertex(ig, nodes)
  vdata <- MSdev::vdata(sub_ig)
  vdata <- dplyr::left_join(
    vdata,
    sub_assign[, c("name", "assign.group", "assign.seed", "TRACE_net_score")],
    by = "name"
  )
  group_compounds <- .trace_lookup_group_compounds(object, pol, sub_assign)
  list(
    demo_id = demo_id,
    polarity = pol,
    conn.component = pick$conn.component,
    assign.groups = pick$assign.groups,
    group_sizes = pick$group_sizes,
    n_cross_edges = pick$n_cross_edges,
    selection_reason = pick$reason,
    group_compounds = group_compounds,
    subtitle = .trace_format_split_demo_subtitle(group_compounds),
    ig = sub_ig,
    vertex_data = vdata,
    edge_data = MSdev::edata(sub_ig),
    assign = sub_assign
  )
}

#' Find a CN network split demo within TRACE assignment results
#'
#' Returns one of the candidates from [`.trace_list_network_split_candidates()`].
#'
#' @noRd
.trace_find_network_split_demo <- function(
    object,
    i.pol = 0L,
    min_group_size = 4L,
    rank = 1L,
    exclude = list()) {
  candidates <- .trace_list_network_split_candidates(object, i.pol, min_group_size)
  if (!length(candidates)) {
    return(NULL)
  }

  is_excluded <- function(conn.component, assign.groups) {
    if (!length(exclude)) {
      return(FALSE)
    }
    any(vapply(exclude, function(ex) {
      identical(as.character(conn.component), as.character(ex$conn.component)) &&
        identical(sort(as.character(assign.groups)), sort(as.character(ex$assign.groups)))
    }, logical(1)))
  }

  if (length(exclude)) {
    keep <- vapply(candidates, function(x) {
      !is_excluded(x$conn.component, x$assign.groups)
    }, logical(1))
    candidates <- candidates[keep]
  }
  if (!length(candidates) || rank > length(candidates)) {
    return(NULL)
  }
  .trace_build_network_split_demo(object, candidates[[rank]], demo_id = rank, i.pol = i.pol)
}

#' Style edges for split-network demo plot
#' @noRd
.trace_split_demo_edge_style <- function(ig, vgrp, edge_data, group_colors) {
  ends <- igraph::ends(ig, igraph::E(ig), names = TRUE)
  eda <- as.data.frame(edge_data)
  elabs <- character(nrow(ends))
  ecols <- character(nrow(ends))
  for (i in seq_len(nrow(ends))) {
    f <- ends[i, 1]
    t <- ends[i, 2]
    row <- eda[as.character(eda$from) == f & as.character(eda$to) == t, , drop = FALSE]
    if (!nrow(row)) {
      row <- eda[as.character(eda$from) == t & as.character(eda$to) == f, , drop = FALSE]
    }
    elabs[i] <- if (nrow(row) && nzchar(row$label[1])) as.character(row$label[1]) else ""
    g1 <- vgrp[f]
    g2 <- vgrp[t]
    ecols[i] <- if (!is.na(g1) && !is.na(g2) && g1 == g2) {
      group_colors[[as.character(g1)]]
    } else {
      "grey55"
    }
  }
  list(label = elabs, color = ecols)
}

#' Build multi-line vertex labels for split-network demo plot
#' @noRd
.trace_split_demo_vertex_label <- function(
    vertex_name,
    formula = NULL,
    type = NULL,
    mz = NULL,
    id_prefix = "FT") {
  formula <- if (is.null(formula) || is.na(formula) || !nzchar(formula)) "" else as.character(formula)
  type <- if (is.null(type) || is.na(type) || !nzchar(type)) "" else as.character(type)
  mz <- if (is.null(mz) || is.na(mz)) "" else formatC(as.numeric(mz), format = "f", digits = 4)

  paste0(formula, "\n", type, "\n", mz)
}

#' Build ggplot network plot for split-network demo
#' @noRd
.trace_ggplot_network_split_demo <- function(
    ig,
    title,
    subtitle = NULL,
    layout = "fr",
    show_edge_label = TRUE,
    ...) {
  if (!requireNamespace("ggraph", quietly = TRUE) ||
      !requireNamespace("tidygraph", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "ggplot network plot requires packages ggraph, tidygraph, and ggplot2.",
      call. = FALSE
    )
  }

  g <- tidygraph::as_tbl_graph(ig)
  lay <- ggraph::create_layout(g, layout = layout)
  vpos <- as.data.frame(lay[, c("x", "y")])
  vpos$name <- as.character(lay$name)

  ends <- igraph::ends(ig, igraph::E(ig), names = TRUE)
  elabs <- igraph::E(ig)$label
  if (!isTRUE(show_edge_label)) {
    elabs <- rep("", length(elabs))
  }
  eda <- data.frame(
    from = ends[, 1],
    to = ends[, 2],
    color = igraph::E(ig)$color,
    label = elabs,
    stringsAsFactors = FALSE
  )
  eda <- merge(eda, vpos, by.x = "from", by.y = "name", suffixes = c("", ".from"))
  eda <- merge(eda, vpos, by.x = "to", by.y = "name", suffixes = c("", ".to"))
  eda$x <- (eda$x + eda$x.to) / 2
  eda$y <- (eda$y + eda$y.to) / 2
  eda$label <- ifelse(nzchar(eda$label), eda$label, NA_character_)

  p <- ggraph::ggraph(lay) +
    ggraph::geom_edge_link(
      ggplot2::aes(color = I(color)),
      arrow = grid::arrow(length = grid::unit(2, "mm"), type = "closed"),
      end_cap = ggraph::circle(4, "mm"),
      linewidth = 0.9
    ) +
    ggplot2::geom_text(
      data = eda,
      ggplot2::aes(x = x, y = y, label = label),
      size = 2.2,
      color = "grey30",
      na.rm = TRUE
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(fill = I(color)),
      shape = 21,
      color = "black",
      size = 8
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = label),
      size = 2.2,
      lineheight = 0.85
    ) +
    ggplot2::theme_void() +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 7, hjust = 0.5, lineheight = 0.95),
      ...
    )
  p
}

#' Plot a TRACE network split demo colored by assignment group
#'
#' Nodes are colored by assignment group. Within-group edges use the group
#' color; cross-group (split) edges are drawn in grey with edge labels shown.
#'
#' @param demo List returned by [`.trace_find_network_split_demo()`] or loaded
#'   from cache.
#' @param file Optional path to save the plot (pdf/png depending on extension).
#' @param group_colors Named vector of colors for assignment groups.
#' @param layout Layout passed to [ggraph::create_layout()].
#' @param ... Passed to [ggplot2::theme()].
#'
#' @returns Invisibly, a `ggplot` object suitable for [MSdev::open_plot_win()].
#' @export
plot_TRACE_network_split_demo <- function(
    demo,
    file = NULL,
    group_colors = c("#00A087", "#4DBBD5"),
    vertex_group = NULL,
    vertex_label_data = NULL,
    vertex_label_id_prefix = "FT",
    show_edge_label = TRUE,
    layout = "fr",
    ...) {
  if (is.null(demo) || is.null(demo$ig)) {
    stop("Invalid network split demo: missing igraph object.")
  }
  ig <- demo$ig
  vnames <- igraph::V(ig)$name

  grp_raw <- if (is.null(vertex_group)) {
    demo$vertex_data$assign.group[match(vnames, demo$vertex_data$name)]
  } else if (is.character(vertex_group) && length(vertex_group) == 1L &&
    vertex_group %in% names(demo$vertex_data)) {
    demo$vertex_data[[vertex_group]][match(vnames, demo$vertex_data$name)]
  } else {
    if (!is.null(names(vertex_group))) {
      unname(vertex_group[as.character(vnames)])
    } else if (length(vertex_group) == length(vnames)) {
      vertex_group
    } else {
      stop(
        "vertex_group must be NULL, a column name in demo$vertex_data, ",
        "a named vector keyed by vertex name, or a vector of same length as vertices."
      )
    }
  }

  grp <- as.character(grp_raw)
  na_idx <- which(is.na(grp_raw) | !nzchar(grp))
  if (!is.null(vertex_group) && length(na_idx)) {
    na_cols <- grDevices::hcl.colors(length(na_idx), palette = "Pastel 1")
    for (k in seq_along(na_idx)) {
      ii <- na_idx[k]
      grp[ii] <- paste0(".na.", vnames[ii])
    }
  } else if (length(na_idx)) {
    grp[na_idx] <- "Unassigned"
  }
  grp_levels <- unique(grp)

  if (is.null(names(group_colors)) || !all(grp_levels %in% names(group_colors))) {
    auto_cols <- grDevices::hcl.colors(length(grp_levels), palette = "Dark 3")
    names(auto_cols) <- grp_levels
    missing <- setdiff(grp_levels, names(group_colors))
    group_colors <- c(group_colors, auto_cols[missing])
  }
  if ("Unassigned" %in% grp_levels && is.null(group_colors[["Unassigned"]])) {
    group_colors[["Unassigned"]] <- "grey70"
  }
  if (!is.null(vertex_group) && length(na_idx)) {
    for (k in seq_along(na_idx)) {
      tag <- paste0(".na.", vnames[na_idx[k]])
      group_colors[[tag]] <- na_cols[k]
    }
  }

  igraph::V(ig)$color <- group_colors[as.character(grp)]

  if (is.null(vertex_label_data)) {
    formula <- if ("TRACE_formula" %in% names(demo$vertex_data)) {
      demo$vertex_data$TRACE_formula[match(vnames, demo$vertex_data$name)]
    } else {
      rep(NA_character_, length(vnames))
    }
    mz <- if ("mz" %in% names(demo$vertex_data)) {
      demo$vertex_data$mz[match(vnames, demo$vertex_data$name)]
    } else {
      rep(NA_real_, length(vnames))
    }
    type <- if ("type" %in% names(demo$vertex_data)) {
      demo$vertex_data$type[match(vnames, demo$vertex_data$name)]
    } else {
      rep(NA_character_, length(vnames))
    }
  } else {
    if (!is.data.frame(vertex_label_data)) {
      stop("vertex_label_data must be a data.frame.")
    }
    id_col <- if ("feature_id" %in% names(vertex_label_data)) "feature_id" else "name"
    if (!id_col %in% names(vertex_label_data)) {
      stop("vertex_label_data must contain `feature_id` or `name`.")
    }
    key <- as.character(vertex_label_data[[id_col]])
    formula <- vertex_label_data$formula[match(as.character(vnames), key)]
    type <- vertex_label_data$type[match(as.character(vnames), key)]
    mz <- vertex_label_data$mz[match(as.character(vnames), key)]
  }

  igraph::V(ig)$label <- mapply(
    .trace_split_demo_vertex_label,
    vertex_name = as.character(vnames),
    formula = formula,
    type = type,
    mz = mz,
    MoreArgs = list(id_prefix = vertex_label_id_prefix),
    USE.NAMES = FALSE
  )

  vgrp <- stats::setNames(grp, as.character(vnames))
  estyle <- .trace_split_demo_edge_style(ig, vgrp, demo$edge_data, group_colors)
  igraph::E(ig)$color <- estyle$color
  igraph::E(ig)$label <- estyle$label

  title <- sprintf(
    "TRACE assignment split (component %s: %s)",
    demo$conn.component,
    paste(unique(demo$assign.groups), collapse = " vs ")
  )
  subtitle <- if (!is.null(demo$subtitle) && nzchar(demo$subtitle)) {
    demo$subtitle
  } else if (!is.null(demo$group_compounds)) {
    .trace_format_split_demo_subtitle(demo$group_compounds)
  } else {
    NULL
  }

  p <- .trace_ggplot_network_split_demo(
    ig,
    title = title,
    subtitle = subtitle,
    layout = layout,
    show_edge_label = show_edge_label,
    ...
  )

  if (!is.null(file)) {
    ext <- tolower(tools::file_ext(file))
    if (ext == "pdf") {
      ggplot2::ggsave(file, p, width = 9, height = 7, device = "pdf")
    } else if (ext %in% c("png", "jpg", "jpeg")) {
      ggplot2::ggsave(file, p, width = 1400 / 120, height = 1000 / 120, dpi = 120)
    } else {
      ggplot2::ggsave(file, p, width = 9, height = 7)
    }
  }

  invisible(p)
}

#' Build and cache a TRACE CN network split demo
#'
#' Searches assignment results for a connectivity component that splits into
#' two subnetworks (each with more than three nodes when available), plots the
#' subgraph colored by assignment group, and saves the demo as an RDS list.
#'
#' @param object MSdev object with TRACE network assignment results.
#' @param i.pol Polarity index (`0` = negative, `1` = positive).
#' @param demo_id Demo index (`1`, `2`, ...). Files are saved as
#'   `cache/demo{n}.rds` and `cache/demo{n}.pdf`.
#' @param cache_dir Directory for cached demo files.
#' @param cache_file Optional explicit RDS output path (overrides `demo_id`).
#' @param plot_file Optional explicit plot output path (overrides `demo_id`).
#' @param min_group_size Minimum nodes per assignment group.
#' @param rank Alias for `demo_id` when selecting the split candidate.
#'
#' @returns Invisibly, the demo list.
#' @export
TRACE_cache_network_split_demo <- function(
    object,
    i.pol = 0L,
    demo_id = 1L,
    cache_dir = NULL,
    cache_file = NULL,
    plot_file = NULL,
    min_group_size = 4L,
    rank = 1L,
    exclude = list()) {
  cache_dir <- if (is.null(cache_dir)) {
    file.path(Sys.getenv("TRACE_PKG_ROOT", unset = getwd()), "cache")
  } else {
    cache_dir
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  if (demo_id > 1L && !length(exclude)) {
    for (i in seq_len(demo_id - 1L)) {
      prev_file <- file.path(cache_dir, paste0("demo", i, ".rds"))
      if (!file.exists(prev_file)) {
        stop("Build demo", i, " first (missing ", prev_file, ").")
      }
      prev <- readRDS(prev_file)
      exclude[[length(exclude) + 1L]] <- list(
        conn.component = prev$conn.component,
        assign.groups = prev$assign.groups
      )
    }
  }

  demo <- .trace_find_network_split_demo(
    object,
    i.pol = i.pol,
    min_group_size = min_group_size,
    rank = rank,
    exclude = exclude
  )
  if (is.null(demo)) {
    stop(
      "No suitable network split demo #", rank,
      " found for polarity ", .trace_get_pol(i.pol), "."
    )
  }
  demo$demo_id <- demo_id

  if (is.null(cache_file)) {
    cache_file <- file.path(cache_dir, paste0("demo", demo_id, ".rds"))
  }
  if (is.null(plot_file)) {
    plot_file <- file.path(cache_dir, paste0("demo", demo_id, ".pdf"))
  }

  saveRDS(demo, cache_file)
  plot_TRACE_network_split_demo(demo, file = plot_file)
  MSdev::message_with_time(
    "Cached TRACE network split demo: ",
    normalizePath(cache_file, winslash = "/", mustWork = FALSE)
  )
  MSdev::message_with_time(
    "Saved plot: ",
    normalizePath(plot_file, winslash = "/", mustWork = FALSE)
  )
  invisible(demo)
}

#' Export all network split demos to cache
#'
#' Finds every connectivity-component split into two assignment groups with
#' more than three nodes per group, then writes `demo1.rds`, `demo1.pdf`, ...
#' plus a `demo_index.rds` manifest under `cache_dir`.
#'
#' @param object MSdev object with TRACE network assignment results.
#' @param i.pol Polarity index (`0` = negative, `1` = positive).
#' @param cache_dir Directory for cached demo files.
#' @param min_group_size Minimum nodes per assignment group (`4` = more than 3).
#'
#' @returns Invisibly, a list of demo objects.
#' @export
TRACE_cache_all_network_split_demos <- function(
    object,
    i.pol = 0L,
    cache_dir = NULL,
    min_group_size = 4L) {
  cache_dir <- if (is.null(cache_dir)) {
    file.path(Sys.getenv("TRACE_PKG_ROOT", unset = getwd()), "cache")
  } else {
    cache_dir
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  candidates <- .trace_list_network_split_candidates(object, i.pol, min_group_size)
  if (!length(candidates)) {
    stop(
      "No two-group network split demos found for polarity ",
      .trace_get_pol(i.pol), " (each group must have >3 nodes)."
    )
  }

  old_files <- list.files(
    cache_dir,
    pattern = "^demo[0-9]+\\.(rds|pdf)$|^demo_index\\.rds$",
    full.names = TRUE
  )
  if (length(old_files)) {
    unlink(old_files)
  }

  demos <- vector("list", length(candidates))
  index <- data.frame(
    demo_id = seq_along(candidates),
    conn.component = vapply(candidates, function(x) x$conn.component, numeric(1)),
    assign.groups = vapply(candidates, function(x) {
      paste(x$assign.groups, collapse = " vs ")
    }, character(1)),
    group_sizes = vapply(candidates, function(x) {
      paste(x$group_sizes, collapse = " vs ")
    }, character(1)),
    n_nodes = vapply(candidates, function(x) x$n_nodes, integer(1)),
    n_cross_edges = vapply(candidates, function(x) x$n_cross_edges, integer(1)),
    rds = file.path(cache_dir, paste0("demo", seq_along(candidates), ".rds")),
    pdf = file.path(cache_dir, paste0("demo", seq_along(candidates), ".pdf")),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(candidates)) {
    demo <- .trace_build_network_split_demo(object, candidates[[i]], demo_id = i, i.pol = i.pol)
    saveRDS(demo, index$rds[i])
    plot_TRACE_network_split_demo(demo, file = index$pdf[i])
    demos[[i]] <- demo
    MSdev::message_with_time(
      "Cached demo", i, "/", length(candidates), ": ",
      normalizePath(index$rds[i], winslash = "/", mustWork = FALSE)
    )
  }

  index$subtitle <- vapply(demos, function(d) {
    if (is.null(d$subtitle) || !nzchar(d$subtitle)) "" else d$subtitle
  }, character(1))

  saveRDS(
    list(
      polarity = .trace_get_pol(i.pol),
      min_group_size = min_group_size,
      n_demo = length(candidates),
      index = index,
      demos = demos
    ),
    file.path(cache_dir, "demo_index.rds")
  )
  MSdev::message_with_time(
    "Exported ", length(candidates), " network split demos to ",
    normalizePath(cache_dir, winslash = "/", mustWork = FALSE)
  )
  invisible(demos)
}

#' Load a cached TRACE network split demo
#'
#' @param demo_id Demo index (`1`, `2`, ...).
#' @param cache_dir Directory containing `demo{n}.rds`.
#'
#' @returns Demo list.
#' @export
TRACE_load_network_split_demo <- function(
    demo_id = 1L,
    cache_dir = file.path("cache")) {
  cache_file <- file.path(cache_dir, paste0("demo", demo_id, ".rds"))
  if (!file.exists(cache_file)) {
    stop("Cache file not found: ", cache_file)
  }
  readRDS(cache_file)
}
