# PAVE-style global network assignment (internal helpers)
#
# Spec defaults (no Matlab source in repo):
# - Graph: full MS candidate net (CN_label + adduct + isotope + fragment), vertices = CN seeds
# - Paths: all simple paths, capped at max_path_length (default 6)
# - Blocking: cumulative adduct chemform along path (fragment/isotope/CN edges ignored);
#   blocked when non-zero net adduct change is not in the adduct mass-diff table
# - Connection score: max(0, 1 - sum|mz.ppm|/ppm.dyn - sum|rt.diff|/rt.tol.dyn) on best non-blocked path
# - Clustering: connected components on pairs with score >= connection_cutoff (default 0.5)

#' Dynamic ppm / RT tolerances from TRACE_temp
#' @noRd
.trace_dyn_tolerance <- function(object, pol) {
  mz.fit <- .trace_get_temp(object, pol, "mz.dyn")
  rt.fit <- .trace_get_temp(object, pol, "rt.dyn")
  ppm <- .trace_get_temp(object, pol, "ppm")
  rt.tol <- .trace_get_temp(object, pol, "rt.tol")
  ppm.dyn <- if (!is.null(mz.fit) && is.finite(mz.fit$sd)) {
    mz.fit$sd * stats::qnorm(0.999)
  } else if (!is.null(ppm)) {
    ppm
  } else {
    5
  }
  rt.tol.dyn <- if (!is.null(rt.fit) && is.finite(rt.fit$sd)) {
    rt.fit$sd * stats::qnorm(0.99999)
  } else if (!is.null(rt.tol)) {
    rt.tol
  } else {
    10
  }
  list(ppm.dyn = ppm.dyn, rt.tol.dyn = rt.tol.dyn)
}

#' Adduct chemform lookup for path blocking
#' @noRd
.trace_adduct_chemforms <- function(i.pol) {
  unique(chemform_simplify(get_adduct_mass_diff(polarity = i.pol)$chemform_diff))
}

#' Build deduplicated candidate edge table (legacy PAVE integration)
#' @noRd
.trace_build_candidate_net <- function(cn.net.hit, ad.net, is.net, fg.net) {
  cn.cn <- data.table::copy(cn.net.hit)
  cn.cn[, type := "CN_label"]
  xcms.net.candidate <- data.table::rbindlist(
    list(ad.net, is.net, fg.net, cn.cn),
    use.names = TRUE,
    fill = TRUE
  )
  xcms.net.candidate[, chemform_diff := chemform_simplify(chemform_diff)]
  xcms.net.candidate[, is_CN := (type == "CN_label")]
  xcms.net.candidate[, temp := data.table::fcase(
    is_CN, TRACE_pattern,
    type == "adduct", paste0(adduct.from, ">>", adduct.to),
    default = chemform_diff
  )]
  xcms.net.candidate[, label := paste0(type, ": ", temp)]
  xcms.net.candidate[, temp := factor(type, levels = c("CN_label", "fragment", "adduct", "isotope"))]
  data.table::setorder(xcms.net.candidate, temp)
  xcms.net.candidate <- xcms.net.candidate[, .SD[1], by = .(ion1, chemform_diff)]
  data.table::setcolorder(xcms.net.candidate, c("from", "to"))
  xcms.net.candidate[, eid := seq_len(.N)]
  xcms.net.candidate
}

#' Per-edge CN consistency retyping (seed subgraph)
#' @noRd
.trace_retype_seed_edges <- function(cn.seed.net, cn.seed.formula) {
  .equal.cn.chemform <- function(cn.diff, chemfrom.diff) {
    m <- MSCC::chemform_parse(c(cn.diff, chemfrom.diff))
    m <- MSdev::get_matrix_value_fill_with_NA(m, colnames_vec = c("C", "N"))
    m[is.na(m)] <- 0
    m1 <- m[seq_along(cn.diff), ]
    m2 <- m[length(cn.diff) + seq_along(cn.diff), ]
    unname(m1[, "C"] == m2[, "C"] & m1[, "N"] == m2[, "N"])
  }

  cn.seed.net %>%
    dplyr::mutate(
      chemform_diff = chemform_simplify(chemform_diff),
      label = dplyr::case_when(
        type == "adduct" ~ paste0(adduct.from, ">>", adduct.to),
        TRUE ~ chemform_diff
      ),
      label = paste0(type, ": ", label),
      from.cn = cn.seed.formula[as.character(from)],
      to.cn = cn.seed.formula[as.character(to)],
      cn.diff = MSCC::chemform_calc(to.cn, from.cn, "-", return = "chemform"),
      cn.equal = cn.diff == "",
      chemform.equal = .equal.cn.chemform(cn.diff, chemform_diff),
      new.type = dplyr::case_when(
        type == "CN_label" ~ "CN_label",
        type == "adduct" & cn.equal ~ "adduct",
        type == "adduct" & !cn.equal & chemform.equal ~ "false",
        type == "adduct" & !cn.equal & !chemform.equal ~ "false",
        type == "fragment" & chemform.equal ~ "fragment",
        type == "fragment" & !chemform.equal ~ "false",
        type == "isotope" & chemform.equal ~ "isotope",
        type == "isotope" & !chemform.equal & element == "[13]C" & cn.diff == "C-2" ~ "isotope",
        TRUE ~ "false"
      ),
      old.type = type,
      type = new.type
    ) %>%
    dplyr::filter(type != "false") %>%
    dplyr::mutate(eid = dplyr::row_number())
}

#' Pre-parse chemform_diff strings and valid adduct lookup into numeric matrices.
#' Returns a list with:
#'   $cd_mat  - named numeric matrix (rows = unique chemform strings, cols = elements)
#'   $valid_adduct_rows - row indices in cd_mat for known adduct formula changes
#'   $eid_from, $eid_to, $eid_cd_row, $eid_is_adduct, $eid_mzppm, $eid_rtdiff
#' @noRd
.trace_build_score_cache <- function(eda, valid_adduct_chemforms) {
  cds <- as.character(eda$chemform_diff)
  cds[is.na(cds)] <- ""
  # parse all unique chemforms at once (single MSCC call)
  all_cds <- unique(c(cds, valid_adduct_chemforms, ""))
  raw_mat <- MSCC::chemform_parse(all_cds)
  # Strip isotope columns: fold their counts into the base element columns so
  # net-path comparison is on a simplified (non-isotope) element basis.
  is_isotope_name <- function(x) {
    grepl("^\\[[0-9]+\\]", x)
  }
  get_uniso_name <- function(x) {
    sub("^\\[[0-9]+\\]", "", x)
  }
  iso_names <- colnames(raw_mat)[is_isotope_name(colnames(raw_mat))]
  if (length(iso_names)) {
    base_names <- get_uniso_name(iso_names)
    missing_base <- setdiff(unique(base_names), colnames(raw_mat))
    if (length(missing_base)) {
      add <- matrix(0, nrow(raw_mat), length(missing_base),
                    dimnames = list(rownames(raw_mat), missing_base))
      raw_mat <- cbind(raw_mat, add)
    }
    for (k in seq_along(iso_names)) {
      raw_mat[, base_names[k]] <- raw_mat[, base_names[k]] + raw_mat[, iso_names[k]]
    }
    raw_mat <- raw_mat[, !is_isotope_name(colnames(raw_mat)), drop = FALSE]
  }
  rownames(raw_mat) <- all_cds
  valid_adduct_rows <- which(all_cds %in% valid_adduct_chemforms)
  eid_cd_row <- match(cds, all_cds)
  eid_type <- as.character(eda$type)
  list(
    cd_mat              = raw_mat,
    all_cds             = all_cds,
    valid_adduct_rows   = valid_adduct_rows,
    eid_cd_row          = eid_cd_row,
    eid_is_adduct       = eid_type == "adduct",
    eid_from            = as.character(eda$from),
    eid_to              = as.character(eda$to),
    eid_mzppm           = abs(as.numeric(eda$mz.ppm)),
    eid_rtdiff          = abs(as.numeric(eda$rt.diff))
  )
}

#' Score one simple path using precomputed numeric cache (no per-path MSCC calls).
#' @noRd
.trace_score_path_fast <- function(vpath, cache, ppm.dyn, rt.tol.dyn) {
  n <- length(vpath)
  if (n < 2L) return(list(score = 1, blocked = FALSE))
  # gather edge rows and directions
  edge_idx  <- integer(n - 1L)
  dirs      <- integer(n - 1L)
  for (k in seq_len(n - 1L)) {
    f <- vpath[k]; t <- vpath[k + 1L]
    fwd <- which(cache$eid_from == f & cache$eid_to == t)
    if (length(fwd)) {
      edge_idx[k] <- fwd[1L]; dirs[k] <- 1L
    } else {
      rev <- which(cache$eid_from == t & cache$eid_to == f)
      if (!length(rev)) return(list(score = 0, blocked = TRUE))
      edge_idx[k] <- rev[1L]; dirs[k] <- -1L
    }
  }
  # cumulative adduct form change: only adduct edges contribute
  cd_rows <- cache$eid_cd_row[edge_idx]
  adduct_wt <- dirs * as.integer(cache$eid_is_adduct[edge_idx])
  net_adduct_vec <- colSums(cache$cd_mat[cd_rows, , drop = FALSE] * adduct_wt)
  if (any(net_adduct_vec != 0)) {
    in_adduct_table <- any(vapply(cache$valid_adduct_rows, function(r) {
      all(cache$cd_mat[r, ] == net_adduct_vec)
    }, logical(1)))
    if (!in_adduct_table) return(list(score = 0, blocked = TRUE))
  }
  mz_sum <- sum(cache$eid_mzppm[edge_idx], na.rm = TRUE)
  rt_sum <- sum(cache$eid_rtdiff[edge_idx], na.rm = TRUE)
  score  <- max(0, 1 - mz_sum / ppm.dyn - rt_sum / rt.tol.dyn)
  list(score = score, blocked = FALSE)
}

#' Best connection score between two nodes (max over simple paths)
#' @noRd
.trace_pair_connection_score <- function(
    ig,
    from,
    to,
    cache,
    ppm.dyn,
    rt.tol.dyn,
    max_path_length = 5L) {
  from <- as.character(from)
  to   <- as.character(to)
  if (from == to) return(list(score = 1, blocked = FALSE))
  vpaths <- igraph::all_simple_paths(ig, from = from, to = to,
                                     mode = "all", cutoff = max_path_length)
  if (!length(vpaths)) return(list(score = 0, blocked = TRUE))
  best <- 0; any_ok <- FALSE
  for (vp in vpaths) {
    vp  <- igraph::as_ids(vp)
    res <- .trace_score_path_fast(vp, cache, ppm.dyn, rt.tol.dyn)
    if (!res$blocked) { any_ok <- TRUE; best <- max(best, res$score) }
  }
  if (!any_ok) return(list(score = 0, blocked = TRUE))
  list(score = best, blocked = FALSE)
}

#' Global assignment within one connectivity component
#' @noRd
.trace_global_assign_component <- function(
    ig,
    nodes,
    valid_adduct_chemforms,
    ppm.dyn,
    rt.tol.dyn,
    max_path_length = 5L,
    connection_cutoff = 0.5) {
  nodes <- as.character(nodes)
  n <- length(nodes)
  if (n == 1L) {
    return(data.frame(
      name = nodes,
      assign.group = 1L,
      assign.seed = nodes,
      TRACE_net_score = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  # Build numeric chemform cache once per component (single MSCC parse call).
  cache <- .trace_build_score_cache(MSdev::edata(ig), valid_adduct_chemforms)
  score_mat <- matrix(0, n, n, dimnames = list(nodes, nodes))
  blocked_mat <- matrix(FALSE, n, n, dimnames = list(nodes, nodes))
  for (i in seq_len(n)) {
    score_mat[i, i] <- 1
    for (j in seq_len(n)) {
      if (i >= j) next
      sc <- .trace_pair_connection_score(
        ig, nodes[i], nodes[j],
        cache, ppm.dyn, rt.tol.dyn, max_path_length
      )
      score_mat[i, j] <- sc$score
      score_mat[j, i] <- sc$score
      blocked_mat[i, j] <- sc$blocked
      blocked_mat[j, i] <- sc$blocked
    }
  }
  adj <- score_mat >= connection_cutoff & !blocked_mat
  diag(adj) <- TRUE
  g_compat <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = TRUE)
  memb <- igraph::components(g_compat)$membership
  assign.group <- as.integer(memb[nodes])
  assign.seed <- vapply(split(nodes, assign.group), function(ns) {
    as.character(min(as.numeric(ns)))
  }, character(1))
  assign.seed <- assign.seed[as.character(assign.group)]
  # TRACE_net_score: best connection score from each node to any *other* member
  # of its connectivity sub-network (this whole component). The diagonal
  # self-term (score 1) is excluded. A single-node component has no partner and
  # is handled by the n == 1 branch above (NA).
  TRACE_net_score <- vapply(seq_len(n), function(i) {
    partners <- seq_len(n)[-i]
    max(score_mat[i, partners])
  }, numeric(1))
  data.frame(
    name = nodes,
    assign.group = assign.group,
    assign.seed = assign.seed,
    TRACE_net_score = TRACE_net_score,
    stringsAsFactors = FALSE
  )
}

#' Run PAVE-style global assignment for all CN seeds
#' @noRd
.trace_run_global_assignment <- function(
    cn.seed.ig.full,
    cn.seed,
    i.pol,
    object,
    pol,
    max_path_length = 5L,
    connection_cutoff = 0.5) {
  tol <- .trace_dyn_tolerance(object, pol)
  valid_adduct_chemforms <- .trace_adduct_chemforms(i.pol)
  conn_comp <- igraph::components(cn.seed.ig.full)$membership
  cn.seed <- as.character(unique(cn.seed))
  assign_list <- lapply(split(cn.seed, conn_comp[cn.seed]), function(nodes) {
    sub_ig <- MSdev::igraph_filter_vertex(cn.seed.ig.full, nodes)
    .trace_global_assign_component(
      sub_ig,
      nodes,
      valid_adduct_chemforms,
      tol$ppm.dyn,
      tol$rt.tol.dyn,
      max_path_length,
      connection_cutoff
    )
  })
  assign_df <- dplyr::bind_rows(assign_list)
  assign_df$conn.component <- conn_comp[assign_df$name]
  assign_df$assign.group <- paste0(assign_df$conn.component, ".", assign_df$assign.group)
  assign_df
}

#' Recompute path-compatibility score matrix for one component
#'
#' This is a debugging/inspection helper that reproduces the internal
#' `score_mat` and `blocked_mat` used to derive `assign.group` within a
#' connectivity component (`conn.component`).
#'
#' @param object MSdev object with TRACE network assignment results.
#' @param i.pol Polarity index (`0` = negative, `1` = positive).
#' @param conn.component Numeric connectivity component id (e.g. `638`).
#' @param max_path_length Maximum simple path length used when scoring (default `5`).
#' @param connection_cutoff Compatibility cutoff used to build the adjacency matrix (default `0.5`).
#'
#' @return A list with `nodes`, `assign` (subset of `cn.seed.assign`),
#'   `score_mat`, `blocked_mat`, `ppm.dyn`, `rt.tol.dyn`, `connection_cutoff`,
#'   and `max_path_length`.
#' @export
TRACE_score_mat_component <- function(
    object,
    i.pol = 0L,
    conn.component,
    max_path_length = 5L,
    connection_cutoff = 0.5) {
  pol <- .trace_get_pol(i.pol)
  assign <- .trace_get_temp(object, pol, "cn.seed.assign")
  ig_full <- .trace_get_temp(object, pol, "cn.seed.ig.full")
  if (is.null(assign) || !nrow(assign) || is.null(ig_full)) {
    stop("Missing TRACE network assignment results in TRACE_temp for ", pol, ".")
  }
  if (missing(conn.component) || length(conn.component) != 1L) {
    stop("`conn.component` must be a single value.")
  }
  conn.component <- as.numeric(conn.component)
  sub_assign <- assign[assign$conn.component == conn.component, , drop = FALSE]
  if (!nrow(sub_assign)) {
    stop("No nodes found for conn.component = ", conn.component, " (", pol, ").")
  }
  nodes <- as.character(sub_assign$name)
  ig <- MSdev::igraph_filter_vertex(ig_full, nodes)

  tol <- .trace_dyn_tolerance(object, pol)
  valid_adduct_chemforms <- .trace_adduct_chemforms(i.pol)
  cache <- .trace_build_score_cache(MSdev::edata(ig), valid_adduct_chemforms)

  n <- length(nodes)
  score_mat <- matrix(0, n, n, dimnames = list(nodes, nodes))
  blocked_mat <- matrix(FALSE, n, n, dimnames = list(nodes, nodes))
  for (i in seq_len(n)) {
    score_mat[i, i] <- 1
    for (j in seq_len(n)) {
      if (i >= j) next
      sc <- .trace_pair_connection_score(
        ig,
        from = nodes[i],
        to = nodes[j],
        cache = cache,
        ppm.dyn = tol$ppm.dyn,
        rt.tol.dyn = tol$rt.tol.dyn,
        max_path_length = max_path_length
      )
      score_mat[i, j] <- sc$score
      score_mat[j, i] <- sc$score
      blocked_mat[i, j] <- sc$blocked
      blocked_mat[j, i] <- sc$blocked
    }
  }

  list(
    polarity = pol,
    conn.component = conn.component,
    nodes = nodes,
    assign = sub_assign,
    score_mat = score_mat,
    blocked_mat = blocked_mat,
    ppm.dyn = tol$ppm.dyn,
    rt.tol.dyn = tol$rt.tol.dyn,
    connection_cutoff = connection_cutoff,
    max_path_length = max_path_length
  )
}
