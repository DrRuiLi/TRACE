# PAVE-style global network assignment (internal helpers)
#
# Spec defaults (no Matlab source in repo):
# - Graph: full MS candidate net (CN_label + adduct + isotope + fragment), vertices = CN seeds
# - Paths: all simple paths, capped at max_path_length (default 6)
# - Blocking: net accumulated chemform (non-iso, simplified) in MS mass-diff lookup -> score 0
# - Connection score: max(0, 1 - sum|mz.ppm|/ppm.dyn - sum|rt.diff|/rt.tol.dyn) on best non-blocked path
# - Clustering: connected components on pairs with score >= connection_cutoff (default 0.5)

#' Dynamic ppm / RT tolerances from TRACE_temp
#' @keywords internal
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

#' MS-derived chemform lookup for path blocking
#' @keywords internal
.trace_ms_mass_diff_chemforms <- function(i.pol) {
  ad <- get_adduct_mass_diff(polarity = i.pol)$chemform_diff
  fg <- get_fragment_mass_diff()$chemform_diff
  is <- get_iso_mass_diff()$chemform_diff
  unique(chemform_simplify(c(ad, fg, is)))
}

#' Build deduplicated candidate edge table (legacy PAVE integration)
#' @keywords internal
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
#' @keywords internal
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

#' Ordered edge rows along a vertex path
#' @keywords internal
.trace_path_edge_rows <- function(eda, vpath) {
  vpath <- as.character(vpath)
  out <- vector("list", length(vpath) - 1L)
  for (k in seq_len(length(vpath) - 1L)) {
    f <- vpath[k]
    t <- vpath[k + 1L]
    fwd <- eda[as.character(eda$from) == f & as.character(eda$to) == t, , drop = FALSE]
    if (nrow(fwd)) {
      out[[k]] <- list(row = fwd[1, , drop = FALSE], dir = 1L)
    } else {
      rev <- eda[as.character(eda$from) == t & as.character(eda$to) == f, , drop = FALSE]
      if (!nrow(rev)) {
        return(NULL)
      }
      out[[k]] <- list(row = rev[1, , drop = FALSE], dir = -1L)
    }
  }
  out
}

#' Score one simple path; blocked if net chemform matches MS mass-diff table
#' @keywords internal
.trace_score_path <- function(eda, vpath, blocked_chemforms, ppm.dyn, rt.tol.dyn) {
  if (length(vpath) < 2L) {
    return(list(score = 1, blocked = FALSE, net_chemform = ""))
  }
  rows <- .trace_path_edge_rows(eda, vpath)
  if (is.null(rows)) {
    return(list(score = 0, blocked = TRUE, net_chemform = NA_character_))
  }
  edge_df <- dplyr::bind_rows(lapply(rows, function(x) x$row))
  dirs <- vapply(rows, function(x) x$dir, integer(1))
  cds <- edge_df$chemform_diff
  cds[is.na(cds)] <- ""
  cds <- MSCC::chemform_multi(cds, dirs, return = "chemform")
  net_cd <- MSCC:::chemform_sum(cds)
  net_cd <- chemform_remove_iso(net_cd)
  net_cd <- chemform_simplify(net_cd)
  blocked <- nzchar(net_cd) && net_cd %in% blocked_chemforms
  if (blocked) {
    return(list(score = 0, blocked = TRUE, net_chemform = net_cd))
  }
  mz_sum <- sum(abs(edge_df$mz.ppm), na.rm = TRUE)
  rt_sum <- sum(abs(edge_df$rt.diff), na.rm = TRUE)
  score <- max(0, 1 - mz_sum / ppm.dyn - rt_sum / rt.tol.dyn)
  list(score = score, blocked = FALSE, net_chemform = net_cd)
}

#' Best connection score between two nodes (max over simple paths)
#' @keywords internal
.trace_pair_connection_score <- function(
    ig,
    from,
    to,
    blocked_chemforms,
    ppm.dyn,
    rt.tol.dyn,
    max_path_length = 6L) {
  from <- as.character(from)
  to <- as.character(to)
  if (from == to) {
    return(list(score = 1, blocked = FALSE))
  }
  eda <- MSdev::edata(ig)
  vpaths <- igraph::all_simple_paths(
    ig,
    from = from,
    to = to,
    mode = "all",
    cutoff = max_path_length
  )
  if (length(vpaths) == 0L) {
    return(list(score = 0, blocked = TRUE))
  }
  best <- 0
  any_ok <- FALSE
  for (vp in vpaths) {
    vp <- as.character(vp)
    res <- .trace_score_path(eda, vp, blocked_chemforms, ppm.dyn, rt.tol.dyn)
    if (!res$blocked) {
      any_ok <- TRUE
      best <- max(best, res$score)
    }
  }
  if (!any_ok) {
    return(list(score = 0, blocked = TRUE))
  }
  list(score = best, blocked = FALSE)
}

#' Global assignment within one connectivity component
#' @keywords internal
.trace_global_assign_component <- function(
    ig,
    nodes,
    blocked_chemforms,
    ppm.dyn,
    rt.tol.dyn,
    max_path_length = 6L,
    connection_cutoff = 0.5) {
  nodes <- as.character(nodes)
  n <- length(nodes)
  if (n == 1L) {
    return(data.frame(
      name = nodes,
      assign.group = 1L,
      assign.seed = nodes,
      stringsAsFactors = FALSE
    ))
  }
  score_mat <- matrix(0, n, n, dimnames = list(nodes, nodes))
  blocked_mat <- matrix(FALSE, n, n, dimnames = list(nodes, nodes))
  for (i in seq_len(n)) {
    score_mat[i, i] <- 1
    for (j in seq_len(n)) {
      if (i >= j) next
      sc <- .trace_pair_connection_score(
        ig, nodes[i], nodes[j],
        blocked_chemforms, ppm.dyn, rt.tol.dyn, max_path_length
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
  data.frame(
    name = nodes,
    assign.group = assign.group,
    assign.seed = assign.seed,
    stringsAsFactors = FALSE
  )
}

#' Run PAVE-style global assignment for all CN seeds
#' @keywords internal
.trace_run_global_assignment <- function(
    cn.seed.ig.full,
    cn.seed,
    i.pol,
    object,
    pol,
    max_path_length = 6L,
    connection_cutoff = 0.5) {
  tol <- .trace_dyn_tolerance(object, pol)
  blocked_chemforms <- .trace_ms_mass_diff_chemforms(i.pol)
  conn_comp <- igraph::components(cn.seed.ig.full)$membership
  cn.seed <- as.character(unique(cn.seed))
  assign_list <- lapply(split(cn.seed, conn_comp[cn.seed]), function(nodes) {
    sub_ig <- MSdev::igraph_filter_vertex(cn.seed.ig.full, nodes)
    .trace_global_assign_component(
      sub_ig,
      nodes,
      blocked_chemforms,
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
