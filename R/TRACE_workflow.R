#' Get polarity label for TRACE workflow
#' @noRd
.trace_get_pol <- function(i.pol) {
  ifelse(i.pol == 0, "Negative", "Positive")
}

#' Initialize TRACE_temp container on MSdev object
#' @noRd
.trace_init_temp <- function(object) {
  if (is.null(object@advancedAna$TRACE_temp)) {
    object@advancedAna$TRACE_temp <- list()
  }
  object
}

#' Read a value from TRACE_temp
#' @noRd
.trace_get_temp <- function(object, pol, key) {
  object@advancedAna$TRACE_temp[[pol]][[key]]
}

#' Write a value to TRACE_temp
#' @noRd
.trace_set_temp <- function(object, pol, key, value) {
  object@advancedAna$TRACE_temp[[pol]][[key]] <- value
  object
}

#' Build initial TRACE CN/adduct/isotope/fragment candidate networks
#'
#' @param object MSdev object.
#' @param i.pol Polarity index. `0` for negative and `1` for positive.
#' @param rt.tol RT tolerance in seconds.
#' @param ppm m/z tolerance in ppm.
#' @param ratio.adjust Numeric vector of length 4 used to adjust the
#'   theoretical CN labelling ratio across `S12C14N`, `S12C15N`,
#'   `S13C14N`, `S13C15N`.
#' @param TRACE_cor_cutoff Minimum `TRACE_cor` to retain a CN hit.
#'
#' @returns Updated MSdev object with intermediate data in
#'   `object@advancedAna$TRACE_temp[[pol]]`.
#' @export
TRACE_get_CN_net <- function(
    object,
    i.pol,
    rt.tol = 10,
    ppm = 5,
    ratio.adjust = c(1, 1, 1, 1),
    TRACE_cor_cutoff = 0.75) {
  object <- .trace_init_temp(object)

  pol <- .trace_get_pol(i.pol)
  MSdev::message_with_time(pol)

  xcms.xcms <- object@xcmsData[[paste0(pol, "MS1")]]
  if (is.null(xcms.xcms)) {
    return(object)
  }

  xcms.xcms <- MSdev::xcms_get_feature_wmean(xcms.xcms)
  xcms.net <- get_xcms_feature_connect(xcms.xcms, rt.tol = rt.tol)
  xcms.fdf <- xcms::featureDefinitions(xcms.xcms)
  xcms.val <- xcms::featureValues(xcms.xcms, missing = 0, value = "maxo")
  xcms.TRACE.sample <- Biobase::pData(xcms.xcms) %>%
    dplyr::filter(sample.type %in% c("S12C14N", "S12C15N", "S13C14N", "S13C15N", "Blank"))

  cn.mass.diff <- get_CN_mass_diff_table(C_max = 99,N_max = 10)[, type := "CN_label"]
  ad.mass.diff <- get_adduct_mass_diff(polarity = i.pol)[, type := "adduct"]
  is.mass.diff <- get_iso_mass_diff()[, type := "isotope"]
  fg.mass.diff <- get_fragment_mass_diff()[, type := "fragment"]

  mass.diff.range <- range(
    cn.mass.diff$mass_diff,
    ad.mass.diff$mass_diff,
    is.mass.diff$mass_diff
  )
  xcms.net <- xcms.net[data.table::between(mz.diff, mass.diff.range[1], mass.diff.range[2])]

  cn.match <- MSdev::match_mz_foverlaps(
    xcms.net$mz.diff, cn.mass.diff$mass_diff,
    ppm.base = xcms.net$mz.mean, ppm = ppm
  )
  ad.match <- MSdev::match_mz_foverlaps(
    xcms.net$mz.diff, ad.mass.diff$mass_diff,
    ppm.base = xcms.net$mz.mean, ppm = ppm
  )
  is.match <- MSdev::match_mz_foverlaps(
    xcms.net$mz.diff, is.mass.diff$mass_diff,
    ppm.base = xcms.net$mz.mean, ppm = ppm
  )
  fg.match <- MSdev::match_mz_foverlaps(
    xcms.net$mz.diff, fg.mass.diff$mass_diff,
    ppm.base = xcms.net$mz.mean, ppm = ppm
  )

  cn.net <- cbind(
    xcms.net[cn.match$ion1, ], cn.match[, c("mz.ppm", "ion1")],
    cn.mass.diff[cn.match$ion2, ]
  )[mass_diff > 0]
  ad.net <- cbind(
    xcms.net[ad.match$ion1, ], ad.match[, c("mz.ppm", "ion1")],
    ad.mass.diff[ad.match$ion2, ]
  )[mass_diff > 0]
  is.net <- cbind(
    xcms.net[is.match$ion1, ], is.match[, c("mz.ppm", "ion1")],
    is.mass.diff[is.match$ion2, ]
  )[mass_diff > 0]
  fg.net <- cbind(
    xcms.net[fg.match$ion1, ], fg.match[, c("mz.ppm", "ion1")],
    fg.mass.diff[fg.match$ion2, ]
  )[mass_diff > 0]

  cn.net <- cn.net %>%
    dplyr::mutate(TRACE_pattern = paste0("C", C_count, "N", N_count))

  cn.net.list <- split(cn.net, cn.net$from)
  prefilt <- vapply(cn.net.list, function(x.cn) 0 %in% x.cn$N_count, FUN.VALUE = logical(1))
  cn.net.list <- cn.net.list[prefilt]

  MSdev::message_with_time("Find CN label pattern...")
  cn.net.list.hit <- BiocParallel::bplapply(
    names(cn.net.list),
    function(x) {
      x.cn <- cn.net.list[[x]]
      possible.c.count <- unique(x.cn$C_count)
      possible.n.count <- unique(x.cn$N_count)
      c.max <- x.cn$from.mz[1] / 14

      possible.c.count <- TRACE_LowC_cutoff %>%
        dplyr::filter(mass_min < x.cn$from.mz[1], mass_max > x.cn$from.mz[1]) %>%
        dplyr::pull(c.count) %>%
        intersect(possible.c.count, .)

      possible.n.count <- possible.n.count[possible.n.count < c.max]
      cn.comb <- expand.grid(C = possible.c.count, N = possible.n.count, p.cor = NA)

      if (!nrow(cn.comb)) {
        return(NULL)
      }
      #message(x)
      for (i.cn in seq_len(nrow(cn.comb))) {
        this.c <- cn.comb$C[i.cn]
        this.n <- cn.comb$N[i.cn]

        all.form <- c(
          paste0("C0N", this.n),
          paste0("C", this.c, "N0"),
          paste0("C", this.c, "N", this.n)
        )
        all.form <- setdiff(all.form, "C0N0")
        if (!all(all.form %in% x.cn$TRACE_pattern)) {
          next
        }

        to.id <- x.cn$to[match(all.form, x.cn$TRACE_pattern)]
        m.detected <- xcms.val[c(x.cn$from[1], to.id), xcms.TRACE.sample$sampleNames]
        colnames(m.detected) <- xcms.TRACE.sample$sample.type
        rownames(m.detected) <- c("C0N0", all.form)

        mean.c0n0 <- mean(
          m.detected[rownames(m.detected) == "C0N0", colnames(m.detected) == "S12C14N"]
        )
        m.detected <- m.detected / mean.c0n0
        m.ideal <- get_ideal_CN_ratio(this.c, this.n, ratio.adjust = ratio.adjust) %>% t()
        m.ideal <- m.ideal[rownames(m.detected), colnames(m.detected)]

        cn.comb$p.cor[i.cn] <- cor(as.vector(m.detected), as.vector(m.ideal))
      }

      p.cor.max <- suppressWarnings(max(cn.comb$p.cor, na.rm = TRUE))
      if (is.infinite(p.cor.max) || p.cor.max < 0) {
        return(NULL)
      }

      cn.comb <- cn.comb %>% dplyr::slice_max(p.cor, with_ties = FALSE)
      all.form <- c(
        paste0("C0N", cn.comb$N),
        paste0("C", cn.comb$C, "N0"),
        paste0("C", cn.comb$C, "N", cn.comb$N)
      )
      all.form <- setdiff(all.form, "C0N0")

      x.cn <- x.cn[match(all.form, x.cn$TRACE_pattern), ]
      x.cn$TRACE_cor <- p.cor.max
      x.cn$TRACE_formula <- paste0("C", cn.comb$C, "N", cn.comb$N)
      x.cn
    },
    BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)
  )

  names(cn.net.list.hit) <- names(cn.net.list)
  cn.net.list.hit <- cn.net.list.hit[!vapply(cn.net.list.hit, is.null, FUN.VALUE = logical(1))]
  if (!length(cn.net.list.hit)) {
    object <- .trace_set_temp(object, pol, "cn.net", cn.net)
    object <- .trace_set_temp(object, pol, "ad.net", ad.net)
    object <- .trace_set_temp(object, pol, "is.net", is.net)
    object <- .trace_set_temp(object, pol, "fg.net", fg.net)
    object <- .trace_set_temp(object, pol, "cn.net.hit", cn.net[0, ])
    object <- .trace_set_temp(object, pol, "CNfinder", cn.net)
    object <- .trace_set_temp(object, pol, "xcms.fdf", xcms.fdf)
    object <- .trace_set_temp(object, pol, "rt.tol", rt.tol)
    object <- .trace_set_temp(object, pol, "ppm", ppm)
    return(object)
  }

  cn.net.hit <- data.table::rbindlist(cn.net.list.hit) %>%
    dplyr::filter(TRACE_cor >= TRACE_cor_cutoff)
  cn.temp <- data.table::rbindlist(cn.net.list.hit)
  cn.finder <- cn.net %>%
    dplyr::mutate(TRACE_cor = cn.temp$TRACE_cor[match(ion1, cn.temp$ion1)])

  object <- .trace_set_temp(object, pol, "cn.net", cn.net)
  object <- .trace_set_temp(object, pol, "ad.net", ad.net)
  object <- .trace_set_temp(object, pol, "is.net", is.net)
  object <- .trace_set_temp(object, pol, "fg.net", fg.net)
  object <- .trace_set_temp(object, pol, "cn.net.hit", cn.net.hit)
  object <- .trace_set_temp(object, pol, "CNfinder", cn.finder)
  object <- .trace_set_temp(object, pol, "xcms.fdf", xcms.fdf)
  object <- .trace_set_temp(object, pol, "rt.tol", rt.tol)
  object <- .trace_set_temp(object, pol, "ppm", ppm)
  object
}

#' Apply dynamic m/z and RT filtering for TRACE network
#'
#' @param object MSdev object.
#' @param i.pol Polarity index. `0` for negative and `1` for positive.
#'
#' @returns Updated MSdev object with dynamic filter outputs in
#'   `object@advancedAna$TRACE_temp[[pol]]`.
#' @export
TRACE_dynamic_filter <- function(object, i.pol) {
  object <- .trace_init_temp(object)
  pol <- .trace_get_pol(i.pol)

  cn.net <- .trace_get_temp(object, pol, "cn.net")
  cn.net.hit <- .trace_get_temp(object, pol, "cn.net.hit")
  ad.net <- .trace_get_temp(object, pol, "ad.net")
  is.net <- .trace_get_temp(object, pol, "is.net")
  fg.net <- .trace_get_temp(object, pol, "fg.net")
  rt.tol <- .trace_get_temp(object, pol, "rt.tol")
  ppm <- .trace_get_temp(object, pol, "ppm")

  if (is.null(cn.net) || is.null(cn.net.hit) || !nrow(cn.net)) {
    return(object)
  }

  cn.net.eval <- cn.net %>%
    dplyr::mutate(cn.hit = ion1 %in% cn.net.hit$ion1) %>%
    dplyr::arrange(cn.hit)
  data.table::setDT(cn.net.eval)

  ppm.fit <- distinct_norm_from_random_backgroud(
    cn.net.eval[(cn.hit), mz.ppm],
    cn.net.eval[!(cn.hit), mz.ppm]
  )
  ppm.dyn <- ppm.fit$sd * stats::qnorm(0.999)

  rt.fit <- distinct_norm_from_random_backgroud(
    cn.net.eval[(cn.hit), rt.diff],
    cn.net.eval[!(cn.hit), rt.diff]
  )
  rt.tol.dyn <- rt.fit$sd * stats::qnorm(0.99999)

  object <- .trace_set_temp(object, pol, "mz.dyn", ppm.fit)
  object <- .trace_set_temp(object, pol, "rt.dyn", rt.fit)
  object <- .trace_set_temp(object, pol, "ppm.dyn", ppm.dyn)
  object <- .trace_set_temp(object, pol, "rt.tol.dyn", rt.tol.dyn)

  cols <- c("TRUE" = "red", "FALSE" = "#888888")

  p <- ggplot2::ggplot() +
    ggrastr::rasterise(
      ggplot2::geom_point(
        data = cn.net.eval,
        ggplot2::aes(x = mz.ppm, y = rt.diff, col = cn.hit),
        pch = 16, alpha = 0.2, size = 0.02
      ),
      dpi = 300
    ) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::labs(x = "mz error (ppm)", y = "rt shift (s)") +
    ggplot2::theme_bw(base_size = 6) +
    ggplot2::theme(legend.position = "none")

  p.r <- ggplot2::ggplot(cn.net.eval) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = rt.diff, x = after_stat(density), fill = cn.hit),
      position = "dodge",
      bins = 20,
      col = "white"
    ) +
    ggplot2::stat_ecdf(ggplot2::aes(y = rt.diff, col = cn.hit), linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = rt.tol.dyn * c(-1, 1), col = "red", linewidth = 0.5, lty = "dashed") +
    ggplot2::annotate(
      geom = "text", x = 0.5, y = rt.tol / 2, label = str_digit(rt.tol.dyn),
      size = 2, col = "red", check_overlap = TRUE
    ) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = c(0, 1)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_classic(base_size = 6) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.5, 0.9),
      axis.ticks = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 5, face = "bold"),
      legend.text = ggplot2::element_text(size = 4),
      legend.background = ggplot2::element_blank(),
      legend.key.size = grid::unit(0.1, "inch"),
      legend.key.spacing = grid::unit(0.02, "inch"),
      legend.title.position = "top"
    )

  p.u <- ggplot2::ggplot(cn.net.eval) +
    ggplot2::geom_histogram(
      ggplot2::aes(x = mz.ppm, y = after_stat(density), fill = cn.hit),
      position = "dodge",
      bins = 20,
      col = "white"
    ) +
    ggplot2::stat_ecdf(ggplot2::aes(x = mz.ppm, col = cn.hit), linewidth = 0.5, show.legend = FALSE) +
    ggplot2::geom_vline(xintercept = ppm.dyn * c(-1, 1), col = "red", linewidth = 0.5, lty = "dashed") +
    ggplot2::annotate(
      geom = "text", x = -ppm / 2, y = 0.8, label = str_digit(ppm.dyn),
      size = 2, col = "red", check_overlap = TRUE
    ) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = c(0, 1)) +
    ggplot2::labs(x = NULL, y = NULL, fill = "CN labeled") +
    ggplot2::theme_classic(base_size = 6) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )

  p.ur <- ggplot2::ggplot(cn.net.eval) +
    ggplot2::geom_bar(ggplot2::aes(y = 0, fill = cn.hit), position = "stack") +
    ggplot2::annotate(
      geom = "text", x = 0, y = 0, size = 2,
      label = num2percent(sum(cn.net.eval$cn.hit) / nrow(cn.net.eval))
    ) +
    ggplot2::scale_fill_manual(values = cols, guide = ggplot2::guide_legend(ncol = 1)) +
    ggplot2::labs(fill = "CN labeled") +
    ggplot2::coord_polar() +
    ggplot2::theme_void(base_size = 6) +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 5, face = "bold"),
      legend.text = ggplot2::element_text(size = 5),
      legend.key.size = grid::unit(0.1, "inch"),
      legend.key.spacing = grid::unit(0.02, "inch"),
      legend.position = "none",
      legend.title.position = "top"
    )

  p.all <- p.u + p.ur + p + p.r +
    patchwork::plot_layout(heights = c(0.2, 0.8), widths = c(0.8, 0.2)) +
    patchwork::plot_annotation(title = paste0(MSdev::get_MSdev_instrument(object), " ", pol)) &
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 6),
      plot.tag.position = "topleft",
      plot.margin = ggplot2::margin(t = 1, r = 1, b = 1, l = 1)
    )

  fo <- paste0(object@projectInfo$MSdevFile, ".TRACE.error.pdf")
  fo <- paste0(
    "C:\\Users\\91879\\OneDrive\\Documents\\YLF_Lab\\Project\\2025.10.10.TRACE\\result/dynamic error/",
    basename(fo)
  )
  MSdev::export_graph2pdf(p.all, file_path = fo, width = 3, height = 3, append = i.pol)

  cn.net.filter <- cn.net.hit[abs(mz.ppm) > ppm.dyn | abs(rt.diff) > rt.tol.dyn]
  cn.net.hit <- cn.net.hit[!from %in% cn.net.filter$from]
  ad.net <- ad.net[abs(mz.ppm) < ppm.dyn & abs(rt.diff) < rt.tol.dyn]
  is.net <- is.net[abs(mz.ppm) < ppm.dyn & abs(rt.diff) < rt.tol.dyn]
  fg.net <- fg.net[abs(mz.ppm) < ppm.dyn & abs(rt.diff) < rt.tol.dyn]

  object <- .trace_set_temp(object, pol, "cn.net.eval", cn.net.eval)
  object <- .trace_set_temp(object, pol, "cn.net.hit", cn.net.hit)
  object <- .trace_set_temp(object, pol, "ad.net", ad.net)
  object <- .trace_set_temp(object, pol, "is.net", is.net)
  object <- .trace_set_temp(object, pol, "fg.net", fg.net)
  object
}

#' Build TRACE seed network assignment
#'
#' Partitions the CN seed graph by MS-derived edge connectivity, then applies
#' PAVE-style global assignment within each component: all simple paths between
#' node pairs are scored using accumulated m/z and RT error; paths whose net
#' chemical formula change matches the MS mass-difference table are blocked.
#' Nodes are clustered into mutually compatible subnetworks using connection
#' scores.
#'
#' @param object MSdev object.
#' @param i.pol Polarity index. `0` for negative and `1` for positive.
#' @param max_path_length Maximum hops for simple-path enumeration (default `6`).
#' @param connection_cutoff Minimum pairwise connection score to link nodes in the
#'   same assignment group (default `0.5`).
#'
#' @returns Updated MSdev object with `cn.seed.ig`, `cn.seed.assign`, and
#'   `xcms.net.candidate` in `object@advancedAna$TRACE_temp[[pol]]`.
#' @export
TRACE_network_assignment <- function(
    object,
    i.pol,
    max_path_length = 6L,
    connection_cutoff = 0.5) {
  object <- .trace_init_temp(object)
  pol <- .trace_get_pol(i.pol)

  cn.net.hit <- .trace_get_temp(object, pol, "cn.net.hit")
  ad.net <- .trace_get_temp(object, pol, "ad.net")
  is.net <- .trace_get_temp(object, pol, "is.net")
  fg.net <- .trace_get_temp(object, pol, "fg.net")
  xcms.fdf <- .trace_get_temp(object, pol, "xcms.fdf")

  if (is.null(cn.net.hit) || is.null(ad.net) || is.null(is.net) || is.null(fg.net) || is.null(xcms.fdf)) {
    return(object)
  }

  if (!nrow(cn.net.hit)) {
    object <- .trace_set_temp(object, pol, "cn.seed.ig", NULL)
    object <- .trace_set_temp(object, pol, "cn.seed.assign", NULL)
    return(object)
  }

  cn.seed <- as.character(unique(cn.net.hit$from))
  cn.seed.formula <- cn.net.hit[, .SD[1], by = from]
  cn.seed.formula <- setNames(cn.seed.formula$TRACE_formula, cn.seed.formula$from)

  xcms.net.candidate <- .trace_build_candidate_net(cn.net.hit, ad.net, is.net, fg.net)
  object <- .trace_set_temp(object, pol, "xcms.net.candidate", xcms.net.candidate)

  cn.seed.net <- xcms.net.candidate[
    as.character(from) %in% cn.seed & as.character(to) %in% cn.seed
  ]
  cn.seed.net <- .trace_retype_seed_edges(
    data.table::as.data.table(cn.seed.net),
    cn.seed.formula
  )

  cn.seed.node <- data.frame(name = cn.seed, stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      color = dplyr::case_when(name %in% cn.net.hit$from ~ "#E64B35", TRUE ~ "#97C2FC"),
      TRACE_formula = cn.seed.formula[name],
      mz = xcms.fdf$mzmed[as.numeric(name)],
      rt = xcms.fdf$rtmed[as.numeric(name)],
      label = paste0(TRACE_formula, "\n", name)
    )

  cn.seed.ig.full <- igraph::graph_from_data_frame(cn.seed.net, vertices = cn.seed.node)
  cn.seed.assign <- .trace_run_global_assignment(
    cn.seed.ig.full,
    cn.seed,
    i.pol = i.pol,
    object = object,
    pol = pol,
    max_path_length = max_path_length,
    connection_cutoff = connection_cutoff
  )

  cn.seed.node <- cn.seed.node %>%
    dplyr::left_join(
      cn.seed.assign %>% dplyr::select(name, assign.group, assign.seed, conn.component),
      by = "name"
    )

  cn.seed.ig <- igraph::graph_from_data_frame(cn.seed.net, vertices = cn.seed.node)

  object <- .trace_set_temp(object, pol, "cn.seed.ig", cn.seed.ig)
  object <- .trace_set_temp(object, pol, "cn.seed.ig.full", cn.seed.ig.full)
  object <- .trace_set_temp(object, pol, "cn.seed.assign", cn.seed.assign)
  object
}

#' Annotate TRACE seed network candidates and finalize TRACE result
#'
#' @param object MSdev object.
#' @param i.pol Polarity index. `0` for negative and `1` for positive.
#' @param cpdb Path to compound table xlsx used for candidate assignment.
#'
#' @returns Updated MSdev object with `object@advancedAna$TRACE[[pol]]`.
#' @export
TRACE_annotate <- function(
    object,
    i.pol,
    cpdb = "d:/data/2025.12.26.PAVE2/trace.cp.db.xlsx") {
  object <- .trace_init_temp(object)
  pol <- .trace_get_pol(i.pol)

  cn.seed.ig <- .trace_get_temp(object, pol, "cn.seed.ig")
  if (is.null(cn.seed.ig)) {
    object@advancedAna$TRACE[[pol]] <- data.frame()
    return(object)
  }

  xcms.xcms <- object@xcmsData[[paste0(pol, "MS1")]]
  if (is.null(xcms.xcms)) {
    object@advancedAna$TRACE[[pol]] <- data.frame()
    return(object)
  }

  cpdb <- openxlsx::read.xlsx(cpdb)
  adducts <- MSCC::adduct.table %>%
    dplyr::filter((sign(Charge) + 1) / 2 == unique(ProtGenerics::polarity(xcms.xcms)), Multi == 1, abs(Charge) == 1)

  cp.adduct <- MSCC::chemform_adduct(cpdb$formula, adducts$Adduct, value = "all") %>%
    dplyr::mutate(compound_id = cpdb$compound_id[id], rt = cpdb$rt[id]) %>%
    dplyr::filter(findInterval(chemform.adduct.mz, MSdev::mzrange(xcms.xcms)) == 1)

  cn.seed.vdata <- MSdev::vdata(cn.seed.ig)
  cn.seed.assign <- .trace_get_temp(object, pol, "cn.seed.assign")
  if (!is.null(cn.seed.assign)) {
    cn.seed.vdata <- cn.seed.vdata %>%
      dplyr::left_join(
        cn.seed.assign %>% dplyr::select(name, assign.group, assign.seed, conn.component),
        by = "name"
      )
  }
  matched.df <- MSdev::match_mz_foverlaps(
    mz1 = cn.seed.vdata$mz,
    mz2 = cp.adduct$chemform.adduct.mz,
    ppm = 10
  )

  cn.seed.vdata$candidate.id <- lapply(seq_len(nrow(cn.seed.vdata)), function(i) {
    idx <- matched.df$ion2[matched.df$ion1 == i]
    cp.adduct$compound_id[as.numeric(idx)]
  })
  cn.seed.vdata$candidate.formula <- lapply(seq_len(nrow(cn.seed.vdata)), function(i) {
    idx <- matched.df$ion2[matched.df$ion1 == i]
    cp.adduct$chemform[as.numeric(idx)]
  })
  cn.seed.vdata$candidate.adduct <- lapply(seq_len(nrow(cn.seed.vdata)), function(i) {
    idx <- matched.df$ion2[matched.df$ion1 == i]
    cp.adduct$adduct[as.numeric(idx)]
  })
  cn.seed.vdata$candidate.rt <- lapply(seq_len(nrow(cn.seed.vdata)), function(i) {
    idx <- matched.df$ion2[matched.df$ion1 == i]
    cp.adduct$rt[as.numeric(idx)]
  })

  for (x in seq_len(nrow(cn.seed.vdata))) {
    x.cn <- cn.seed.vdata$TRACE_formula[x]
    x.candi.formula <- cn.seed.vdata$candidate.formula[[x]]
    x.candi.cn <- extract_formula_CN(x.candi.formula)
    id.cn.match <- x.candi.cn %in% x.cn

    cn.seed.vdata$candidate.id[[x]] <- cn.seed.vdata$candidate.id[[x]][id.cn.match]
    cn.seed.vdata$candidate.formula[[x]] <- cn.seed.vdata$candidate.formula[[x]][id.cn.match]
    cn.seed.vdata$candidate.adduct[[x]] <- cn.seed.vdata$candidate.adduct[[x]][id.cn.match]
    cn.seed.vdata$candidate.rt[[x]] <- cn.seed.vdata$candidate.rt[[x]][id.cn.match]
  }

  node.annos <- lapply(seq_len(nrow(cn.seed.vdata)), function(x) {
    x.name <- cn.seed.vdata$name[x]
    x.candi.id <- cn.seed.vdata$candidate.id[[x]]
    x.candi.formula <- cn.seed.vdata$candidate.formula[[x]]
    x.candi.adduct <- cn.seed.vdata$candidate.adduct[[x]]
    x.candi.rt <- cn.seed.vdata$candidate.rt[[x]]
    x.candi.rt[is.na(x.candi.rt)] <- Inf
    x.candi.rtd <- x.candi.rt - cn.seed.vdata$rt[x]

    x.from <- MSdev::edata(cn.seed.ig) %>%
      dplyr::filter(type != "isotope", from == x.name) %>%
      dplyr::select(type, eid, adduct = adduct.from, seed = to, fragment, element)

    x.to <- MSdev::edata(cn.seed.ig) %>%
      dplyr::filter(type != "fragment", to == x.name) %>%
      dplyr::select(type, eid, adduct = adduct.to, seed = from, fragment, element)

    anno <- dplyr::bind_rows(x.from, x.to)
    ad.score <- ifelse(x.candi.adduct %in% anno$adduct, 1, 0)
    rt.score <- 1 - abs(x.candi.rtd) / 1000
    rt.score[rt.score < 0] <- 0
    score <- (ad.score + rt.score) / 2

    if (any(score> 0) ) {
      idx <- which.max(score)
      return(data.frame(
        name = x.name,
        type = "metabolite",
        seed = x.name,
        score = score[idx],
        compound.id = x.candi.id[idx],
        compound.formula = x.candi.formula[idx],
        compound.adduct = x.candi.adduct[idx],
        compound.rt = x.candi.rt[idx]
      ))
    }

    if (any(c("fragment", "isotope") %in% anno$type)) {
      x.anno <- anno %>%
        dplyr::filter(type %in% c("fragment", "isotope")) %>%
        dplyr::slice_head(n = 1)
      x.seed <- x.anno$seed
      if (!is.null(cn.seed.assign) && x.name %in% cn.seed.assign$name) {
        x.seed <- cn.seed.assign$assign.seed[match(x.name, cn.seed.assign$name)]
      }
      return(data.frame(
        name = x.name,
        type = x.anno$type,
        seed = x.seed,
        score = 0,
        compound.id = NA,
        compound.formula = NA,
        compound.adduct = NA,
        compound.rt = NA
      ))
    }

    if (length(score) > 0) {
      return(data.frame(
        name = x.name,
        type = "metabolite",
        seed = x.name,
        score = 0,
        compound.id = x.candi.id[1],
        compound.formula = x.candi.formula[1],
        compound.adduct = x.candi.adduct[1],
        compound.rt = x.candi.rt[1]
      ))
    }

    data.frame(
      name = x.name,
      type = "unknown",
      seed = NA,
      score = 0,
      compound.id = NA,
      compound.formula = NA,
      compound.adduct = NA,
      compound.rt = NA
    )
  })

  node.anno <- data.table::rbindlist(node.annos)
  cn.seed.vdata.anned <- cn.seed.vdata %>%
    dplyr::left_join(node.anno, by = "name")

  adducts.score <- adducts %>%
    dplyr::ungroup() %>%
    dplyr::arrange(abs(Mass)) %>%
    dplyr::mutate(score = (dplyr::row_number()), score = 1 - score / max(score)) %>%
    dplyr::pull(score, name = Adduct)

  cn.seed.vdata.known <- cn.seed.vdata.anned %>%
    dplyr::filter(type %in% "metabolite") %>%
    dplyr::mutate(
      score.ad = adducts.score[compound.adduct],
      score.rt = 1 - abs(rt - compound.rt) / 1e5,
      score.rt = ifelse(score.rt < 0, 0, score.rt),
      score = score.rt + score.ad
    ) %>%
    dplyr::group_by(compound.id) %>%
    dplyr::arrange(dplyr::desc(score)) %>%
    dplyr::mutate(
      temp = seq_len(dplyr::n()),
      type = dplyr::case_when(temp == 1 ~ type, TRUE ~ "adduct"),
      seed = dplyr::first(name)
    )

  cn.seed.vdata.fi <- cn.seed.vdata.anned %>%
    dplyr::filter(type %in% c("fragment", "isotope"))

  cn.seed.vdata.unknown <- cn.seed.vdata.anned %>%
    dplyr::filter(type %in% "unknown") %>%
    dplyr::mutate(type = ifelse(name %in% cn.seed.vdata.fi$seed, "metabolite", type), seed = name)

  if (!is.null(cn.seed.assign) && "assign.group" %in% names(cn.seed.vdata.unknown)) {
    cn.seed.vdata.unknown <- cn.seed.vdata.unknown %>%
      dplyr::group_by(assign.group) %>%
      dplyr::arrange(type) %>%
      dplyr::mutate(
        temp = seq_len(dplyr::n()),
        type = dplyr::case_when(temp == 1L ~ "metabolite", TRUE ~ "adduct"),
        seed = dplyr::first(assign.seed)
      ) %>%
      dplyr::ungroup()
  } else {
    uk.m <- MSdev::igraph_filter_vertex(cn.seed.ig, dplyr::pull(cn.seed.vdata.unknown, name)) %>%
      MSdev::get_igraph_membership()
    cn.seed.vdata.unknown <- cn.seed.vdata.unknown %>%
      dplyr::mutate(membership = uk.m[name]) %>%
      dplyr::group_by(membership) %>%
      dplyr::arrange(type) %>%
      dplyr::mutate(
        temp = seq_len(dplyr::n()),
        type = dplyr::case_when(temp == 1L ~ "metabolite", TRUE ~ "adduct"),
        seed = dplyr::first(name)
      ) %>%
      dplyr::ungroup()
  }

  cn.seed.vdata3 <- rbind(cn.seed.vdata.known, cn.seed.vdata.unknown, cn.seed.vdata.fi)
  if (!is.null(cn.seed.assign)) {
    cn.seed.vdata3 <- cn.seed.vdata3 %>%
      dplyr::left_join(
        cn.seed.assign %>% dplyr::select(name, assign.group, assign.seed, conn.component),
        by = "name"
      )
  }
  cn.seed.vdata3 <- cn.seed.vdata3 %>%
    dplyr::select(
      feature_id = name,
      TRACE_formula,
      mz, rt, type, seed,
      dplyr::any_of(c("assign.group", "assign.seed", "conn.component")),
      compound_id = compound.id, compound.formula,
      compound.adduct, compound.rt
    )

  trace.res <- dplyr::left_join(
    cn.seed.vdata3,
    cpdb[, c("compound_id", "name", "kegg_id")],
    by = "compound_id"
  )
  object@advancedAna$TRACE[[pol]] <- trace.res
  object
}

#' Run the refactored TRACE workflow
#'
#' @param object MSdev object.
#' @param rt.tol RT tolerance in seconds.
#' @param ppm m/z tolerance in ppm.
#' @param cpdb Path to compound table xlsx used for candidate assignment.
#' @param eval_top Proportion in `(0, 1]` passed to
#'   `TRACE_CN_labelling_ratio_adjust()` for top-`TRACE_cor` ratio evaluation.
#' @param ratio.plot Logical. If `TRUE`, draw the CN labelling-ratio plot in
#'   `TRACE_CN_labelling_ratio_adjust()`.
#' @param ratio.reconstruct Logical. If `TRUE`, re-run `TRACE_get_CN_net()`
#'   with evaluated group ratios after CN net construction.
#'
#' @returns Updated MSdev object with final outputs in `object@advancedAna$TRACE`
#'   and intermediate outputs in `object@advancedAna$TRACE_temp`.
#' @export
TRACE_workflow <- function(
    object,
    rt.tol = 10,
    ppm = 5,
    cpdb = "d:/data/2025.12.26.PAVE2/trace.cp.db.xlsx",
    eval_top = 0.2,
    ratio.plot = FALSE,
    ratio.reconstruct = TRUE) {
  object <- .trace_init_temp(object)
  for (i.pol in 0:1) {
    object <- TRACE_get_CN_net(object, i.pol = i.pol, rt.tol = rt.tol, ppm = ppm)
  }
  object <- TRACE_CN_labelling_ratio_adjust(
    object,
    eval_top = eval_top,
    plot = ratio.plot,
    reconstruct = ratio.reconstruct
  )
  for (i.pol in 0:1) {
    object <- TRACE_dynamic_filter(object, i.pol = i.pol)
    object <- TRACE_network_assignment(object, i.pol = i.pol)
    object <- TRACE_annotate(object, i.pol = i.pol, cpdb = cpdb)
  }
  object
}

#' Get CN labelling ratio for TRACE seeds
#'
#' @param object MSdev object.
#' @param eval_top Proportion in `(0, 1]` used to select top `TRACE_cor`
#'   entries for ratio-adjust evaluation.
#' @param plot Logical. If `TRUE`, draw labelling-ratio jitter/crossbar plot
#'   for the top `eval_top` `TRACE_cor` entries.
#'
#' @returns A data.frame with columns `TRACE_seed`, `TRACE_formula`,
#'   `TRACE_cor`, `S12C14N`, `S12C15N`, `S13C14N`, `S13C15N`.
#' @export
get_TRACE_CN_labelling_ratio <- function(object, eval_top = 0.2, plot = FALSE) {
  if (!is.numeric(eval_top) || length(eval_top) != 1 || is.na(eval_top) || eval_top <= 0 || eval_top > 1) {
    stop("eval_top must be a numeric scalar in (0, 1].")
  }
  if (!is.logical(plot) || length(plot) != 1 || is.na(plot)) {
    stop("plot must be a single TRUE/FALSE value.")
  }

  sample.types <- c("S12C14N", "S12C15N", "S13C14N", "S13C15N")
  out.list <- list()

  for (i.pol in 0:1) {
    pol <- .trace_get_pol(i.pol)
    cn.hit <- object@advancedAna$TRACE_temp[[pol]]$cn.net.hit
    xcms.xcms <- object@xcmsData[[paste0(pol, "MS1")]]

    if (is.null(cn.hit) || !nrow(cn.hit) || is.null(xcms.xcms)) {
      next
    }

    cn.hit <- data.table::as.data.table(cn.hit)
    cn.seed <- cn.hit[, .SD[1], by = from]
    xcms.val <- xcms::featureValues(xcms.xcms, missing = 0, value = "maxo")
    xcms.pda <- Biobase::pData(xcms.xcms)
    xcms.trace.sample <- xcms.pda %>%
      dplyr::filter(sample.type %in% sample.types)

    if (!nrow(xcms.trace.sample)) {
      next
    }

    .mean_feature_in_group <- function(fid, sample.type) {
      if (is.na(fid)) {
        return(NA_real_)
      }
      fid <- suppressWarnings(as.numeric(fid))
      if (!is.finite(fid) || fid < 1 || fid > nrow(xcms.val)) {
        return(NA_real_)
      }
      this.sample <- xcms.trace.sample$sampleNames[xcms.trace.sample$sample.type == sample.type]
      if (!length(this.sample)) {
        return(NA_real_)
      }
      x <- xcms.val[fid, this.sample, drop = TRUE]
      x.mean <- mean(x, na.rm = TRUE)
      if (!is.finite(x.mean)) {
        return(NA_real_)
      }
      x.mean
    }

    out.list[[pol]] <- lapply(seq_len(nrow(cn.seed)), function(i.seed) {
      this.seed <- cn.seed[i.seed, ]
      this.hit <- cn.hit[from == this.seed$from]
      this.formula <- as.character(this.seed$TRACE_formula)

      # TRACE_formula is "CXNY"; we need Y for C0NY and X for CXN0.
      this.cn <- regmatches(this.formula, regexec("^C([0-9]+)N([0-9]+)$", this.formula))[[1]]
      this.c <- if (length(this.cn) == 3) as.numeric(this.cn[2]) else NA_real_
      this.n <- if (length(this.cn) == 3) as.numeric(this.cn[3]) else NA_real_

      fid.c0n0 <- this.seed$from
      fid.c0ny <- this.hit[C_count == 0 & N_count == this.n, to][1]
      fid.cxn0 <- this.hit[C_count == this.c & N_count == 0, to][1]
      fid.cxny <- this.hit[C_count == this.c & N_count == this.n, to][1]

      # For CXN0 formulas (N == 0), C0NY == C0N0 and CXNY == CXN0.
      if (is.finite(this.n) && this.n == 0) {
        fid.c0ny <- fid.c0n0
        fid.cxny <- fid.cxn0
      }

      intensity <- c(
        S12C14N = .mean_feature_in_group(fid.c0n0, "S12C14N"),
        S12C15N = .mean_feature_in_group(fid.c0ny, "S12C15N"),
        S13C14N = .mean_feature_in_group(fid.cxn0, "S13C14N"),
        S13C15N = .mean_feature_in_group(fid.cxny, "S13C15N")
      )

      baseline <- intensity[["S12C14N"]]
      baseline <- mean(intensity)
      ratio <- if (is.finite(baseline) && baseline > 0) {
        intensity / baseline
      } else {
        rep(NA_real_, length(intensity))
      }

      data.frame(
        TRACE_seed = this.seed$from,
        TRACE_formula = this.seed$TRACE_formula,
        TRACE_cor = this.seed$TRACE_cor,
        S12C14N = ratio[["S12C14N"]],
        S12C15N = ratio[["S12C15N"]],
        S13C14N = ratio[["S13C14N"]],
        S13C15N = ratio[["S13C15N"]],
        row.names = NULL
      )
    })
  }

  if (!length(out.list)) {
    return(data.frame(
      TRACE_seed = numeric(),
      TRACE_formula = character(),
      TRACE_cor = numeric(),
      S12C14N = numeric(),
      S12C15N = numeric(),
      S13C14N = numeric(),
      S13C15N = numeric()
    ))
  }

  ratio.df <- data.table::rbindlist(unlist(out.list, recursive = FALSE), fill = TRUE) %>%
    as.data.frame()

  eval.df <- ratio.df %>%
    dplyr::filter(!is.na(TRACE_cor)) %>%
    dplyr::arrange(dplyr::desc(TRACE_cor))

  if (nrow(eval.df) > 0) {
    n.keep <- max(1, ceiling(nrow(eval.df) * eval_top))
    eval.df <- utils::head(eval.df, n.keep)
    med <- apply(eval.df[, sample.types, drop = FALSE], 2, stats::median, na.rm = TRUE)
    med[!is.finite(med)] <- NA_real_
  } else {
    med <- setNames(rep(NA_real_, length(sample.types)), sample.types)
  }

  xx <- formatC(eval_top * 100, format = "f", digits = 2)
  msg.adjust <- paste(
    formatC(med[["S12C14N"]], format = "f", digits = 4),
    formatC(med[["S12C15N"]], format = "f", digits = 4),
    formatC(med[["S13C14N"]], format = "f", digits = 4),
    formatC(med[["S13C15N"]], format = "f", digits = 4),
    sep = ","
  )
  message(
    "Evaluate the adjusted labelling ratio with top ", xx,
    "% TRACE_cor, run TRACE_get_CN_net(ratio.adjust = c(",
    msg.adjust,
    ")) to get the best CN label."
  )

  if (plot) {
    plot.df <- eval.df %>%
      tidyr::pivot_longer(dplyr::all_of(sample.types), names_to = "name", values_to = "value")
    p <- ggplot2::ggplot(plot.df, ggplot2::aes(x = name, y = value, col = TRACE_cor)) +
      ggplot2::geom_jitter(alpha = 0.5) +
      ggplot2::stat_summary(
        fun = "median",
        fun.min = "median",
        fun.max = "median",
        geom = "crossbar",
        width = 0.5,
        color = "black",
        size = 0.5
      ) +
      ggplot2::geom_hline(yintercept = 1) +
      ggplot2::scale_color_gradient(low = "yellow", high = "red")
    print(p)
  }

  attr(ratio.df, "cn.group.ratio") <- med
  ratio.df
}

#' Adjust CN labelling ratio and optionally reconstruct CN hits
#'
#' This function calculates the labelling ratio of each CN hit across the
#' 4 TRACE sample groups, stores the evaluated group-level ratio in
#' `object@advancedAna$TRACE_temp$cn.group.ratio` and per-seed ratios in
#' `object@advancedAna$TRACE_temp$cn.ratio.df`, and optionally reconstructs
#' the CN hit network by re-running `TRACE_get_CN_net()` with the new
#' `ratio.adjust`.
#'
#' @param object MSdev object.
#' @param eval_top Proportion in `(0, 1]` used to select top `TRACE_cor`
#'   entries for ratio-adjust evaluation.
#' @param plot Logical. If `TRUE`, draw labelling-ratio jitter/crossbar plot
#'   for the top `eval_top` `TRACE_cor` entries.
#' @param reconstruct Logical. If `TRUE`, re-run `TRACE_get_CN_net()` using
#'   the evaluated 4-group labelling ratio as `ratio.adjust`.
#'
#' @returns Updated MSdev object.
#' @export
TRACE_CN_labelling_ratio_adjust <- function(object, eval_top = 0.2, plot = FALSE, reconstruct = TRUE) {
  if (!is.logical(reconstruct) || length(reconstruct) != 1 || is.na(reconstruct)) {
    stop("reconstruct must be a single TRUE/FALSE value.")
  }

  ratio.df <- get_TRACE_CN_labelling_ratio(object, eval_top = eval_top, plot = plot)
  cn.group.ratio <- attr(ratio.df, "cn.group.ratio", exact = TRUE)
  sample.types <- c("S12C14N", "S12C15N", "S13C14N", "S13C15N")

  if (is.null(object@advancedAna$TRACE_temp)) {
    object@advancedAna$TRACE_temp <- list()
  }
  if (is.null(object@advancedAna$TRACE_temp$cn.group.ratio)) {
    object@advancedAna$TRACE_temp$cn.group.ratio <- list()
  }
  object@advancedAna$TRACE_temp$cn.group.ratio[["all"]] <- cn.group.ratio
  object@advancedAna$TRACE_temp$cn.ratio.df <- ratio.df

  if (reconstruct) {
    ratio.adjust <- cn.group.ratio[sample.types]
    ratio.adjust <- as.numeric(ratio.adjust)
    ratio.adjust[!is.finite(ratio.adjust)] <- 1
    if (length(ratio.adjust) != 4) {
      ratio.adjust <- rep(1, 4)
    }

    for (i.pol in 0:1) {
      pol <- .trace_get_pol(i.pol)
      rt.tol <- .trace_get_temp(object, pol, "rt.tol")
      ppm <- .trace_get_temp(object, pol, "ppm")
      if (is.null(rt.tol)) rt.tol <- 10
      if (is.null(ppm)) ppm <- 5

      object <- TRACE_get_CN_net(
        object,
        i.pol = i.pol,
        rt.tol = rt.tol,
        ppm = ppm,
        ratio.adjust = ratio.adjust
      )
    }
  }

  object
}

#' Export TRACE results to xlsx
#'
#' Writes TRACE results to an Excel workbook, including one sheet per polarity
#' with **all xcms features** and their TRACE annotations (if present).
#' Additional sheets include the raw TRACE annotation tables and (when
#' available) CN labelling-ratio outputs from `object@advancedAna$TRACE_temp`.
#'
#' @param object MSdev object with TRACE results in `object@advancedAna$TRACE`.
#' @param file Output xlsx path. Default is
#'   `file.path(object@projectInfo$projectDir, "Statistic", "TRACE.xlsx")`
#'   when `projectDir` is set; otherwise
#'   `paste0(object@projectInfo$MSdevFile, ".TRACE.xlsx")`.
#'
#' @returns Invisibly, the output file path.
#' @export
TRACE_export <- function(object, file = NULL) {
  trace <- object@advancedAna$TRACE
  if (is.null(trace)) {
    stop("No TRACE results found in object@advancedAna$TRACE. Run TRACE_workflow() first.")
  }

  if (is.null(file)) {
    if (!is.null(object@projectInfo$projectDir)) {
      file <- file.path(object@projectInfo$projectDir, "Statistic", "TRACE.xlsx")
    } else if (!is.null(object@projectInfo$MSdevFile)) {
      file <- paste0(object@projectInfo$MSdevFile, ".TRACE.xlsx")
    } else {
      stop("file must be provided when object@projectInfo$projectDir and MSdevFile are NULL.")
    }
  }

  df.list <- list()
  for (pol in c("Negative", "Positive")) {
    x <- trace[[pol]]

    # All features (from xcms) + left-join TRACE annotations by feature_id.
    xcms.xcms <- object@xcmsData[[paste0(pol, "MS1")]]
    if (!is.null(xcms.xcms)) {
      fdf <- xcms::featureDefinitions(xcms.xcms)
      all.df <- as.data.frame(fdf)
      all.df$feature_id <- as.character(seq_len(nrow(all.df)))

      if (!is.null(x) && nrow(x)) {
        x.df <- as.data.frame(x)
        x.df$feature_id <- as.character(x.df$feature_id)
        all.df <- dplyr::left_join(all.df, x.df, by = "feature_id")
      }

      df.list[[paste0("AllFeatures_", pol)]] <- all.df
    }

    # Raw TRACE annotation table (only features in seed network)
    if (!is.null(x) && nrow(x)) {
      df.list[[paste0("TRACE_", pol)]] <- as.data.frame(x)
    }
  }

  ratio.df <- object@advancedAna$TRACE_temp$cn.ratio.df
  if (!is.null(ratio.df)) {
    df.list[["CN_labelling_ratio"]] <- as.data.frame(ratio.df)
  }

  cn.group.ratio <- object@advancedAna$TRACE_temp$cn.group.ratio[["all"]]
  if (!is.null(cn.group.ratio) && length(cn.group.ratio)) {
    df.list[["CN_group_ratio"]] <- data.frame(
      sample.type = names(cn.group.ratio),
      ratio = as.numeric(cn.group.ratio),
      row.names = NULL,
      check.names = FALSE
    )
  }

  if (!length(df.list)) {
    stop("No TRACE data available to export.")
  }

  out.dir <- dirname(file)
  if (nzchar(out.dir) && out.dir != ".") {
    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
  }

  MSdev::xlsx.write.list(df.list, file = file)
  MSdev::message_with_time(
    "Exported TRACE results to: ",
    normalizePath(file, winslash = "/", mustWork = FALSE)
  )
  invisible(file)
}

#' Plot TRACE CN labeling fraction per metabolite and sample
#'
#' For each TRACE metabolite seed, counts all CN-labeled features in the seed
#' network (C0N0, C0Ny, CxN0, CxNy, and other matched CxNy patterns) and
#' calculates the labeling fraction in each sample:
#' \itemize{
#'   \item `S12C14N`: C0N0 / total
#'   \item `S12C15N`: C0Ny / total
#'   \item `S13C14N`: CxN0 / total
#'   \item `S13C15N`: CxNy / total
#' }
#'
#' @param object MSdev object with TRACE results in `object@advancedAna$TRACE`
#'   and CN hits in `object@advancedAna$TRACE_temp`.
#' @param i.pol Polarity index. `0` for negative and `1` for positive.
#' @param TRACE_cor_cutoff Minimum `TRACE_cor` used to select CN seeds when
#'   TRACE metabolite annotations are unavailable.
#' @param return.data Logical. If `TRUE`, return a list with `plot` and `data`.
#'
#' @returns A patchwork object combining the box/point plot and density plot,
#'   or a list with `plot`, `p1`, `p2`, and `data` when `return.data = TRUE`.
#' @export
plot_TRACE_labeling_fraction <- function(
    object,
    i.pol = 0,
    TRACE_cor_cutoff = 0.7,
    return.data = FALSE) {
  sample.types <- c("S12C14N", "S12C15N", "S13C14N", "S13C15N")
  pol <- .trace_get_pol(i.pol)

  cn.hit <- object@advancedAna$TRACE_temp[[pol]]$cn.net.hit
  trace.res <- object@advancedAna$TRACE[[pol]]
  xcms.xcms <- object@xcmsData[[paste0(pol, "MS1")]]

  if (is.null(cn.hit) || !nrow(cn.hit) || is.null(xcms.xcms)) {
    stop("TRACE CN hits or xcms data not found for polarity ", pol, ".")
  }

  cn.hit <- data.table::as.data.table(cn.hit)
  xcms.val <- xcms::featureValues(xcms.xcms, missing = 0, value = "maxo")
  xcms.pda <- Biobase::pData(xcms.xcms)
  xcms.trace.sample <- xcms.pda %>%
    dplyr::filter(sample.type %in% sample.types)

  if (!nrow(xcms.trace.sample)) {
    stop("No TRACE sample groups found in xcms sample metadata.")
  }

  seed.df <- NULL
  if (!is.null(trace.res) && nrow(trace.res)) {
    seed.df <- trace.res %>%
      dplyr::filter(type == "metabolite") %>%
      dplyr::transmute(
        seed = as.character(seed),
        metabolite = dplyr::coalesce(as.character(name), as.character(compound_id), seed)
      ) %>%
      dplyr::distinct(seed, .keep_all = TRUE)
  }

  if (is.null(seed.df) || !nrow(seed.df)) {
    cn.seed <- cn.hit[, .SD[1], by = from]
    cn.seed <- cn.seed[TRACE_cor >= TRACE_cor_cutoff]
    seed.df <- cn.seed %>%
      dplyr::transmute(
        seed = as.character(from),
        metabolite = paste0(TRACE_formula, "_", from)
      )
  }

  .feature_intensity <- function(fid, sample.name) {
    if (is.na(fid)) {
      return(NA_real_)
    }
    fid <- suppressWarnings(as.numeric(fid))
    if (!is.finite(fid) || fid < 1 || fid > nrow(xcms.val)) {
      return(NA_real_)
    }
    x <- xcms.val[fid, sample.name, drop = TRUE]
    x <- x[is.finite(x)]
    if (!length(x)) {
      return(NA_real_)
    }
    mean(x)
  }

  frac.list <- lapply(seq_len(nrow(seed.df)), function(i.seed) {
    this.seed <- seed.df$seed[i.seed]
    this.met <- seed.df$metabolite[i.seed]
    this.hit <- cn.hit[as.character(from) == this.seed]
    if (!nrow(this.hit)) {
      return(NULL)
    }

    this.formula <- as.character(this.hit$TRACE_formula[1])
    this.cn <- regmatches(this.formula, regexec("^C([0-9]+)N([0-9]+)$", this.formula))[[1]]
    this.c <- if (length(this.cn) == 3) as.numeric(this.cn[2]) else NA_real_
    this.n <- if (length(this.cn) == 3) as.numeric(this.cn[3]) else NA_real_

    fid.c0n0 <- this.seed
    fid.c0ny <- this.hit[C_count == 0 & N_count == this.n, to][1]
    fid.cxn0 <- this.hit[C_count == this.c & N_count == 0, to][1]
    fid.cxny <- this.hit[C_count == this.c & N_count == this.n, to][1]

    if (is.finite(this.n) && this.n == 0) {
      fid.c0ny <- fid.c0n0
      fid.cxny <- fid.cxn0
    }

    all.fids <- unique(c(fid.c0n0, this.hit$to))
    all.fids <- all.fids[!is.na(all.fids)]

    lapply(seq_len(nrow(xcms.trace.sample)), function(i.sample) {
      sample.name <- xcms.trace.sample$sampleNames[i.sample]
      sample.type <- xcms.trace.sample$sample.type[i.sample]

      total <- sum(vapply(all.fids, .feature_intensity, numeric(1), sample.name = sample.name), na.rm = TRUE)
      if (!is.finite(total) || total <= 0) {
        return(NULL)
      }

      fid.num <- switch(
        as.character(sample.type),
        S12C14N = fid.c0n0,
        S12C15N = fid.c0ny,
        S13C14N = fid.cxn0,
        S13C15N = fid.cxny,
        NA
      )
      numerator <- .feature_intensity(fid.num, sample.name)
      if (!is.finite(numerator)) {
        return(NULL)
      }

      data.frame(
        metabolite = this.met,
        seed = this.seed,
        TRACE_formula = this.formula,
        sample = sample.name,
        sample.type = sample.type,
        labeling_fraction = numerator / total,
        total_intensity = total,
        numerator_intensity = numerator,
        n_cn_features = length(all.fids),
        row.names = NULL
      )
    })
  })

  plot.df <- data.table::rbindlist(
    unlist(frac.list, recursive = FALSE),
    fill = TRUE
  ) %>%
    as.data.frame()

  if (!nrow(plot.df)) {
    stop("No labeling-fraction values could be calculated.")
  }

  sample.labels <- c(
    S12C14N = "S12C14N\nC0N0/total",
    S12C15N = "S12C15N\nC0Ny/total",
    S13C14N = "S13C14N\nCxN0/total",
    S13C15N = "S13C15N\nCxNy/total"
  )
  cn.pattern.colors <- stats::setNames(
    c("#C9352D", "#FF7F0E", "#1F77B4", "#9467BD", "#CCCCCC"),
    c("C0N0", "CxN0", "C0Ny", "CxNy", "noise")
  )
  sample.pattern <- c(
    S12C14N = "C0N0",
    S13C14N = "CxN0",
    S12C15N = "C0Ny",
    S13C15N = "CxNy"
  )
  sample.colors <- cn.pattern.colors[sample.pattern]
  names(sample.colors) <- names(sample.pattern)

  plot.df$sample.type <- factor(plot.df$sample.type, levels = sample.types)
  plot.df$sample.label <- sample.labels[as.character(plot.df$sample.type)]
  plot.df$sample.label <- factor(plot.df$sample.label, levels = sample.labels)

  p1 <- ggplot2::ggplot(
    plot.df,
    ggplot2::aes(
      x = sample.label,
      y = labeling_fraction,
      color = sample.type,
      fill = sample.type
    )
  ) +
    ggplot2::geom_boxplot(
      outlier.shape = NA,
      alpha = 0.55
    ) +
    #ggplot2::geom_point(
    #  position = ggplot2::position_jitter(width = 0.12, height = 0),
    #  alpha = 0.55,
    #  size = 1.2
    #) +
    ggplot2::scale_color_manual(
      values = sample.colors,
      labels = sample.pattern,
      name = "CN pattern"
    ) +
    ggplot2::scale_fill_manual(
      values = sample.colors,
      labels = sample.pattern,
      name = "CN pattern"
    ) +
    ggplot2::labs(
      x = "Sample group (numerator/total)",
      y = "Labeling fraction",
      title = "Labeling fraction by metabolite"
    ) +
    ggplot2::theme_bw(base_size = 6) +
    ggplot2::theme(
     # axis.text.x = ggplot2::element_text(size = 9),
      legend.position = "none"
    )

  p2 <- ggplot2::ggplot(
    plot.df,
    ggplot2::aes(
      x = labeling_fraction,
      color = sample.type,
      fill = sample.type
    )
  ) +
    ggplot2::geom_density(alpha = 0.25) +
    ggplot2::scale_color_manual(
      values = sample.colors,
      labels = sample.pattern,
      name = "CN pattern"
    ) +
    ggplot2::scale_fill_manual(
      values = sample.colors,
      labels = sample.pattern,
      name = "CN pattern"
    ) +
    ggplot2::labs(
      x = "Labeling fraction",
      y = "Density",
      title = "Labeling fraction distribution"
    ) +
    ggplot2::theme_bw(base_size = 6) +
    ggplot2::theme(legend.position = c(0.1,0.7),
                   legend.key.size  = unit(0.1,"inch"))

  p <- p1 / p2 +
    patchwork::plot_layout(widths = c(1.2, 1)) +
    patchwork::plot_annotation(tag_levels = "A",
      #title = paste0("TRACE labeling fraction (", pol, ")")
    )

  if (return.data) {
    return(list(plot = p, p1 = p1, p2 = p2, data = plot.df))
  }
  p
}

#' Plot feature RSD distribution by sample group
#'
#' For each sample group (`group` in xcms sample metadata, also referred to as
#' `sample.group` in MSdev sample sheets), computes the relative standard
#' deviation (RSD = SD / mean) of feature intensities across biological
#' replicates within that group, then plots the RSD distribution across features.
#'
#' @param object MSdev object with processed xcms data in `object@xcmsData`.
#' @param i.pol Polarity index. `0` for negative and `1` for positive.
#' @param sample.groups Character vector of sample groups to include. Default
#'   uses all `group` levels with at least `min.replicates` samples, optionally
#'   excluding `Blank` and `QC`.
#' @param exclude.blank Logical. If `TRUE`, exclude `Blank` and `QC` groups when
#'   `sample.groups` is `NULL`.
#' @param min.replicates Minimum number of samples required in a group to
#'   calculate RSD.
#' @param binwidth Histogram bin width for the RSD histogram.
#' @param rsd.max Upper x-axis limit for RSD plots.
#' @param return.data Logical. If `TRUE`, return a list with `plot` and `data`.
#'
#' @returns A patchwork object combining the histogram and ECDF plots, or a
#'   list with `plot`, `p1`, `p2`, and `data` when `return.data = TRUE`.
#' @export
plot_TRACE_RSD <- function(
    object,
    i.pol = 0,
    sample.groups = NULL,
    exclude.blank = TRUE,
    min.replicates = 2,
    binwidth = 0.05,
    rsd.max = 1,
    return.data = FALSE) {
  pol <- .trace_get_pol(i.pol)
  xcms.xcms <- object@xcmsData[[paste0(pol, "MS1")]]

  if (is.null(xcms.xcms)) {
    stop("xcms data not found for polarity ", pol, ".")
  }

  xcms.val <- xcms::featureValues(xcms.xcms, missing = 0, value = "maxo")
  pda <- Biobase::pData(xcms.xcms)

  if (is.null(pda$group)) {
    stop("Sample metadata column 'group' (sample.group) not found in xcms pData.")
  }

  if (is.null(sample.groups)) {
    sample.groups <- unique(as.character(pda$group))
    sample.groups <- sample.groups[nzchar(sample.groups)]
    if (exclude.blank) {
      sample.groups <- sample.groups[!sample.groups %in% c("Blank", "QC")]
    }
  }

  .feature_rsd <- function(x) {
    x <- x[is.finite(x) & x > 0]
    if (length(x) < min.replicates) {
      return(NA_real_)
    }
    stats::sd(x) / mean(x)
  }

  rsd.list <- lapply(sample.groups, function(this.group) {
    samps <- rownames(pda)[as.character(pda$group) == this.group]
    if (length(samps) < min.replicates) {
      return(NULL)
    }

    rsd <- apply(xcms.val[, samps, drop = FALSE], 1, .feature_rsd)
    data.frame(
      feature_id = seq_along(rsd),
      sample.group = this.group,
      rsd = rsd,
      n_replicates = length(samps),
      row.names = NULL
    )
  })

  plot.df <- as.data.frame(data.table::rbindlist(rsd.list, fill = TRUE))
  plot.df <- plot.df[is.finite(plot.df$rsd), , drop = FALSE]

  if (!nrow(plot.df)) {
    stop("No RSD values could be calculated for the selected sample groups.")
  }

  trace.types <- c("S12C14N", "S12C15N", "S13C14N", "S13C15N")
  group.levels <- unique(c(
    intersect(trace.types, sample.groups),
    setdiff(unique(plot.df$sample.group), trace.types)
  ))
  plot.df$sample.group <- factor(plot.df$sample.group, levels = group.levels)

  cn.pattern.colors <- stats::setNames(
    c("#C9352D", "#FF7F0E", "#1F77B4", "#9467BD", "#CCCCCC"),
    c("C0N0", "CxN0", "C0Ny", "CxNy", "noise")
  )
  sample.pattern <- c(
    S12C14N = "C0N0",
    S12C15N = "C0Ny",
    S13C14N = "CxN0",
    S13C15N = "CxNy",
    Blank = "noise",
    QC = "noise"
  )
  group.colors <- cn.pattern.colors[sample.pattern[as.character(group.levels)]]
  names(group.colors) <- as.character(group.levels)
  missing.groups <- is.na(group.colors)
  if (any(missing.groups)) {
    extra.colors <- grDevices::colorRampPalette(c("#4DAF4A", "#984EA3", "#A65628"))(sum(missing.groups))
    group.colors[missing.groups] <- extra.colors
    names(group.colors) <- as.character(group.levels)
  }

  p1 <- ggplot2::ggplot(
    plot.df,
    ggplot2::aes(
      x = rsd,
      fill = sample.group
    )
  ) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(count)),
      position = ggplot2::position_dodge(width = binwidth),
      binwidth = binwidth,
      color = "black",
      alpha = 0.75
    ) +
    ggplot2::scale_fill_manual(values = group.colors, name = "Sample group") +
    ggplot2::scale_x_continuous(limits = c(0, rsd.max), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(
      x = "RSD",
      y = "Feature count",
      title = "Feature RSD distribution by sample group"
    ) +
    ggplot2::theme_bw(base_size = 6) +
    ggplot2::theme(legend.position = "none")

  p2 <- ggplot2::ggplot(
    plot.df,
    ggplot2::aes(
      x = rsd,
      color = sample.group
    )
  ) +
    ggplot2::stat_ecdf(linewidth = 0.6) +
    ggplot2::scale_color_manual(values = group.colors, name = "Sample group") +
    ggplot2::scale_x_continuous(limits = c(0, rsd.max), expand = c(0, 0)) +
    ggplot2::labs(
      x = "RSD",
      y = "Cumulative probability",
      title = "Feature RSD cumulative distribution"
    ) +
    ggplot2::theme_bw(base_size = 6) +
    ggplot2::theme(legend.position = c(0.15, 0.7),
                   legend.key.size = grid::unit(0.1, "inch"))

  p <- p1 / p2 +
    patchwork::plot_layout(widths = c(1.2, 1)) +
    patchwork::plot_annotation(
      tag_levels = "A",
      title = paste0("TRACE feature RSD (", pol, ")")
    )

  if (return.data) {
    return(list(plot = p, p1 = p1, p2 = p2, data = plot.df))
  }
  p
}
