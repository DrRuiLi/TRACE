

# Mon Jun  1 20:01:43 2026 ------------------------------
{


  {
    trace.cor <- 0.5
    trace.demo <- TRACE_get_CN_net(trace.demo,0,TRACE_cor_cutoff = trace.cor,
                                   ratio.adjust = c(1.0842,0.7641,1.0566,1.0707))
    trace.demo <- TRACE_get_CN_net(trace.demo,1,TRACE_cor_cutoff = trace.cor,
                                   ratio.adjust = c(1.0842,0.7641,1.0566,1.0707))
  }

  a <- get_TRACE_CN_labelling_ratio(trace.demo,eval_top = 1,plot = T)
  a <- trace.demo@advancedAna$TRACE_temp$cn.ratio.df

  #trace.demo <- TRACE_CN_labelling_ratio_adjust(trace.demo,eval_top = 0.3,plot = T,reconstruct = T)

  ratio.adj <- c(1.0688,0.6341,1.0480,1.1131)
  names(ratio.adj) <- names(a)[4:7]

  df <- a%>%
    pivot_longer(4:7)%>%
    dplyr::slice_max(TRACE_cor,prop = 1 )%>%
    dplyr::mutate(
      regeion = cut(TRACE_cor, breaks = seq(0, 1, 0.05)),
      ratio.bench = ratio.adj[name],
      ratio.error = value - ratio.bench
    )



  p <- ggplot(df,aes(x = name , y = value, col = TRACE_cor))+
    geom_jitter(
      alpha = 0.2
    )+
    stat_summary(
      fun = "median",
      fun.min = "median",
      fun.max = "median",
      geom = "crossbar",
      width = 0.5,
      color = "black",
      size = 0.5
    ) +
    geom_hline(yintercept = 1)+
    scale_color_gradient(low = "yellow",high = "red")+
    labs(x = NULL, y = "Ratio")+
    theme_bw()

  open_plot_win(p,5,3)


  p <- ggplot(df,aes(x = TRACE_cor , y = ratio.error,col = TRACE_cor))+
    geom_point(alpha = 0.3)+
    scale_color_gradient(low = "yellow",high = "red")+
    labs(y = "Ratio shift")+
    theme_bw()

  open_plot_win(p)

  library(ggridges)

  p <- ggplot(df,aes(x = ratio.error , y = regeion,fill = regeion))+
    geom_density_ridges(   )+
    scale_fill_manual(values =  MSdev:::colramp()(seq(0.1,1,0.1)))+
    labs(x = "Ratio shift", y = NULL, fill = expression( rho ~ "range"))+
    theme_bw(base_size =  6 )+
    theme(legend.key.size = unit(0.1,"inch") )

  open_plot_win(p,2.5,3)



  object <- TRACE_network_assignment(object ,i.pol = 0)
  object <- TRACE_annotate(object,i.pol = 0)

  object <- TRACE_network_assignment(object ,i.pol = 1)
  object <- TRACE_annotate(object,i.pol = 1)


}
# Thu Jun  4 16:02:11 2026 ------------------------------
{


  obj <- MSdev_load("d:/data/2025.12.26.PAVE2/PAVE_With_Params/OE480_120k_ppm10_sn10.rdata")
  obj <- MSdev:::.update_MSdev_object(obj)
  obj <- TRACE_workflow(obj)

}

# Fri Jun  5 14:10:00 2026 ------------------------------
# Compare PAVE-matlab and TRACE annotation results
{
  library(dplyr)
  library(ggplot2)
  library(MSdev)

  pave.file <- "c:/Users/91879/OneDrive/Documents/YLF_Lab/Project/2025.10.10.PAVE/PAVE-Matlab/example/pks_features.xlsx"
  trace.file <- "d:/temp/TRACE.xlsx"
  trace.sheet <- "AllFeatures_Negative"

  # PAVE xlsx has 800+ columns (scoremat_*); read only cols 1:17 for comparison.
  pave <- openxlsx::read.xlsx(
    pave.file,
    cols = 1:17
  ) %>%
    dplyr::select(feature_id, feature, C_num, N_num, score) %>%
    dplyr::mutate(
      feature_id_num = as.integer(gsub("[^0-9]+", "", feature_id)),
      pave_type = dplyr::case_when(
        is.na(feature) | feature == "Background" ~ "blank/noise",
        TRUE ~ as.character(feature)
      ),
      pave_formula = paste0("C", C_num, "N", N_num)
    )

  # TRACE AllFeatures sheet: read only feature_id, TRACE_formula, type, seed.
  trace <- openxlsx::read.xlsx(
    trace.file,
    sheet = trace.sheet,
    cols = c(1, 26, 29, 30)
  ) %>%
    dplyr::mutate(
      feature_id_num = as.integer(feature_id),
      trace_type = dplyr::case_when(
        is.na(type) | type == "" ~ "blank/noise",
        TRUE ~ as.character(type)
      ),
      trace_formula = as.character(TRACE_formula)
    )

  ratio.df <- tryCatch(
    openxlsx::read.xlsx(
      trace.file,
      sheet = "CN_labelling_ratio",
      cols = c(1, 2, 3)
    ),
    error = function(e) NULL
  )

  cmp <- dplyr::inner_join(
    pave %>%
      dplyr::select(
        feature_id_num, pave_type, pave_formula, score
      ),
    trace %>%
      dplyr::select(
        feature_id_num, trace_type, trace_formula, seed
      ),
    by = "feature_id_num"
  )
  cmp$cn_match <- !is.na(cmp$trace_formula) &
    cmp$trace_formula != "" &
    cmp$trace_type != "blank/noise" &
    cmp$pave_type != "blank/noise" &
    cmp$pave_formula == cmp$trace_formula

  if (!is.null(ratio.df) && all(c("TRACE_seed", "TRACE_cor") %in% names(ratio.df))) {
    trace.cor.map <- stats::setNames(ratio.df$TRACE_cor, as.character(ratio.df$TRACE_seed))
    cmp$trace_cor <- as.numeric(trace.cor.map[as.character(cmp$seed)])
  } else {
    cmp$trace_cor <- NA_real_
  }
  cmp$trace_cor <- ifelse(is.finite(cmp$trace_cor), cmp$trace_cor, cmp$score)

  type.plot.df <- cmp %>%
    dplyr::group_by(trace_type, pave_type) %>%
    dplyr::summarize(
      n = dplyr::n(),
      trace_cor = mean(trace_cor, na.rm = TRUE),
      .groups = "drop"
    )
  type.plot.df <- type.plot.df %>%
    dplyr::mutate(n_plot = pmin(n, 2000L))

  p1 <- ggplot2::ggplot(
    type.plot.df,
    ggplot2::aes(
      x = trace_type,
      y = pave_type,
      color = trace_cor,
      size = n_plot
    )
  ) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::scale_color_gradient(low = "yellow", high = "red", name = "TRACE cor") +
    ggplot2::scale_size_area(name = "Features (capped at 2000)", max_size = 12) +
    ggplot2::labs(
      x = "TRACE annotation type",
      y = "PAVE-matlab annotation type",
      title = "TRACE vs PAVE-matlab annotation comparison"
    ) +
    ggplot2::theme_bw()

  if (interactive()) {
    open_plot_win(p1, 5, 4)
  } else {
    print(p1)
  }

  trace.cn.set <- cmp %>%
    dplyr::filter(trace_type != "blank/noise", !is.na(trace_formula), trace_formula != "") %>%
    dplyr::pull(feature_id_num) %>%
    unique()
  pave.cn.set <- cmp %>%
    dplyr::filter(pave_type != "blank/noise", is.finite(score), score > 0.7) %>%
    dplyr::pull(feature_id_num) %>%
    unique()

  venn.sets <- list(TRACE = trace.cn.set, `PAVE-matlab` = pave.cn.set)

  if (requireNamespace("VennDiagram", quietly = TRUE)) {
    if (interactive()) {
      VennDiagram::venn.diagram(
        x = venn.sets,
        filename = NULL,
        imagetype = "png",
        fill = c("steelblue", "grey80"),
        alpha = c(0.4, 0.4),
        cex = 1.2,
        cat.cex = 1.0,
        cat.pos = c(-20, 20),
        margin = 0.1
      ) %>% grid::grid.draw()
    } else {
      # In non-interactive runs, draw without opening a new device.
      g <- VennDiagram::venn.diagram(
        x = venn.sets,
        filename = NULL,
        imagetype = "png",
        fill = c("steelblue", "grey80"),
        alpha = c(0.4, 0.4),
        cex = 1.2,
        cat.cex = 1.0,
        cat.pos = c(-20, 20),
        margin = 0.1
      )
      grid::grid.newpage()
      grid::grid.draw(g)
    }
  }

  cn.eval <- cmp %>%
    dplyr::filter(
      trace_type != "blank/noise",
      pave_type != "blank/noise",
      !is.na(trace_formula),
      trace_formula != ""
    )

  message(
    "Aligned ", nrow(cmp), " features by feature_id; CN formula agreement (non-blank/noise) ",
    formatC(100 * mean(cn.eval$cn_match, na.rm = TRUE), format = "f", digits = 1), "%"
  )
}


