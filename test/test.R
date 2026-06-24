

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

# Fri Jun  6 15:00:00 2026 ------------------------------
# CN net compare: exclusive CN-assigned features (TRACE vs PAVE-matlab)
{
  library(dplyr)
  library(ggplot2)
  library(MSdev)

  pave.file <- "c:/Users/91879/OneDrive/Documents/YLF_Lab/Project/2025.10.10.PAVE/PAVE-Matlab/example/pks_features.xlsx"
  trace.file <- "d:/temp/TRACE.xlsx"
  trace.sheet <- "AllFeatures_Negative"

  pave <- openxlsx::read.xlsx(pave.file, cols = 1:17) %>%
    dplyr::select(feature_id, feature, C_num, N_num, score) %>%
    dplyr::mutate(
      feature_id_num = as.integer(gsub("[^0-9]+", "", feature_id)),
      pave_cn = !is.na(feature) &
        feature != "Background" &
        is.finite(score) &
        score > 0.7
    )

  trace <- openxlsx::read.xlsx(
    trace.file,
    sheet = trace.sheet,
    cols = c(1, 26, 29, 30)
  ) %>%
    dplyr::mutate(
      feature_id_num = as.integer(feature_id),
      trace_cn = !is.na(type) &
        type != "" &
        !is.na(TRACE_formula) &
        TRACE_formula != ""
    )

  cmp <- dplyr::inner_join(
    pave %>% dplyr::select(feature_id_num, pave_cn),
    trace %>% dplyr::select(feature_id_num, trace_cn),
    by = "feature_id_num"
  )

  trace.cn.ids <- cmp$feature_id_num[cmp$trace_cn]
  pave.cn.ids <- cmp$feature_id_num[cmp$pave_cn]
  trace.only.ids <- setdiff(trace.cn.ids, pave.cn.ids)
  pave.only.ids <- setdiff(pave.cn.ids, trace.cn.ids)
  both.cn.ids <- intersect(trace.cn.ids, pave.cn.ids)

  cn.dist.df <- data.frame(
    category = c(
      "TRACE CN only",
      "PAVE-matlab CN only",
      "Both CN assigned"
    ),
    n = c(
      length(trace.only.ids),
      length(pave.only.ids),
      length(both.cn.ids)
    ),
    stringsAsFactors = FALSE
  )
  cn.dist.df$ratio <- cn.dist.df$n / sum(cn.dist.df$n)

  cn.exclusive.df <- data.frame(
    source = c("TRACE", "PAVE-matlab"),
    exclusive_n = c(length(trace.only.ids), length(pave.only.ids)),
    total_cn = c(length(trace.cn.ids), length(pave.cn.ids)),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      exclusive_ratio = exclusive_n / total_cn,
      label = paste0(
        formatC(100 * exclusive_ratio, format = "f", digits = 1),
        "% (n=", exclusive_n, ")"
      )
    )

  p.cn.exclusive <- ggplot2::ggplot(
    cn.exclusive.df,
    ggplot2::aes(x = source, y = exclusive_ratio, fill = source)
  ) +
    ggplot2::geom_col(width = 0.6, alpha = 0.85) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      vjust = -0.3,
      size = 3.5
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1.05),
      expand = c(0, 0)
    ) +
    ggplot2::scale_fill_manual(values = c("TRACE" = "steelblue", "PAVE-matlab" = "grey60")) +
    ggplot2::labs(
      x = NULL,
      y = "Ratio of CN-assigned features not shared",
      title = "Exclusive CN assignment: TRACE vs PAVE-matlab",
      subtitle = paste0(
        "TRACE CN: ", length(trace.cn.ids),
        "; PAVE CN: ", length(pave.cn.ids),
        "; overlap: ", length(both.cn.ids)
      ),
      fill = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")

  if (interactive()) {
    open_plot_win(p.cn.exclusive, 5, 4)
  } else {
    print(p.cn.exclusive)
  }

  message(
    "CN exclusive ratio — TRACE only: ",
    formatC(100 * cn.exclusive.df$exclusive_ratio[1], format = "f", digits = 1),
    "% (", cn.exclusive.df$exclusive_n[1], "/", cn.exclusive.df$total_cn[1], "); ",
    "PAVE only: ",
    formatC(100 * cn.exclusive.df$exclusive_ratio[2], format = "f", digits = 1),
    "% (", cn.exclusive.df$exclusive_n[2], "/", cn.exclusive.df$total_cn[2], ")"
  )
}

# Fri Jun  6 15:30:00 2026 ------------------------------
# Merge PAVE-matlab and TRACE results by feature_id
{
  library(dplyr)

  pave.file <- "c:/Users/91879/OneDrive/Documents/YLF_Lab/Project/2025.10.10.PAVE/PAVE-Matlab/example/pks_features.xlsx"
  trace.file <- "d:/temp/TRACE.xlsx"
  trace.sheet <- "AllFeatures_Negative"

  pave.df.raw <- openxlsx::read.xlsx(pave.file)
  pave.df <- pave.df.raw %>%
    dplyr::mutate(
      feature_id = as.integer(gsub("[^0-9]+", "", feature_id)),
      formula_pave = paste0("C", C_num, "N", N_num),
      type_pave = dplyr::case_when(
        is.na(feature) | feature == "Background" ~ "blank/noise",
        TRUE ~ as.character(feature)
      ),
      type_pave = case_when(
        type_pave %in% c("Metabolite","Isotope","Fragment","Adduct")~type_pave,
        T~"blank/noise"
      ),
      score_pave = score,
      pave_cn = formula_pave != "C0N0"
    ) %>%
    dplyr::select(
      feature_id,
      mz_pave = mz,
      rt_pave = rt,
      formula_pave,
      type_pave,
      score_pave,
      pave_cn
    )

  trace.df <- openxlsx::read.xlsx(
    trace.file,
    sheet = trace.sheet,
    cols = c(1, 26, 27, 28, 29, 30)
  ) %>%
    dplyr::mutate(
      feature_id = as.integer(feature_id),
      formula_trace = as.character(TRACE_formula),
      type_trace = dplyr::case_when(
        is.na(type) | type == "" ~ "blank/noise",
        TRUE ~ as.character(type)
      ),
      trace_cn =  !is.na(TRACE_formula),
      seed_trace = seed
    ) %>%
    dplyr::select(
      feature_id,
      mz_trace = mz,
      rt_trace = rt,
      formula_trace,
      type_trace,
      seed_trace,
      trace_cn
    )

  ratio.df <- tryCatch(
    openxlsx::read.xlsx(
      trace.file,
      cols = c(1, 2, 3)
    ),
    error = function(e) NULL
  )

  pave.trace.merged <- dplyr::full_join(pave.df, trace.df, by = "feature_id") %>%
    dplyr::mutate(
      type_pave = if_else(is.na(type_pave),"blank/noise",type_pave),
      cn_match = trace_cn &
        pave_cn &
        !is.na(formula_trace) &
        formula_trace != "" &
        formula_pave == formula_trace,
      cn_assignment = dplyr::case_when(
        trace_cn & pave_cn ~ "both",
        trace_cn & !pave_cn ~ "TRACE_only",
        !trace_cn & pave_cn ~ "PAVE_only",
        TRUE ~ "neither"
      )
    )

  if (!is.null(ratio.df) && all(c("TRACE_seed", "TRACE_cor") %in% names(ratio.df))) {
    trace.cor.map <- stats::setNames(ratio.df$TRACE_cor, as.character(ratio.df$TRACE_seed))
    pave.trace.merged$trace_cor <- as.numeric(
      trace.cor.map[as.character(pave.trace.merged$seed_trace)]
    )
  } else {
    pave.trace.merged$trace_cor <- NA_real_
  }
  pave.trace.merged$trace_cor <- ifelse(
    is.finite(pave.trace.merged$trace_cor),
    pave.trace.merged$trace_cor,
    pave.trace.merged$score_pave
  )

  pave.trace.merged <- pave.trace.merged %>%
    dplyr::select(
      feature_id,
      formula_pave,
      formula_trace,
      type_pave,
      type_trace,
      score_pave,
      trace_cor,
      pave_cn,
      trace_cn,
      cn_match,
      cn_assignment,
      mz_pave,
      rt_pave,
      mz_trace,
      rt_trace,
      seed_trace
    ) %>%
    dplyr::arrange(feature_id)


  ggplot(pave.trace.merged,aes(x = type_trace,y = type_pave))+
    geom_count()+
    scale_size(range = c(1,20),
               transform = "sqrt")+
    stat_sum(
      aes(label = after_stat(n)),
      geom = "text",
      color = "white",
      size = 3.5,
      show.legend = FALSE
    )


}


# Fri Jun 12 10:43:59 2026 Nutrients ------------------------------
{
  # TRACE
  {

    library(MSdev)
    library(tidyverse)
    library(SummarizedExperiment)
    obj <- MSdev_load("d:/data/2026.01.07.PAVE.Nutrition/MSdev_2026_01_07.Rdata")

    ## positive
    {

      xcms.pos <- obj@xcmsData$PositiveMS1
      pave.seed.pos <- obj@statData$TRACE$Positive%>%
        dplyr::filter(type== "metabolite")

      xcms.se.raw <- get_xcms_feature_se(xcms.pos,missing = 1)
      xcms.se.pos <- xcms.se.raw[,grepl("NS",xcms.se.raw$sample.type)]
      xcms.se.pos <- se_adjuset_by_weight(xcms.se.pos)
      xcms.fdf.pos <- as.data.frame(rowData(xcms.se.pos))
      xcms.fdf.pos$peakMaxo <- rowMeans(assay(xcms.se.raw[,xcms.se.raw$group == "S12C14N"]))
      xcms.mat.pos <- assay(xcms.se.pos)


      pave.label.ratio <- list()
      pb <- MSdev:::get_progress_bar(nrow(pave.seed.pos))
      for (i in 1:nrow(pave.seed.pos)) {

        pb$tick()
        i.pave.formula <- pave.seed.pos$pave_formula[i]
        i.fid <- pave.seed.pos$feature_id[i]%>%as.numeric()
        i.mz <- xcms.fdf.pos$mzmed[i.fid]

        #i.mz <- 148.0604
        #i.pave.formula <- "C5N1"
        i.cn.count <- MSCC::chemform_parse(i.pave.formula)
        i.cn.diff <- get_CN_mass_diff_table(i.cn.count[1,"C"],i.cn.count[1,"N"])
        i.cn.table <- i.cn.diff[,mz:=mass_diff + i.mz]
        i.cn.match <- MSdev:::match_mz_foverlaps(i.cn.diff$mz,xcms.fdf.pos$mzmed,ppm = 3)
        i.cn.match[,fid := ion2] [,rt := xcms.fdf.pos$rtmed[fid]]
        i.cn.table <- cbind(i.cn.match,i.cn.table[i.cn.match$ion1]  )
        i.cn.table[,label_pattern := paste0("C",C_count ,"N",N_count )][
          !is.na(fid)]
        #d <- density(i.cn.table$rt,bw = 5)
        # plot(d)
        #i.rt <- d$x[which.max(d$y)]
        i.rt <-pave.seed.pos$rt[i]
        i.cn.table <- i.cn.table[abs(rt -i.rt) < 5]%>%
          dplyr::group_by(label_pattern)%>%
          dplyr::slice_min(mz.ppm,n = 1)%>%
          dplyr::ungroup()%>%
          dplyr::arrange(N_count,C_count)

        if(nrow(i.cn.table) <= 1) next
        i.cn.exp.mat <- xcms.mat.pos[i.cn.table$fid,xcms.se.pos$sample.name]
        rownames(i.cn.exp.mat) <-i.cn.table$label_pattern
        #norm.to <- i.cn.exp.mat[,xcms.se.pos$group == "NSA"]%>%sum(na.rm = T)/10
        norm.to <- i.cn.exp.mat%>%mean(na.rm = T)
        i.cn.exp.mat <- (i.cn.exp.mat/norm.to)

        ### C/N label ratio
        if(T){
          i.cn.label.ratio <- i.cn.exp.mat%>%
            as.data.frame()%>%
            rownames_to_column("label_pattern")%>%
            pivot_longer(-"label_pattern" ,names_to = "sample")%>%
            left_join(i.cn.table,by = "label_pattern")%>%
            dplyr::group_by(sample)%>%
            dplyr::mutate(
              group = setNames(nm = xcms.se.pos$sample.name,object = xcms.se.pos$sample.type)[sample],
              int.sum = sum(value),
              c.cum = sum(C_count * value),
              c.ratio = c.cum/int.sum,
              n.cum = sum(N_count * value),
              n.ratio = n.cum/int.sum
            )%>%
            #dplyr::filter(grepl("NS",group))%>%
            #dplyr::filter(!grepl("Blank",group))%>%
            dplyr::ungroup()%>%
            dplyr::mutate(x = c.ratio + rnorm(n(),sd = 0.0001),
                          y = n.ratio + rnorm(n(),sd = 0.0001))

          res <- i.cn.label.ratio%>%
            dplyr::filter(mz == max(mz))%>%
            dplyr::mutate(fid = i.fid,
                          label_pattern = i.pave.formula,
                          mz = i.mz)

          pave.label.ratio[[i]] <- res




          if(F){

            p1 <- ggplot(i.cn.label.ratio,aes(x = x,y= y,col = group , fill = group))+
              geom_point(size = 1,alpha = 0.5)+
              stat_ellipse(geom = "polygon",alpha = 0.3)+
              labs(x = "Labeled fraction of C",
                   y = "Labeled fraction of N",
                   fill = "Treatment",col ="Treatment"
              )+
              scale_fill_manual(values = unname( ggsci:::ggsci_db$rickandmorty$schwifty),
                                labels =c("Control","+AA/uracil\n/adenine","+leucine","+threonine",
                                          "+tryptophan","+adenine","+uracil","+acetate",
                                          "12C14N","13C14N","12C15N","13C15N")  )+
              scale_color_manual(values = unname( ggsci:::ggsci_db$rickandmorty$schwifty),
                                 labels =c("Control","+AA/uracil\n/adenine","+leucine","+threonine",
                                           "+tryptophan","+adenine","+uracil","+acetate",
                                           "12C14N","13C14N","12C15N","13C15N")  )+
              theme_bw(base_size = 6)+
              theme(
                plot.tag = element_text(size = 12,face = "plain",vjust  = 0,hjust = 1 ),
                plot.tag.position = "topleft",
                plot.margin = margin(t = 2, r = 0, b = 2, l = 3 ),
                legend.key.size = unit(0.1,"inch"),
                legend.text = element_text(size = 4))

            open_plot_win(p1,2.3,2)
          }

          if(F){


            plot.data <- i.cn.label.ratio%>%
              dplyr::filter(sample == c("NS__1_1","NS__2_1"))%>%
              dplyr::group_by(sample)%>%
              dplyr::mutate(
                c.labeled = C_count * value ,
                c.unlabeled = max(C_count) * value - c.labeled) %>%
              pivot_longer(c(c.labeled,c.unlabeled),names_to = "label", values_to = "v")%>%
              dplyr::mutate(
                C_count = ifelse(label == "c.unlabeled",0,C_count)
              )%>%
              dplyr::mutate(v = v/sum(v) *max(C_count) )%>%
              dplyr::group_by(sample,C_count)%>%
              dplyr::mutate(v = sum(v))%>%
              dplyr::distinct(sample,C_count,.keep_all = T)%>%
              dplyr::arrange(C_count)%>%
              dplyr::ungroup()%>%
              dplyr::mutate(x = factor(group),
                            x = as.numeric(x))

            p.label.compose <- ggplot(plot.data)+
              geom_bar(aes(x = group , y = v, fill = C_count),
                       col = "#333333", stat = "identity")+
              scale_fill_gradient2(low = "#FFFFFF",high = "#D13B2E")+
              scale_x_discrete(labels = c( "Control\n\n","Nutrition\n\n"))+
              coord_flip()+
              labs(fill = "Count of\nLabeled C",tag = "A")+
              theme_void(base_size = 6)+
              theme(
                plot.tag = element_text(size = 12,face = "plain",vjust  = 0,hjust = 1 ),
                axis.text.y = element_text(size = 6),
                plot.margin = margin(t = 5, r = 0, b = 20, l = 0 ),
                legend.title = element_text(size = 6),
                legend.text = element_text(size = 4),
                legend.key.size = unit(0.05,"inch"))

            open_plot_win(p.label.compose,2,1.2)
            open_plot_ppt(p.label.compose,2,1.2)


            absorb <- i.cn.label.ratio%>%
              dplyr::mutate(
                c.absorb = 1 - c.ratio/max(C_count),
                n.absorb = 1 - n.ratio/max(N_count),
                x = c.absorb +  rnorm(n(),sd = 0.00001),
                y = n.absorb + rnorm(n(),sd = 0.0001)

              )%>%
              dplyr::filter(grepl("NS",group))

            p1 <- ggplot(absorb,aes(x = x,y= y,col = group , fill = group))+
              geom_point(size = 1,alpha = 0.5)+
              stat_ellipse(geom = "polygon",alpha = 0.3)+
              labs(tag = "B",x = "Absorption fraction of C",
                   y = "Absorption fraction of N",
                   fill = "Treatment",col ="Treatment"
              )+
              scale_fill_manual(values = col.map,
                                labels = labs.map )+
              scale_color_manual(values =col.map,
                                 labels =labs.map )+
              theme_bw(base_size = 6)+
              theme(
                legend.position = c(0.02,0.98),
                legend.justification = c(0,1),
                legend.background = element_rect(fill = "transparent"),
                plot.tag = element_text(size = 12,face = "plain",vjust  = 0,hjust = 1 ),
                plot.tag.position = "topleft",
                plot.margin = margin(t = 0, r = 2, b = 0, l = 0 ),
                legend.key.size = unit(0.08,"inch"),
                legend.text = element_text(size = 6))

            open_plot_win(p1,1.6,1.6)
            open_plot_ppt(p1,1.6,1.6)

            absorb.amount   <- absorb%>%
              dplyr::mutate(
                c.absorb = c.absorb * int.sum,
                n.absorb =n.absorb * int.sum,
                x = c.absorb +  rnorm(n(),sd = 0.00001),
                y = n.absorb + rnorm(n(),sd = 0.00001)
              )

            p2 <- ggplot(absorb.amount,aes(x = x,y= y,col = group , fill = group))+
              geom_point(size = 1,alpha = 0.5)+
              stat_ellipse(geom = "polygon",alpha = 0.3)+
              labs(tag = "C",x = "Relative absorbed C",
                   y = "Relative absorbed N",
                   fill = "Treatment",col ="Treatment"
              )+
              scale_fill_manual(values = unname( ggsci:::ggsci_db$npg$nrc),
                                labels =c("Control","+AA/uracil\n/adenine","+leucine","+threonine",
                                          "+tryptophan","+adenine","+uracil","+acetate",
                                          "12C14N","13C14N","12C15N","13C15N")  )+
              scale_color_manual(values = unname( ggsci:::ggsci_db$npg$nrc),
                                 labels =c("Control","+AA/uracil\n/adenine","+leucine","+threonine",
                                           "+tryptophan","+adenine","+uracil","+acetate",
                                           "12C14N","13C14N","12C15N","13C15N")  )+
              theme_bw(base_size = 6)
            #open_plot_win(p2,2,1.4)

            p <- ((p1/p2+plot_layout(guides = "collect")) )&
              theme(
                legend.position = "none",
                plot.tag = element_text(size = 12,face = "plain",vjust  = 0,hjust = 1 ),
                plot.tag.position = "topleft",
                plot.margin = margin(t = 0, r = 2, b = 0, l = 0 ),
                legend.key.size = unit(0.1,"inch"),
                legend.text = element_text(size = 4))
            open_plot_win(p, 1.5,3)

          }

        }

        #if(nrow(i.cn.exp.mat) > 4 ){
        if(F){
          n.ratio <- res %>%
            dplyr::pull(n.ratio,name = sample)
          c.ratio <- res %>%
            dplyr::pull(c.ratio,name = sample)
          int.sum <- res %>%
            dplyr::pull(int.sum,name = sample)


          #hm.mat <- log2(i.cn.exp.mat)
          hm.mat <- (i.cn.exp.mat)
          hm.mat <- apply(hm.mat,2, function(x){ x/sum(x) })
          hm.mat <- hm.mat[(rowSums(hm.mat) > 0.1),]
          hm <- Heatmap(hm.mat,
                        #col = colramp(breaks = c(min(min(hm.mat),-1),0,
                        #                         max((hm.mat))),colors = c("#259644","white","#D84704")),
                        col = MSdev:::colramp(breaks = c(0,max(hm.mat)/2,max(hm.mat))),

                        row_names_gp = gpar(fontsize = 6),
                        heatmap_legend_param = list(title =paste0("Abundance"),
                                                    grid_width  = unit(0.1, "inch"),
                                                    title_gp  =gpar(fontsize = 6),
                                                    labels_gp = gpar(fontsize = 6)),

                        #column_title = c("Control","+AA/uracil/adenine","+leucine","+threonine",
                        #                 "+tryptophan","+adenine","+uracil","+acetate",
                        #                 "12C14N","13C14N","12C15N","13C15N"
                        #),
                        column_title_side = "bottom",
                        column_title_rot = -30,
                        column_title_gp = gpar(fontsize = 6),

                        show_heatmap_legend = F,

                        cluster_rows = F,cluster_columns = F,cluster_column_slices = F,
                        column_split = xcms.se.pos$group,
                        show_column_names = F,
                        row_names_side = "left",
                        rect_gp = gpar(col = "black"),
                        top_annotation = columnAnnotation(
                          #intsum = anno_barplot(int.sum,gp = gpar(fill = "#DF3A2D"),
                          #                      axis_param = list(at = ceiling(max(int.sum)) * c(0.5,1),
                          #                                        gp = gpar(fontsize = 6),
                          #                                        side = "right")),
                          cratio = anno_barplot(c.ratio,gp = gpar(fill = "#FF7F0E"),
                                                axis_param = list(at = ceiling(max(c.ratio)) * c(0.5,1),
                                                                  gp = gpar(fontsize = 6),
                                                                  side = "right")),
                          nratio = anno_barplot(n.ratio,gp = gpar(fill = "#1F77B4"),
                                                axis_param = list(at = ceiling(max(n.ratio)) * c(0.5,1),
                                                                  gp = gpar(fontsize = 6),
                                                                  side = "right")),
                          annotation_label = c("Labeled\nFraction of C","Labeled\nFraction of N"),
                          annotation_name_gp  = gpar(fontsize = 6),
                          annotation_name_side  = "left",
                          annotation_name_rot  = 0,
                          height = unit(0.7,"inch")
                        )

          )
          hm
          open_plot_win(hm, 4,1.3 + nrow(hm.mat) * 0.12)
          #export_graph2pdf(hm,file_path = "d:/temp.pdf",append = T,width = 5,height = 3)

          if(F){
            ### i = 24
            p <- get_ggplot_from_heatmap(hm)+
              labs(tag = "G")+
              theme(
                plot.tag = element_text(size = 12,face = "plain",vjust  = 0),
                plot.tag.position = "topleft",
                plot.margin = margin(t = 0, r = 0, b = 0, l = 0 ))
            open_plot_win(p,4,4)


          }

        }
        #readline()

      }

      pave.label.ratio.df.pos <- rbindlist(pave.label.ratio)%>%
        dplyr::group_by(fid)%>%
        dplyr::mutate(
          kegg.id = setNames(pave.seed.pos$kegg_id,
                             pave.seed.pos$feature_id)[as.character(fid)],
          name = setNames(pave.seed.pos$name,
                          pave.seed.pos$feature_id)[as.character(fid)],
          peakMaxo = xcms.fdf.pos$peakMaxo[fid],
          fid = paste0("pos",fid),
          c.absorb = 1 - c.ratio/max(C_count),
          n.absorb = 1 - n.ratio/max(N_count),
          x = c.absorb +  rnorm(n(),sd = 0.00001),
          y = n.absorb + rnorm(n(),sd = 0.00001),
          c.ratio = c.ratio/max(C_count),
          n.ratio = c.ratio/max(N_count))%>%
        dplyr::filter(grepl("NS",group))
    }


    ## negative
    {

      xcms.xcms.neg <- obj@xcmsData$NegativeMS1
      pave.seed.neg <- obj@statData$TRACE$Negative%>%
        dplyr::filter(type== "metabolite")

      xcms.se.raw <- get_xcms_feature_se(xcms.xcms.neg,missing = 1)
      xcms.se.neg <- xcms.se.raw[,grepl("NS",xcms.se.raw$sample.type)]
      xcms.se.neg <- se_adjuset_by_weight(xcms.se.neg)

      xcms.fdf.neg <- as.data.frame(rowData(xcms.se.neg))
      xcms.fdf.neg$peakMaxo <- rowMeans(assay(xcms.se.raw[,xcms.se.raw$group == "S12C14N"]))
      xcms.mat.neg <- assay(xcms.se.neg)


      pave.label.ratio <- list()
      pb <- MSdev::: get_progress_bar(nrow(pave.seed.neg))
      for (i in 1:nrow(pave.seed.neg)) {

        pb$tick()
        i.pave.formula <- pave.seed.neg$pave_formula[i]
        i.fid <- pave.seed.neg$feature_id[i]%>%as.numeric()
        i.mz <- xcms.fdf.neg$mzmed[i.fid]

        #i.mz <- 664.116394
        #i.pave.formula <- "C21N7"
        i.cn.count <- MSCC::chemform_parse(i.pave.formula)
        i.cn.diff <- get_CN_mass_diff_table(i.cn.count[1,"C"],i.cn.count[1,"N"])
        i.cn.table <- i.cn.diff[,mz:=mass_diff + i.mz]
        i.cn.match <- MSdev:::match_mz_foverlaps(i.cn.diff$mz,xcms.fdf.neg$mzmed)
        i.cn.match[,fid := ion2] [,rt := xcms.fdf.neg$rtmed[fid]]
        i.cn.table <- cbind(i.cn.match,i.cn.table[i.cn.match$ion1]  )
        i.cn.table[,label_pattern := paste0("C",C_count ,"N",N_count )][
          !is.na(fid)]
        #d <- density(i.cn.table$rt,bw = 5)
        # plot(d)
        #i.rt <- d$x[which.max(d$y)]
        i.rt <-pave.seed.neg$rt[i]
        i.cn.table <- i.cn.table[abs(rt -i.rt) < 5]%>%
          dplyr::group_by(label_pattern)%>%
          dplyr::slice_min(mz.ppm,n = 1)%>%
          dplyr::ungroup()%>%
          dplyr::arrange(N_count,C_count)

        if(nrow(i.cn.table) <= 1) next
        i.cn.exp.mat <- xcms.mat.neg[i.cn.table$fid,xcms.se.neg$sample.name]
        rownames(i.cn.exp.mat) <-i.cn.table$label_pattern
        norm.to <- i.cn.exp.mat[,xcms.se.neg$group == "NSA"]%>%sum(na.rm = T)/10
        #norm.to <- i.cn.exp.mat%>%mean(na.rm = T)
        i.cn.exp.mat <- (i.cn.exp.mat/norm.to)

        ### C/N label ratio
        if(T){

          i.cn.label.ratio <- i.cn.exp.mat%>%
            as.data.frame()%>%
            rownames_to_column("label_pattern")%>%
            pivot_longer(-"label_pattern" ,names_to = "sample")%>%
            left_join(i.cn.table,by = "label_pattern")%>%
            dplyr::group_by(sample)%>%
            dplyr::mutate(
              group = setNames(nm = xcms.se.neg$sample.name,object = xcms.se.neg$sample.type)[sample],
              int.sum = sum(value),
              c.cum = sum(C_count * value),
              c.ratio = c.cum/int.sum,
              n.cum = sum(N_count * value),
              n.ratio = n.cum/int.sum
            )%>%
            #dplyr::filter(grepl("NS",group))%>%
            #dplyr::filter(!grepl("Blank",group))%>%
            dplyr::ungroup()%>%
            dplyr::mutate(x = c.ratio + rnorm(n(),sd = 0.0001),
                          y = n.ratio + rnorm(n(),sd = 0.0001))

          res <- i.cn.label.ratio%>%
            dplyr::filter(mz == max(mz))%>%
            dplyr::mutate(fid = i.fid,
                          mz = i.mz)
          pave.label.ratio[[i]] <- res




          if(F){

            p1 <- ggplot(i.cn.label.ratio,aes(x = x,y= y,col = group , fill = group))+
              geom_point(size = 1,alpha = 0.5)+
              stat_ellipse(geom = "polygon",alpha = 0.3)+
              labs(x = "Labeled fraction of C",
                   y = "Labeled fraction of N",
                   fill = "Treatment",col ="Treatment"
              )+
              scale_fill_manual(values = unname( ggsci:::ggsci_db$rickandmorty$schwifty),
                                labels =c("Control","+AA/uracil\n/adenine","+leucine","+threonine",
                                          "+tryptophan","+adenine","+uracil","+acetate",
                                          "12C14N","13C14N","12C15N","13C15N")  )+
              scale_color_manual(values = unname( ggsci:::ggsci_db$rickandmorty$schwifty),
                                 labels =c("Control","+AA/uracil\n/adenine","+leucine","+threonine",
                                           "+tryptophan","+adenine","+uracil","+acetate",
                                           "12C14N","13C14N","12C15N","13C15N")  )+
              theme_bw(base_size = 6)+
              theme(
                plot.tag = element_text(size = 12,face = "plain",vjust  = 0,hjust = 1 ),
                plot.tag.position = "topleft",
                plot.margin = margin(t = 2, r = 0, b = 2, l = 3 ),
                legend.key.size = unit(0.1,"inch"),
                legend.text = element_text(size = 4))

            open_plot_win(p1,2.3,2)
          }

          if(F){


            plot.data <- i.cn.label.ratio%>%
              dplyr::filter(sample == c("NS__1_1","NS__2_1"))%>%
              dplyr::group_by(sample)%>%
              dplyr::mutate(
                c.labeled = C_count * value ,
                c.unlabeled = max(C_count) * value - c.labeled) %>%
              pivot_longer(c(c.labeled,c.unlabeled),names_to = "label", values_to = "v")%>%
              dplyr::mutate(
                C_count = ifelse(label == "c.unlabeled",0,C_count)
              )%>%
              dplyr::mutate(v = v/sum(v) *max(C_count) )%>%
              dplyr::group_by(sample,C_count)%>%
              dplyr::mutate(v = sum(v))%>%
              dplyr::distinct(sample,C_count,.keep_all = T)%>%
              dplyr::arrange(C_count)%>%
              dplyr::ungroup()%>%
              dplyr::mutate(x = factor(group),
                            x = as.numeric(x))

            p.label.compose <- ggplot(plot.data)+
              geom_bar(aes(x = group , y = v, fill = C_count),
                       col = "#333333", stat = "identity")+
              scale_fill_gradient2(low = "#FFFFFF",high = "#D13B2E")+
              scale_x_discrete(labels = c( "Control\n\n","Nutrition\n\n"))+
              coord_flip()+
              labs(fill = "Count of\nLabeled C",tag = "A")+
              theme_void(base_size = 6)+
              theme(
                plot.tag = element_text(size = 12,face = "plain",vjust  = 0,hjust = 1 ),
                axis.text.y = element_text(size = 6),
                plot.margin = margin(t = 5, r = 0, b = 20, l = 0 ),
                legend.title = element_text(size = 6),
                legend.text = element_text(size = 4),
                legend.key.size = unit(0.05,"inch"))

            open_plot_win(p.label.compose,2,1.2)


            absorb <- i.cn.label.ratio%>%
              dplyr::mutate(
                c.absorb = 1 - c.ratio/max(C_count),
                n.absorb = 1 - n.ratio/max(N_count),
                x = c.absorb +  rnorm(n(),sd = 0.00001),
                y = n.absorb + rnorm(n(),sd = 0.00001)

              )%>%
              dplyr::filter(grepl("NS",group))

            p1 <- ggplot(absorb,aes(x = x,y= y,col = group , fill = group))+
              geom_point(size = 1,alpha = 0.5)+
              stat_ellipse(geom = "polygon",alpha = 0.3)+
              labs(tag = "B",x = "Absorption fraction of C",
                   y = "Absorption fraction of N",
                   fill = "Treatment",col ="Treatment"
              )+
              scale_fill_manual(values = unname( ggsci:::ggsci_db$npg$nrc),
                                labels =c("Control","+AA/uracil\n/adenine","+leucine","+threonine",
                                          "+tryptophan","+adenine","+uracil","+acetate",
                                          "12C14N","13C14N","12C15N","13C15N")  )+
              scale_color_manual(values = unname( ggsci:::ggsci_db$npg$nrc),
                                 labels =c("Control","+AA/uracil\n/adenine","+leucine","+threonine",
                                           "+tryptophan","+adenine","+uracil","+acetate",
                                           "12C14N","13C14N","12C15N","13C15N")  )+
              theme_bw(base_size = 6)+
              theme(
                legend.position = "none",
                plot.tag = element_text(size = 12,face = "plain",vjust  = 0,hjust = 1 ),
                plot.tag.position = "topleft",
                plot.margin = margin(t = 0, r = 2, b = 0, l = 0 ),
                legend.key.size = unit(0.1,"inch"),
                legend.text = element_text(size = 4))

            open_plot_win(p1,1.8,1.8)

            absorb.amount   <- absorb%>%
              dplyr::mutate(
                c.absorb = c.absorb * int.sum,
                n.absorb =n.absorb * int.sum,
                x = c.absorb +  rnorm(n(),sd = 0.00001),
                y = n.absorb + rnorm(n(),sd = 0.00001)
              )

            p2 <- ggplot(absorb.amount,aes(x = x,y= y,col = group , fill = group))+
              geom_point(size = 1,alpha = 0.5)+
              stat_ellipse(geom = "polygon",alpha = 0.3)+
              labs(tag = "C",x = "Relative absorbed C",
                   y = "Relative absorbed N",
                   fill = "Treatment",col ="Treatment"
              )+
              scale_fill_manual(values = unname( ggsci:::ggsci_db$npg$nrc),
                                labels =c("Control","+AA/uracil\n/adenine","+leucine","+threonine",
                                          "+tryptophan","+adenine","+uracil","+acetate",
                                          "12C14N","13C14N","12C15N","13C15N")  )+
              scale_color_manual(values = unname( ggsci:::ggsci_db$npg$nrc),
                                 labels =c("Control","+AA/uracil\n/adenine","+leucine","+threonine",
                                           "+tryptophan","+adenine","+uracil","+acetate",
                                           "12C14N","13C14N","12C15N","13C15N")  )+
              theme_bw(base_size = 6)
            #open_plot_win(p2,2,1.4)

            p <- ((p1/p2+plot_layout(guides = "collect")) )&
              theme(
                legend.position = "none",
                plot.tag = element_text(size = 12,face = "plain",vjust  = 0,hjust = 1 ),
                plot.tag.position = "topleft",
                plot.margin = margin(t = 0, r = 2, b = 0, l = 0 ),
                legend.key.size = unit(0.1,"inch"),
                legend.text = element_text(size = 4))
            open_plot_win(p, 1.5,3)

          }

        }

        #if(nrow(i.cn.exp.mat) > 4 ){
        if(F){
          n.ratio <- res %>%
            dplyr::pull(n.ratio,name = sample)
          c.ratio <- res %>%
            dplyr::pull(c.ratio,name = sample)
          int.sum <- res %>%
            dplyr::pull(int.sum,name = sample)


          #hm.mat <- log2(i.cn.exp.mat)
          hm.mat <- (i.cn.exp.mat)
          hm.mat <- apply(hm.mat,2, function(x){ x/sum(x) })
          hm.mat <- hm.mat[(rowSums(hm.mat) > 0.1),]
          hm <- Heatmap(hm.mat,
                        #col = colramp(breaks = c(min(min(hm.mat),-1),0,
                        #                         max((hm.mat))),colors = c("#259644","white","#D84704")),
                        col = colramp(breaks = c(0,max(hm.mat)/2,max(hm.mat))),

                        row_names_gp = gpar(fontsize = 6),
                        heatmap_legend_param = list(title =paste0("Abundance"),
                                                    grid_width  = unit(0.1, "inch"),
                                                    title_gp  =gpar(fontsize = 6),
                                                    labels_gp = gpar(fontsize = 6)),

                        column_title = c("Control","+AA/uracil/adenine","+leucine","+threonine",
                                         "+tryptophan","+adenine","+uracil","+acetate"
                                         #"12C14N","13C14N","12C15N","13C15N"
                        ),
                        column_title_side = "bottom",
                        column_title_rot = -30,
                        column_title_gp = gpar(fontsize = 6),


                        cluster_rows = F,cluster_columns = F,cluster_column_slices = F,
                        column_split = xcms.se.neg$group,
                        show_column_names = F,
                        row_names_side = "left",
                        rect_gp = gpar(col = "black"),
                        top_annotation = columnAnnotation(
                          cratio = anno_barplot(c.ratio,gp = gpar(fill = "#FF7F0E"),
                                                axis_param = list(at = ceiling(max(c.ratio)) * c(0.5,1),
                                                                  gp = gpar(fontsize = 6),
                                                                  side = "right")),
                          nratio = anno_barplot(n.ratio,gp = gpar(fill = "#1F77B4"),
                                                axis_param = list(at = ceiling(max(n.ratio)) * c(0.5,1),
                                                                  gp = gpar(fontsize = 6),
                                                                  side = "right")),
                          annotation_label = c("Labeled\nFraction of C","Labeled\nFraction of N"),
                          annotation_name_gp  = gpar(fontsize = 6),
                          annotation_name_side  = "left",
                          annotation_name_rot  = 0,
                          height = unit(1,"inch")
                        )

          )
          hm
          open_plot_win(hm, 4,1.3 + nrow(hm.mat) * 0.2)
          #export_graph2pdf(hm,file_path = "d:/temp.pdf",append = T,width = 5,height = 3)

          if(F){
            ### i = 24
            p <- get_ggplot_from_heatmap(hm)+
              labs(tag = "G")+
              theme(
                plot.tag = element_text(size = 12,face = "plain",vjust  = 0),
                plot.tag.position = "topleft",
                plot.margin = margin(t = 0, r = 0, b = 0, l = 0 ))
            open_plot_win(p,4,4)


          }

        }
        #readline()

      }

      pave.label.ratio.df.neg <- rbindlist(pave.label.ratio)%>%
        dplyr::group_by(fid)%>%
        dplyr::mutate(
          kegg.id = setNames(pave.seed.neg$kegg_id,
                             pave.seed.neg$feature_id)[as.character(fid)],
          name = setNames(pave.seed.neg$name,
                          pave.seed.neg$feature_id)[as.character(fid)],
          peakMaxo = xcms.fdf.neg$peakMaxo[fid],
          fid = paste0("neg",fid),
          c.absorb = 1 - c.ratio/max(C_count),
          n.absorb = 1 - n.ratio/max(N_count),
          x = c.absorb +  rnorm(n(),sd = 0.00001),
          y = n.absorb + rnorm(n(),sd = 0.00001),
          c.ratio = c.ratio/max(C_count),
          n.ratio = c.ratio/max(N_count))%>%
        dplyr::filter(grepl("NS",group))


    }


    pave.label.ratio.df <- rbind(pave.label.ratio.df.pos,
                                 pave.label.ratio.df.neg)




    ### Label ratio
    {

      library(patchwork)
      col.map <- c("#E64B35","#4DBBD5","#00A087","#3C5488",
                   "#F39B7F","#8491B4","#91D1C2","#7E6148")
      labs.map <- c("Control","+mixture","+leucine","+threonine",
                    "+tryptophan","+adenine","+uracil","+acetate")
      p1 <- ggplot(pave.label.ratio.df)+
        #geom_point(aes(x = c.ratio,y = n.ratio, col = group))+
        geom_histogram(aes(fill = group , x = c.absorb),col = "black",bins = 20,show.legend = F)+
        scale_fill_manual(values = col.map,
                          labels = labs.map)+
        labs(x = "Absorption fraction of C", y = "Count",
             fill = "Treatment", tag = "C")+
        theme_bw(base_size = 6)
      #p1

      p2 <- ggplot(pave.label.ratio.df)+
        geom_histogram(aes(fill = group , x = n.absorb),col = "black",bins = 20,show.legend = F)+
        scale_fill_manual(values = col.map,
                          labels = labs.map )+
        labs(x = "Absorption fraction of N", y = "Count",
             fill = "Treatment", tag = "D")+
        theme_bw(base_size = 6)
      #p2

      p <- (p1/p2+plot_layout(guides = "collect"))&
        theme(
          plot.tag = element_text(size = 12,face = "plain",vjust  = 0,hjust = 1 ),
          plot.tag.position = "topleft",
          plot.margin = margin(t = 2, r = 0, b = 2, l = 3 ),
          legend.key.size = unit(0.1,"inch"),
          legend.text = element_text(size = 4))
      #open_plot_win(p, 4,3)


      p3 <- ggplot(pave.label.ratio.df%>%filter(c.absorb > 0.02) )+
        #geom_point(aes(x = c.ratio,y = n.ratio, col = group))+
        geom_boxplot(aes(fill = group,x = group , y = c.absorb),col = "black")+
        scale_fill_manual(values =col.map,
                          labels = labs.map )+
        scale_x_discrete(labels  =labs.map)+
        labs(x = NULL, y = "Absorption fraction of C",
             fill = "Treatment", tag = "E")+
        theme_bw(base_size = 6)+
        theme(#axis.text.x = element_text(angle = -45,hjust = 0)
          #axis.text.x = element_blank()
        )
      #p3
      p4 <- ggplot(pave.label.ratio.df%>%filter(n.absorb > 0.02)  )+
        #geom_point(aes(x = c.ratio,y = n.ratio, col = group))+
        geom_boxplot(aes(fill = group,x = group , y = n.absorb),col = "black")+
        scale_fill_manual(values = col.map,
                          labels = labs.map  )+
        scale_x_discrete(labels  =labs.map)+
        labs(x = NULL, y = "Absorption fraction of N",
             fill = "Treatment", tag = "F")+
        theme_bw(base_size = 6)+
        theme(#axis.text.x = element_text(angle = -45,hjust = 0)
          axis.text.x = element_blank()
        )
      #p4

      p <- (p1+p3+p2+p4 + plot_layout(widths = c(0.6,0.4),guides = "collect",axis_titles ="keep") )&
        theme(
          plot.tag = element_text(size = 12,face = "plain",vjust  = 0,hjust = 1 ),
          plot.tag.position = "topleft",
          plot.margin = margin(t = 2, r = 0, b = 0, l = 3 ),
          legend.key.size = unit(0.1,"inch"),
          legend.text = element_text(size = 4))
      p
      open_plot_win(p, 4.5,3)



      p3 <- p3+labs(tag =  "C")+theme(legend.position = "none",
                                      axis.text.x = element_text(angle = -30,hjust = 0.2,vjust = 0.6))

      p4 <- p4+labs(tag =  NULL)+theme(legend.position = "none",
                                       axis.text.x = element_text(angle = -30,hjust = 0.2,vjust = 0.6))
      p <- (p3/p4 + plot_layout(guides = "collect",axis_titles ="keep") )&
        theme(
          plot.tag = element_text(size = 12,face = "plain",vjust  = 0,hjust = 1 ),
          plot.tag.position = "topleft",
          plot.margin = margin(t = 0.5, r = 5, b = 1, l = 3 ),
          legend.key.size = unit(0.1,"inch"),
          legend.text = element_text(size = 4))
      p
      open_plot_win(p, 2,3)

      p <- (p3+p4 + plot_layout(guides = "collect",axis_titles ="keep") )&
        theme(
          plot.tag = element_text(size = 12,face = "plain",vjust  = 0,hjust = 1 ),
          plot.tag.position = "topleft",
          plot.margin = margin(t = 0.5, r = 4, b = 1, l = 1 ),
          #legend.margin =margin(t = 0, r = 0, b = 0 , l = 0 ),
          legend.key.size = unit(0.1,"inch"),
          legend.text = element_text(size = 4))
      p
      open_plot_win(p, 3,1.6)
      #open_plot_ppt(p, 3,1.6)
    }

    pave.label.ratio.df <- pave.label.ratio.df%>%
      dplyr::filter(!is.na(kegg.id),
                    peakMaxo > 1e4 ) %>%
      dplyr::group_by(sample,kegg.id  ) %>%
      dplyr::slice_max(peakMaxo)

    ### Pathway
    {


      pave.label.ratio.stat <- pave.label.ratio.df%>%
        dplyr::group_by(fid )%>%
        dplyr::mutate(
          c.absorb.base = median(c.absorb[group == "NSA"]),
          n.absorb.base = median(n.absorb[group == "NSA"]),
          c.absorb.change = c.absorb - c.absorb.base,
          n.absorb.change = n.absorb - n.absorb.base
        )%>%
        dplyr::group_by(fid,group)%>%
        dplyr::mutate(#c.absorb.change = log2(mean(c.absorb.change)),
          #n.absorb.change = log2(mean(n.absorb.change))
        )



      ### C
      {
        path.res.c <- list()
        for (i.group in unique(pave.label.ratio.stat$group)) {

          up <- pave.label.ratio.stat%>%
            dplyr::filter(group %in% i.group)
          up.path <- analyzePathwayHyperTest(up$kegg.id)%>%
            dplyr::mutate(change = "up")



          up.path$fc.mean <-  sapply(up.path$compounds,function(x){
            kid <- stringr::str_split(x , ";")[[1]]%>%na.omit()
            up$c.absorb.change[up$kegg.id %in%kid ]%>%mean(na.rm = T)
          })



          path.res.c[[i.group]] <- up.path




        }


        }

      ### N
      {
        path.res.n <- list()
        for (i.group in unique(pave.label.ratio.stat$group)) {

          up <- pave.label.ratio.stat%>%
            dplyr::filter(group %in% i.group)
          up.path <- analyzePathwayHyperTest(up$kegg.id)%>%
            dplyr::mutate(change = "up")


          up.path$fc.mean <-  sapply(up.path$compounds,function(x){
            kid <- stringr::str_split(x , ";")[[1]]%>%na.omit()
            up$n.absorb.change[up$kegg.id %in%kid ]%>%mean(na.rm = T)
          })


          path.res.n[[i.group]] <- up.path




        }



      }



      ### DATA MERGE
      {

        min.metabolite <- 5

        kegg_to_short <- c(
          "hsa00470" = "D-Amino acid metabolism",
          "hsa00270" = "Cys & Met metabolism",
          "hsa00350" = "Tyr metabolism",
          "hsa00920" = "Sulfur metabolism",
          "hsa00340" = "His metabolism",
          "hsa00310" = "Lys degradation",
          "hsa00330" = "Arg & Pro metabolism",
          "hsa00260" = "Gly, Ser & Thr metabolism",
          "hsa00380" = "Trp metabolism",
          "hsa00230" = "Purine metabolism",
          "hsa00400" = "Phe, Tyr & Trp biosynthesis",
          "hsa00410" = "Ala metabolism"
        )

        path.res.df.c <- rbindlist(path.res.c,idcol = "group")%>%
          dplyr::filter(group != "NSA",
                        change =="up",
                        grepl("Metabolism", pathway.class))%>%
          dplyr::mutate(
            pathway.abbr = kegg_to_short[pathway.id],
            lp = -log10(p.value),
            Absorption = "C",
            Absorption = ifelse(change =="up", Absorption,NA),
            Absorption = factor(Absorption,levels = c("C","N")))%>%
          dplyr::group_by(pathway.id,group) %>%
          #dplyr::slice_max(lp,n=1,with_ties = F)%>%
          #dplyr::filter(lp != 0)%>%
          dplyr::group_by(pathway.id) %>%
          dplyr::mutate(lp.ave = mean(lp),
                        fc.sd = sd(fc.mean),
                        fc.max = max(fc.mean),
                        hit.ave = mean(Hit))%>%
          dplyr::ungroup()%>%
          dplyr::filter( hit.ave > min.metabolite)%>%
          dplyr::arrange(desc(fc.max))%>%
          dplyr::slice_head(n = 70)%>%
          dplyr::arrange(fc.max)%>%
          dplyr::mutate( pathway.abbr = factor(pathway.abbr,levels = unique(pathway.abbr)) )

        path.res.df.n <- rbindlist(path.res.n,idcol = "group")%>%
          dplyr::filter(group != "NSA",
                        change =="up",
                        grepl("Metabolism", pathway.class))%>%
          dplyr::mutate(lp = -log10(p.value),
                        pathway.abbr = kegg_to_short[pathway.id],
                        #fc.mean = ifelse(),
                        Absorption = "N")%>%
          dplyr::group_by(pathway.id,group) %>%
          #dplyr::slice_max(lp,n=1,with_ties = F)%>%
          #dplyr::filter(lp != 0)%>%
          dplyr::group_by(pathway.id) %>%
          dplyr::mutate(lp.ave = mean(lp),
                        # fc.mean.mean = mean(fc.mean),
                        fc.sd = sd(fc.mean),
                        #  fc.rsd = fc.sd/fc.mean.mean,
                        fc.max = max(fc.mean),
                        hit.ave = mean(Hit))%>%
          dplyr::ungroup()%>%
          dplyr::filter( hit.ave >= min.metabolite)%>%
          dplyr::arrange(desc(fc.max))%>%
          dplyr::slice_head(n = 70)%>%
          dplyr::arrange(hit.ave)%>%
          dplyr::arrange(fc.max)%>%
          dplyr::mutate( pathway.abbr = factor(pathway.abbr,levels = unique(pathway.abbr)) )

      }


      ### plot
      {

        size.max <- max(path.res.df.c$fc.mean, path.res.df.n$fc.mean,na.rm = T) %>% ceiling()
        size.min <- min(path.res.df.c$fc.mean, path.res.df.n$fc.mean,na.rm = T) %>% floor()
        min.to.show <- 0.1
        p1 <- ggplot(path.res.df.c)+
          geom_point(data = data.frame(group = path.res.df.c$group[1],
                                       pathway.abbr = path.res.df.c$pathway.abbr[1],
                                       lp = c(0.01,0.01),Absorption = c("C","N")),
                     aes(x = group,y = pathway.abbr,
                         size = lp, fill = Absorption),pch = 21)+
          geom_point(aes(x = group,y = pathway.abbr,
                         size = fc.mean, fill = Absorption),pch = 21)+
          #scale_fill_manual(values = c("up" = "#E33F31","down" = "white"),guide = "none")+
          scale_fill_manual(values = c("C" = "#FF7F0E","N" = "#1F77B4"),
                            drop = F,na.value = "white",na.translate = FALSE)+
          scale_radius(breaks = ( size.max * seq(0,1,0.2)),
                       range = c(0.2, 3),limits = c(min.to.show, size.max))+
          scale_x_discrete(labels  = c("+mixture","+leucine","+threonine",
                                       "+tryptophan","+adenine","+uracil","+acetate"))+
          theme_bw(base_size = 6)+
          labs(x = NULL, y = NULL,size = "Abosorption\nIncrease",fill = "Abosorption",tag = NULL)+
          theme(legend.position = "right",
                legend.title.position = "top",
                legend.title = element_text(margin = margin(b = 2)),
                legend.box.spacing = unit(0.01, "inch"),
                legend.key.height = unit(0.1, "inch"),
                axis.text.y = element_text(size = 6),
                axis.text.x = element_text(size = 4,
                                           angle = -30,hjust = 0)
          )

        p2 <- ggplot(path.res.df.n)+
          geom_point(aes(x = group,y = pathway.abbr,
                         size = fc.mean, fill = Absorption),pch = 21,show.legend = F)+
          #scale_fill_manual(values = c("up" = "#E33F31","down" = "white"),guide = "none")+
          scale_fill_manual(values = c("C" = "#FF7F0E","N" = "#1F77B4"),
                            drop = F,na.value = "white",na.translate = FALSE)+
          scale_radius(breaks = size.max * seq(0,1,0.2),
                       range = c(0.2, 3),limits = c(min.to.show, size.max))+
          scale_x_discrete(labels  = c("+mixture","+leucine","+threonine",
                                       "+tryptophan","+adenine","+uracil","+acetate"))+
          theme_bw(base_size = 6)+
          labs(x = NULL, y = NULL,size = "Abosorption\nIncrease",fill = "Abosorption")+
          theme(legend.position = "right",
                legend.title.position = "top",
                legend.title = element_text(margin = margin(b = 2)),
                legend.box.spacing = unit(0.01, "inch"),
                legend.key.height = unit(0.1, "inch"),,
                axis.text.y = element_text(size = 6),
                axis.text.x = element_text(size= 4,
                                           angle = -30,hjust = 0))

        p1 <- p1+labs(tag =  "D")
        p <- p1/p2+
          plot_layout(guides = "collect",axes = "collect")&
          theme(
            plot.tag = element_text(size = 12,face = "plain",vjust  = 0),
            plot.tag.position = "topleft",
            plot.margin = margin(t = 1, r =0, b = 2, l = 0 ))
        open_plot_win(p,2.7,3.2)
        open_plot_ppt(p,2.7,3.2)

      }




      }


    ### show selected pathway
    {

      # edit_df_in_excel(path.cp <- path.res.c$NSE)
      path <- path.res.c$NSE%>%
        #dplyr::filter(pathway.id == "hsa00920")%>% ### SULFUR
        # dplyr::filter(pathway.id == "hsa00230")%>% ### purine
        # dplyr::filter(pathway.id == "hsa00380")%>% ### TRP
        dplyr::filter(pathway.id == "hsa00350")%>% ### tyr
        # dplyr::filter(pathway.id == "hsa00480")%>% ### GSH
        # dplyr::filter(pathway.id == "hsa00400")%>%#Phenylalanine, tyrosine and tryptophan biosynthesis
        #  dplyr::filter(pathway.id == "hsa00410")%>% ### beta-Alanine metabolism
        # dplyr::filter(pathway.id == "hsa00310")%>% ### Lysine degradation
        #  dplyr::filter(pathway.id == "hsa00770")%>% ### Pantothenate and CoA biosynthesis
        #    dplyr::filter(pathway.id == "hsa00760")%>% ### Nicotinate and nicotinamide metabolism
        dplyr::mutate()

      path.cp <- path%>%
        dplyr::pull(compounds)%>%
        stringr::str_split( ";")%>%
        `[[`(1)

      cp.name.map <- c("C04751" = "CAIR")

      path.df <- pave.label.ratio.stat%>%
        dplyr::filter( kegg.id %in% path.cp)%>%
        pivot_longer(c("c.absorb.change","n.absorb.change"),
                     names_to = "abs",
                     values_to = "change")%>%
        dplyr::mutate(abs = ifelse(abs == "c.absorb.change","C","N"),
                      name = case_when(
                        kegg.id %in% names(cp.name.map)~cp.name.map[kegg.id],
                        T~name
                      ))%>%
        dplyr::group_by(group,abs,kegg.id)%>%
        dplyr::mutate(change = mean(change,na.rm = T))%>%
        dplyr::ungroup()%>%
        dplyr::filter(group!= "NSA",
                      ! kegg.id %in%c("C00955","C00130") )

      path.df%>%
        dplyr::distinct(kegg.id,.keep_all = T)%>%
        dplyr::pull(name,name = kegg.id)

      mix <- path.df %>%
        dplyr::group_by(name,group)%>%
        dplyr::mutate(d = max(change)-min(change),
                      n = n())%>%
        dplyr::filter(d < 0.1)

      p <- ggplot(path.df)+
        geom_point(#data = subset(path.df, abs == "N"),
          aes(x = name ,
              y = group,
              col = abs,
              size = change),stroke  = 0.8,alpha = 0.5 ,pch = 1,
          #position = position_nudge(y = 0.2),
          fill = "transparent")+
        geom_point(data = mix,
                   aes(x = name ,
                       y = group,
                       size = change),stroke  = 0.8,alpha = 0.5 ,pch = 1,
                   col = "#8D5397",
                   #position = position_nudge(y = -0.2),
                   fill = "transparent")+
        scale_shape_manual(values = c("C"= 1  ,"N"=1))+
        scale_y_discrete(labels  = c("+AA/uracil\n/adenine","+leucine","+threonine",
                                     "+tryptophan","+adenine","+uracil","+acetate"))+
        scale_radius(breaks = seq(0,1,0.2),
                     range = c(0.2, 3),limits = c(0.05, 1.2))+
        scale_color_manual(values = c("C" = "#FF7F0E","N" = "#1F77B4"),
                           drop = F,na.value = "white",na.translate = FALSE)+
        labs(title = path$pathway.name,
             x = NULL, y = NULL,size = "Abosorption\nIncrease",col = "Abosorption")+
        coord_flip()+
        theme_bw(base_size = 6)+
        theme(
          legend.position = "right",
          legend.title.position = "top",
          legend.title = element_text(margin = margin(b = 2)),
          legend.box.spacing = unit(0.01, "inch"),
          legend.key.height = unit(0.1, "inch"),
          #axis.text.y = element_text(size = 4),
          axis.text.x = element_text(angle = -30,hjust = 0)
        )
      open_plot_win(p,2.6,2)
      open_plot_ppt(p,2.6,1.8)



      cpd <-path.df%>%
        dplyr::filter(group=="NSB")%>%
        dplyr::distinct(kegg.id,abs,.keep_all = T)%>%
        pivot_wider(names_from   = "abs",values_from ="change" )%>%
        column_to_rownames("kegg.id")%>%
        dplyr::select(C,N)%>%
        as.matrix()

      plot_cpd_pathview(cpd =  cpd,pathway.id = "hsa00760",
                        split.group = T,multi.state = TRUE,
                        limit = list(
                          cpd = c(0, 1), gene = 1),
                        low = list(gene = "grey", cpd =
                                     "grey"),
                        mid = list(gene = "#FF964E", cpd =
                                     "#FF964E"),
                        dir.to.save = "d:/temp/")



      ### heatmap pos
      {


        ### QA
        {
          i.mz <- 168.0288
          i.pave.formula <- "C7N1"
          i.rt <- 955.8462
        }

        ## NAD
        {
          i.mz <- 664.1148
          i.pave.formula <- "C21N7"
          i.rt <- 872.6339
        }


        ## Nicotinamide
        {
          i.mz <- 123.055290
          i.pave.formula <- "C6N2"
          i.rt <- 874.8
        }

        ## GSH
        {
          i.mz <- 308.09054
          i.pave.formula <- "C10N3"
          i.rt <- 875.4
        }

        ## ATP
        {
          i.mz <- 508.0023445
          i.pave.formula <- "C10N5"
          i.rt <- 934.0222374

        }

        ## QA
        {
          i.mz <- 168.0288
          i.pave.formula <- "C7N1"
          i.rt <- 955.8462

        }

        ## GLU
        {
          i.mz <- 148.0601
          i.pave.formula <- "C5N1"
          i.rt <- 911.2694

        }

        ## GLYCINE
        {
          i.mz <- 76.03916
          i.pave.formula <- "C2N1"
          i.rt <- 843.8066

        }

        ## CYS
        {
          i.mz <- 122.027031
          i.pave.formula <- "C3N1"
          i.rt <- 871.8

        }

        ## glutamine
        {
          i.mz <- 147.075903566073
          i.pave.formula <- "C5N2"
          i.rt <- 838.2881214


        }


        ### heatmap
        {

          xcms.mat<- xcms.mat.pos
          xcms.se <- xcms.se.pos

          i.cn.count <- MSCC::chemform_parse(i.pave.formula)
          i.cn.diff <- get_CN_mass_diff_table(i.cn.count[1,"C"],i.cn.count[1,"N"])
          i.cn.table <- i.cn.diff[,mz:=mass_diff + i.mz]
          i.cn.match <- match_mz_foverlaps(i.cn.diff$mz,xcms.fdf.pos$mzmed)
          i.cn.match[,fid := ion2] [,rt := xcms.fdf.pos$rtmed[fid]]
          i.cn.table <- cbind(i.cn.match,i.cn.table[i.cn.match$ion1]  )
          i.cn.table[,label_pattern := paste0("C",C_count ,"N",N_count )][
            !is.na(fid)]
          i.cn.table <- i.cn.table[abs(rt -i.rt) < 5]%>%
            dplyr::group_by(label_pattern)%>%
            dplyr::slice_min(mz.ppm,n = 1)%>%
            dplyr::ungroup()%>%
            dplyr::arrange(N_count,C_count)

          i.cn.exp.mat <- xcms.mat[i.cn.table$fid,xcms.se$sample.name]
          rownames(i.cn.exp.mat) <-i.cn.table$label_pattern
          norm.to <- i.cn.exp.mat[,xcms.se$group == "NSA"]%>%sum(na.rm = T)/10
          i.cn.exp.mat <- (i.cn.exp.mat/norm.to)

          ### C/N label ratio
          if(T){
            i.cn.label.ratio <- i.cn.exp.mat%>%
              as.data.frame()%>%
              rownames_to_column("label_pattern")%>%
              pivot_longer(-"label_pattern" ,names_to = "sample")%>%
              left_join(i.cn.table,by = "label_pattern")%>%
              dplyr::group_by(sample)%>%
              dplyr::mutate(
                group = setNames(nm = xcms.se$sample.name,
                                 object = xcms.se$sample.type)[sample],
                int.sum = sum(value),
                c.cum = sum(C_count * value),
                c.ratio = c.cum/int.sum,
                n.cum = sum(N_count * value),
                n.ratio = n.cum/int.sum
              )%>%
              #dplyr::filter(grepl("NS",group))%>%
              #dplyr::filter(!grepl("Blank",group))%>%
              dplyr::ungroup()%>%
              dplyr::mutate(x = c.ratio + rnorm(n(),sd = 0.0001),
                            y = n.ratio + rnorm(n(),sd = 0.0001))

            res <- i.cn.label.ratio%>%
              dplyr::filter(mz == max(mz))%>%
              dplyr::mutate(fid = i.fid,
                            mz = i.mz)



          }

          #if(nrow(i.cn.exp.mat) > 4 ){
          if(T){

            n.ratio <- res %>%
              dplyr::pull(n.ratio,name = sample)
            c.ratio <- res %>%
              dplyr::pull(c.ratio,name = sample)
            int.sum <- res %>%
              dplyr::pull(int.sum,name = sample)


            #hm.mat <- log2(i.cn.exp.mat)
            hm.mat <- (i.cn.exp.mat)
            hm.mat <- apply(hm.mat,2, function(x){ x/sum(x) })
            hm.mat <- hm.mat[(rowSums(hm.mat) > 0.1),  ]
            hm <- Heatmap(hm.mat,
                          #col = colramp(breaks = c(min(min(hm.mat),-1),0,
                          #                         max((hm.mat))),colors = c("#259644","white","#D84704")),
                          col = colramp(breaks = c(0,max(hm.mat)/2,max(hm.mat))),
                          #col = colramp(breaks = c(0,0.5 ,1)),
                          row_names_gp = gpar(fontsize = 6),
                          heatmap_legend_param = list(title =paste0("Abundance"),
                                                      grid_width  = unit(0.1, "inch"),
                                                      title_gp  =gpar(fontsize = 6),
                                                      labels_gp = gpar(fontsize = 6)),

                          column_title = c("Control","+mixture","+leucine","+threonine",
                                           "+tryptophan","+adenine","+uracil","+acetate"
                                           #"12C14N","13C14N","12C15N","13C15N"
                          ),
                          column_title_side = "bottom",
                          column_title_rot = -30,
                          column_title_gp = gpar(fontsize = 6),

                          show_heatmap_legend = F,

                          cluster_rows = F,cluster_columns = F,cluster_column_slices = F,
                          column_split = xcms.se$group,
                          show_column_names = F,
                          row_names_side = "left",
                          rect_gp = gpar(col = "black"),
                          top_annotation = columnAnnotation(
                            #intsum = anno_barplot(int.sum,gp = gpar(fill = "#DF3A2D"),
                            #                      axis_param = list(at = ceiling(max(int.sum)) * c(0.5,1),
                            #                                        gp = gpar(fontsize = 6),
                            #                                        side = "right")),
                            cratio = anno_barplot(c.ratio,gp = gpar(fill = "#FF7F0E"),
                                                  axis_param = list(
                                                    #at = ceiling(max(c.ratio)) * c(0.5,1),
                                                    at = 0,
                                                    gp = gpar(fontsize = 0),
                                                    side = "right")),
                            nratio = anno_barplot(n.ratio,gp = gpar(fill = "#1F77B4"),
                                                  axis_param = list(#at = ceiling(max(n.ratio)) * c(0.5,1),
                                                    at = 0,
                                                    gp = gpar(fontsize = 0),
                                                    side = "right")),
                            annotation_label = c("Labeled\nFraction of C","Labeled\nFraction of N"),
                            annotation_name_gp  = gpar(fontsize = 6),
                            annotation_name_side  = "left",
                            annotation_name_rot  = 0,
                            height = unit(0.4,"inch")
                          )

            )
            hm
            #open_plot_win(hm, 2.5 , 1.12)

            open_plot_win(hm, 3,1.3 + nrow(hm.mat) * 0.06)
          }
          #readline()

        }

        open_plot_ppt(hm, 2.5 , 1.5)


      }


      ### heatmap neg
      {


        ### R5P
        {
          i.mz <- 229.011887
          i.pave.formula <- "C5N0"
          i.rt <- 910.2
        }




        ### heatmap
        {

          xcms.mat<- xcms.mat.neg
          xcms.se <- xcms.se.neg

          i.cn.count <- MSCC::chemform_parse(i.pave.formula)
          i.cn.diff <- get_CN_mass_diff_table(i.cn.count[1,"C"],i.cn.count[1,"N"])
          i.cn.table <- i.cn.diff[,mz:=mass_diff + i.mz]
          i.cn.match <- match_mz_foverlaps(i.cn.diff$mz,xcms.fdf.neg$mzmed)
          i.cn.match[,fid := ion2] [,rt := xcms.fdf.neg$rtmed[fid]]
          i.cn.table <- cbind(i.cn.match,i.cn.table[i.cn.match$ion1]  )
          i.cn.table[,label_pattern := paste0("C",C_count ,"N",N_count )][
            !is.na(fid)]
          i.cn.table <- i.cn.table[abs(rt -i.rt) < 5]%>%
            dplyr::group_by(label_pattern)%>%
            dplyr::slice_min(mz.ppm,n = 1)%>%
            dplyr::ungroup()%>%
            dplyr::arrange(N_count,C_count)

          i.cn.exp.mat <- xcms.mat[i.cn.table$fid,xcms.se$sample.name,drop = F]
          rownames(i.cn.exp.mat) <-i.cn.table$label_pattern
          norm.to <- i.cn.exp.mat[,xcms.se$group == "NSA"]%>%sum(na.rm = T)/10
          i.cn.exp.mat <- (i.cn.exp.mat/norm.to)

          ### C/N label ratio
          if(T){
            i.cn.label.ratio <- i.cn.exp.mat%>%
              as.data.frame()%>%
              rownames_to_column("label_pattern")%>%
              pivot_longer(-"label_pattern" ,names_to = "sample")%>%
              left_join(i.cn.table,by = "label_pattern")%>%
              dplyr::group_by(sample)%>%
              dplyr::mutate(
                group = setNames(nm = xcms.se$sample.name,
                                 object = xcms.se$sample.type)[sample],
                int.sum = sum(value),
                c.cum = sum(C_count * value),
                c.ratio = c.cum/int.sum,
                n.cum = sum(N_count * value),
                n.ratio = n.cum/int.sum
              )%>%
              #dplyr::filter(grepl("NS",group))%>%
              #dplyr::filter(!grepl("Blank",group))%>%
              dplyr::ungroup()%>%
              dplyr::mutate(x = c.ratio + rnorm(n(),sd = 0.0001),
                            y = n.ratio + rnorm(n(),sd = 0.0001))

            res <- i.cn.label.ratio%>%
              dplyr::filter(mz == max(mz))%>%
              dplyr::mutate(fid = i.fid,
                            mz = i.mz)



          }

          #if(nrow(i.cn.exp.mat) > 4 ){
          if(T){

            n.ratio <- res %>%
              dplyr::pull(n.ratio,name = sample)
            c.ratio <- res %>%
              dplyr::pull(c.ratio,name = sample)
            int.sum <- res %>%
              dplyr::pull(int.sum,name = sample)


            #hm.mat <- log2(i.cn.exp.mat)
            hm.mat <- (i.cn.exp.mat)
            hm.mat <- apply(hm.mat,2, function(x){ x/sum(x) })
            hm.mat <- hm.mat[(rowSums(hm.mat) > 0.001),  ]
            hm <- Heatmap(hm.mat,
                          #col = colramp(breaks = c(min(min(hm.mat),-1),0,
                          #                         max((hm.mat))),colors = c("#259644","white","#D84704")),
                          col = colramp(breaks = c(0,max(hm.mat)/2,max(hm.mat))),
                          row_names_gp = gpar(fontsize = 6),
                          heatmap_legend_param = list(title =paste0("Abundance"),
                                                      grid_width  = unit(0.1, "inch"),
                                                      title_gp  =gpar(fontsize = 6),
                                                      labels_gp = gpar(fontsize = 6)),

                          column_title = c("Control","+AA/uracil/adenine","+leucine","+threonine",
                                           "+tryptophan","+adenine","+uracil","+acetate"
                                           #"12C14N","13C14N","12C15N","13C15N"
                          ),
                          column_title_side = "bottom",
                          column_title_rot = -30,
                          column_title_gp = gpar(fontsize = 6),

                          show_heatmap_legend = F,

                          cluster_rows = F,cluster_columns = F,cluster_column_slices = F,
                          column_split = xcms.se$group,
                          show_column_names = F,
                          row_names_side = "left",
                          rect_gp = gpar(col = "black"),
                          top_annotation = columnAnnotation(
                            #intsum = anno_barplot(int.sum,gp = gpar(fill = "#DF3A2D"),
                            #                      axis_param = list(at = ceiling(max(int.sum)) * c(0.5,1),
                            #                                        gp = gpar(fontsize = 6),
                            #                                        side = "right")),
                            cratio = anno_barplot(c.ratio,gp = gpar(fill = "#FF7F0E"),
                                                  axis_param = list(
                                                    #at = ceiling(max(c.ratio)) * c(0.5,1),
                                                    at = 0,
                                                    gp = gpar(fontsize = 0),
                                                    side = "right")),
                            nratio = anno_barplot(n.ratio,gp = gpar(fill = "#1F77B4"),
                                                  axis_param = list(#at = ceiling(max(n.ratio)) * c(0.5,1),
                                                    at = 0,
                                                    gp = gpar(fontsize = 0),
                                                    side = "right")),
                            annotation_label = c("Labeled\nFraction of C","Labeled\nFraction of N"),
                            annotation_name_gp  = gpar(fontsize = 6),
                            annotation_name_side  = "left",
                            annotation_name_rot  = 0,
                            height = unit(0.5,"inch")
                          )

            )
            hm
            open_plot_win(hm, 3,1.3 + nrow(hm.mat) * 0.1)


          }
          #readline()

        }

      }

    }







  }


  {
    p <- ggplot(absorb.top20, aes(x = nutrient, y = sum_absorb_ratio, fill = absorb_display)) +
      geom_col(color = "grey30", linewidth = 0.15, width = 0.75) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
      scale_fill_manual(
        values = c(
          setNames(
            viridisLite::viridis(length(top20.absorb$absorb_pattern), option = "turbo", direction = -1),
            top20.absorb$absorb_pattern
          ),
          Other = "#CCCCCC"
        )
      ) +
      labs(
        #title = "Top 20 absorb forms + Other across nutrient conditions",
        #subtitle = "Equal bar height per condition; lower top-20 slice in mixture reflects broader absorb-form usage",
        x = NULL,
        y = "Absorb form frequency",
        fill = "Absorb form"
      ) +
      theme_bw(base_size = 6) +
      theme(
        axis.text.x = element_text(angle = -35, hjust = 0),
        legend.key.size = unit(0.1, "cm"),
        legend.text = element_text(size = 4)
      )
    open_plot_win(p,3,2)


    ggplot(absorb.top20.hm, aes(x = nutrient, y = absorb_pattern, fill = relative_absorption)) +
      geom_tile(color = "white", linewidth = 0.3) +
      scale_fill_gradient2(
        low = "#FFFFFF",
        mid = "#D13B2E",
        high = "#D13B2E",
        midpoint = 0.3,
        limits = c(0, 0.4),
        name = "Relative\nabsorption"
      ) +
      labs(
        #title = "Relative absorption of top 20 absorb forms",
        #subtitle = "Row-normalized sum absorb ratio per nutrient condition",
        x = NULL,
        y = "Absorb form"
      ) +
      theme_bw(base_size = 6) +
      theme(
        axis.text.x = element_text(angle = -35, hjust = 0),
        legend.key.size = unit(0.2, "cm"),
        panel.grid = element_blank()
      )->p

    open_plot_win(p,3,2.5)

    col_fun.top20 <- circlize::colorRamp2(
      c(0, 0.15, 0.2),
      c("#FFFFFF", "#FFCD99", "#D13B2E")
    )

    ht.top20 <- Heatmap(
      hm_mat.top20,
      name = "Relative\nabsorption",
      col = col_fun.top20,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      #column_title = "Relative absorption of top 20 absorb forms",
      column_title_gp = grid::gpar(fontsize = 6, fontface = "bold"),
      column_labels = nutrient.labels[nutrient.groups],
      column_names_rot = -30,
      row_names_gp = grid::gpar(fontsize = 6),
      column_names_gp = grid::gpar(fontsize = 6),
      rect_gp = grid::gpar(fill = "white",col = "black"),
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.circle(
          x = x,
          y = y,
          r = min(unit.c(width, height)) * 0.42,
          gp = grid::gpar(fill = fill, col = "black")
        )
      },
      heatmap_legend_param = list(
        title_gp = grid::gpar(fontsize = 6),
        labels_gp = grid::gpar(fontsize = 6)
      )
    )

    open_plot_win(ht.top20,3,2.5)


    p_scatter <- ggplot(global.absorb, aes(x = c_absorb_mean, y = n_absorb_mean, color = group, fill = group)) +
      geom_point(size = 0.2, alpha = 1) +
      #stat_ellipse(geom = "polygon", alpha = 0.15, linewidth = 0.4) +
      scale_color_manual(values = nutrient.colors, labels = nutrient.labels) +
      scale_fill_manual(values = nutrient.colors, labels = nutrient.labels) +
      labs(
        #title = "Absorb-form preference across metabolites",
       # subtitle = "Mean absorbed C and N counts (C0N0 excluded)",
        x = "Mean absorbed C count",
        y = "Mean absorbed N count",
        color = "Treatment",
        fill = "Treatment"
      ) +
      theme_bw(base_size = 6) +
      theme(legend.position = "none")
    open_plot_win(p_scatter,2,1.5)

  }
}

{

  # Simulate 3-edge paths built from adduct, fragment, and isotope transitions.
  # Columns: path_id, edge1, edge2, edge3, net_form, net_form_adduct, net_mass_shift.
  # Generate 20 paths covering mixed edge-type combinations and net formula changes.

  iso.mass.diff <- get_iso_mass_diff()
  demo_adduct_names <- c(
    "[M-H]-", "[M+Na-2H]-", "[M+K-2H]-",
    "[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+"
  )
  ad.mass.diff <- data.table::rbindlist(list(
    get_adduct_mass_diff(polarity = 0),
    get_adduct_mass_diff(polarity = 1)
  ))[
    mass_diff != 0 &
      adduct.from %in% demo_adduct_names &
      adduct.to %in% demo_adduct_names
  ]
  fg.mass.diff <- get_fragment_mass_diff()

  # Adduct edges chain: a later adduct must start from the last adduct endpoint
  # (even if fragment/isotope edges lie in between).
  .pick_adduct_edge <- function(idx, from_adduct = NULL) {
    pool <- ad.mass.diff
    if (!is.null(from_adduct)) {
      pool <- pool[adduct.from == from_adduct]
      if (!nrow(pool)) {
        stop("No adduct transition from ", from_adduct, call. = FALSE)
      }
    }
    idx <- ((idx - 1L) %% nrow(pool)) + 1L
    row <- pool[idx]
    list(
      type = "adduct",
      adduct.from = row$adduct.from,
      adduct.to = row$adduct.to,
      label = paste0("adduct: ", row$adduct.from, " > ", row$adduct.to),
      chemform_diff = row$chemform_diff,
      mass_diff = row$mass_diff
    )
  }

  .pick_transition_edge <- function(type, idx, adduct_from = NULL) {
    if (type == "adduct") {
      return(.pick_adduct_edge(idx, from_adduct = adduct_from))
    }
    dt <- switch(
      type,
      isotope = iso.mass.diff,
      fragment = fg.mass.diff,
      stop("Unknown transition type: ", type)
    )
    idx <- ((idx - 1L) %% nrow(dt)) + 1L
    row <- dt[idx]
    detail <- switch(
      type,
      isotope = row$element,
      fragment = chemform_simplify(row$chemform_diff)
    )
    list(
      type = type,
      label = paste0(type, ": ", detail),
      chemform_diff = row$chemform_diff,
      mass_diff = row$mass_diff
    )
  }

  .pick_path_edges <- function(spec) {
    adduct_cursor <- NULL
    edges <- vector("list", 3L)
    types <- c(spec$type1, spec$type2, spec$type3)
    idxs <- c(spec$idx1, spec$idx2, spec$idx3)
    for (k in seq_along(types)) {
      edges[[k]] <- .pick_transition_edge(
        types[k],
        idxs[k],
        adduct_from = if (types[k] == "adduct") adduct_cursor else NULL
      )
      if (types[k] == "adduct") {
        adduct_cursor <- edges[[k]]$adduct.to
      }
    }
    edges
  }

  .path_net <- function(edges) {
    net_form <- chemform_simplify(MSCC:::chemform_sum(vapply(edges, `[[`, "", "chemform_diff")))
    adduct_forms <- vapply(edges[vapply(edges, function(e) e$type == "adduct", logical(1))],
                           `[[`, "", "chemform_diff")
    net_form_adduct <- if (length(adduct_forms)) {
      chemform_simplify(MSCC:::chemform_sum(adduct_forms))
    } else {
      ""
    }
    net_mass_shift <- sum(vapply(edges, `[[`, numeric(1), "mass_diff"))
    list(
      net_form = net_form,
      net_form_adduct = net_form_adduct,
      net_mass_shift = net_mass_shift
    )
  }

  path_specs <- data.frame(
    path_id = 1:20,
    type1 = c(
      "adduct", "adduct", "adduct", "fragment", "fragment", "fragment",
      "isotope", "isotope", "isotope", "adduct", "adduct", "fragment",
      "fragment", "isotope", "isotope", "adduct", "fragment", "isotope",
      "adduct", "fragment"
    ),
    idx1 = c(1, 2, 3, 1, 3, 8, 1, 2, 5, 1, 4, 1, 4, 1, 3, 2, 2, 4, 5, 6),
    type2 = c(
      "fragment", "isotope", "adduct", "fragment", "isotope", "adduct",
      "fragment", "adduct", "isotope", "isotope", "fragment", "adduct",
      "isotope", "fragment", "adduct", "fragment", "isotope", "adduct",
      "isotope", "adduct"
    ),
    idx2 = c(2, 3, 1, 5, 1, 1, 3, 3, 2, 4, 6, 2, 7, 2, 5, 9, 1, 6, 3, 4),
    type3 = c(
      "isotope", "fragment", "fragment", "adduct", "adduct", "isotope",
      "adduct", "fragment", "fragment", "adduct", "isotope", "isotope",
      "adduct", "adduct", "fragment", "isotope", "adduct", "fragment",
      "fragment", "isotope"
    ),
    idx3 = c(3, 1, 5, 1, 2, 4, 2, 5, 6, 1, 2, 4, 1, 3, 8, 1, 7, 3, 5, 2),
    stringsAsFactors = FALSE
  )

  x <- do.call(
    rbind,
    lapply(seq_len(nrow(path_specs)), function(i) {
      spec <- path_specs[i, ]
      path_edges <- .pick_path_edges(spec)
      e1 <- path_edges[[1L]]
      e2 <- path_edges[[2L]]
      e3 <- path_edges[[3L]]
      net <- .path_net(list(e1, e2, e3))
      data.frame(
        path_id = spec$path_id,
        edge1 = e1$label,
        edge2 = e2$label,
        edge3 = e3$label,
        net_form = net$net_form,
        net_form_adduct = net$net_form_adduct,
        net_mass_shift = net$net_mass_shift,
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(x) <- NULL
  x
}

