

# Mon Jun  1 20:01:43 2026 ------------------------------
{


  {
    trace.cor <- 0.5
    trace.demo <- TRACE_get_CN_net(trace.demo,0,TRACE_cor_cutoff = trace.cor,
                                   ratio.adjust = c(1.0842,0.7641,1.0566,1.0707))
    trace.demo <- TRACE_get_CN_net(trace.demo,1,TRACE_cor_cutoff = trace.cor,
                                   ratio.adjust = c(1.0842,0.7641,1.0566,1.0707))
  }

  a <- TRACE_get_CN_labelling_ratio(trace.demo,eval_top = 1,plot = T,)

  trace.demo <- TRACE_CN_labelling_ratio_adjust(trace.demo,eval_top = 0.3,plot = T,reconstruct = T)

  ratio.adj <- c(1.0782,0.6663,1.0357,1.0992)
  names(ratio.adj) <- names(a)[4:7]

  df <- a%>%
    pivot_longer(4:7)%>%
    dplyr::slice_max(TRACE_cor,prop = 1 )%>%
    dplyr::mutate(
      regeion = cut(TRACE_cor, breaks = seq(0, 1, 0.05)),
      ratio.bench = ratio.adj[name],
      ratio.error = value - ratio.bench
    )



  ggplot(df,aes(x = name , y = value, col = TRACE_cor))+
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
    scale_color_gradient(low = "yellow",high = "red")




  ggplot(df,aes(x = TRACE_cor , y = ratio.error,col = TRACE_cor))+
    geom_point(alpha = 0.3)+
    scale_color_gradient(low = "yellow",high = "red")

  library(ggridges)

  p <- ggplot(df,aes(x = ratio.error , y = regeion,fill = regeion))+
    geom_density_ridges(   )+
    scale_fill_manual(values =  MSdev:::colramp()(seq(0.1,1,0.1)))+
    labs(x = "Ratio shift", y = NULL, fill = expression( rho ~ "range"))+
    theme_bw(base_size =  6 )+
    theme(legend.key.size = unit(0.1,"inch") )

  open_plot_win(p,3,3)



}
